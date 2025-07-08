"""
Module for handling GTF annotation files.

Classes
-------
GTF: A class to handle GTF annotation files.
"""

from pathlib import Path

import numpy as np
import pandas as pd
from pyranges import PyRanges as pr
from pyranges import read_gtf
from tqdm import tqdm

from .gene import Gene, GeneSets


class GTF:
    """Base class for handling GTF annotation files"""

    def __init__(self, gtf_path, exclude_chrs=["chrM", "chrY"]):
        """Initialize GTF reader with a path to GTF file

        Parameters
        ----------
        gtf_path: Path to GTF file (can be gzipped)
        exclude_chrs: List of chromosomes to exclude, defaults to ["chrM", "chrY"]
        """
        self.gtf_path = Path(gtf_path)

        if not self.gtf_path.exists():
            raise FileNotFoundError(f"GTF file not found: {self.gtf_path}")

        self.gtf = self._load_gtf(exclude_chrs)

    def _load_gtf(self, exclude_chrs):
        """Load and process GTF file into standardized format.

        We also filter out the excluded chromosomes and remove genes with multiple chromosomes.

        Parameters
        ----------
        exclude_chrs: List of chromosomes to exclude, defaults to ["chrM", "chrY"]

        Returns
        -------
        pd.DataFrame
            The GTF data in a standardized format.
        """
        gtf_df = read_gtf(str(self.gtf_path)).as_df()
        self._original_gtf = gtf_df
        # Process transcripts for positive and negative strands
        positive = gtf_df[(gtf_df.Feature == "transcript") & (gtf_df.Strand == "+")][
            [
                "Chromosome",
                "Start",
                "Start",
                "Strand",
                "gene_name",
                "gene_id",
                "gene_type",
                "transcript_id",
            ]
        ]
        negative = gtf_df[(gtf_df.Feature == "transcript") & (gtf_df.Strand == "-")][
            [
                "Chromosome",
                "End",
                "End",
                "Strand",
                "gene_name",
                "gene_id",
                "gene_type",
                "transcript_id",
            ]
        ]

        columns = [
            "Chromosome",
            "Start",
            "End",
            "Strand",
            "gene_name",
            "gene_id",
            "gene_type",
            "transcript_id",
        ]
        positive.columns = columns
        negative.columns = columns

        gtf_df = pd.concat([positive, negative], axis=0).drop_duplicates().reset_index()
        gtf_df["gene_id"] = gtf_df.gene_id.str.split(".").str[0]

        # Filter excluded chromosomes
        gtf_df = gtf_df[~gtf_df.Chromosome.isin(exclude_chrs)]

        # Remove genes with multiple chromosomes
        gtf_df["chrom_count"] = gtf_df.groupby("gene_name")["Chromosome"].transform(
            "nunique"
        )
        gtf_df = gtf_df[gtf_df.chrom_count == 1]

        return gtf_df

    def get_gene(self, gene_name) -> Gene:
        """Get a Gene object for the given gene name.

        Parameters
        ----------
        gene_name: str
            The gene name to get the Gene object for.

        Returns
        -------
        Gene
            The Gene object.
        """
        df = self.gtf[self.gtf.gene_name == gene_name]
        return Gene(
            name=df.gene_name.iloc[0],
            id=df.gene_id.iloc[0],
            chrom=df.Chromosome.iloc[0],
            strand=df.Strand.iloc[0],
            tss_list=df[
                [
                    "Chromosome",
                    "Start",
                    "End",
                    "Strand",
                    "gene_name",
                    "gene_id",
                    "gene_type",
                ]
            ],
        )

    def get_genes(self, gene_names) -> GeneSets:
        """Get a list of Gene objects for the given gene names.

        Parameters
        ----------
        gene_names: List of gene names

        Returns
        -------
        GeneSets
            A list of Gene objects.
        """
        gene_names = np.intersect1d(gene_names, np.unique(self.gtf.gene_name.values))
        return GeneSets([self.get_gene(gene_name) for gene_name in tqdm(gene_names)])

    def get_gene_id(self, gene_id) -> Gene:
        """Get a Gene object for the given gene ID.

        Parameters
        ----------
        gene_id: str
            The gene ID to get the Gene object for.

        Returns
        -------
        Gene
            The Gene object.
        """
        df = self.gtf[self.gtf.gene_id.str.startswith(gene_id)]
        return Gene(
            name=df.gene_name.iloc[0],
            id=df.gene_id.iloc[0],
            chrom=df.Chromosome.iloc[0],
            strand=df.Strand.iloc[0],
            tss_list=df[
                [
                    "Chromosome",
                    "Start",
                    "End",
                    "Strand",
                    "gene_name",
                    "gene_id",
                    "gene_type",
                ]
            ],
        )

    def get_genebodies(self, gene_names=None) -> pd.DataFrame:
        """Get the gene bodies for the given gene names.

        Parameters
        ----------
        gene_names: List of gene names, optional
            The gene names to get the gene bodies for, defaults to None.

        Returns
        -------
        pd.DataFrame
            The gene bodies.
        """
        genebodies = self.gtf.query('gene_type == "protein_coding"')
        genebodies["Start"] = genebodies.groupby(["Chromosome", "gene_name"])[
            "Start"
        ].transform("min")
        genebodies["End"] = genebodies.groupby(["Chromosome", "gene_name"])[
            "End"
        ].transform("max")
        genebodies = genebodies.drop_duplicates(subset=["Chromosome", "gene_name"])
        genebodies = genebodies[
            [
                "Chromosome",
                "Start",
                "End",
                "Strand",
                "gene_name",
                "gene_id",
                "gene_type",
            ]
        ]
        if gene_names is not None:
            genebodies = genebodies[genebodies.gene_name.isin(gene_names)]
        return genebodies

    def get_exp_feather(self, peaks, extend_bp=300) -> pd.DataFrame:
        """Get the expression data for the given peaks. Only for backwards compatibility.

        Parameters
        ----------
        peaks: pd.DataFrame
            The peaks to query.
        extend_bp: int, optional
            The number of base pairs to extend the peaks, defaults to 300.

        Returns
        -------
        pd.DataFrame
            The expression data.
        """
        exp = (
            pr(peaks, int64=True)
            .join(pr(self.gtf, int64=True).extend(extend_bp), how="left")
            .as_df()
        )
        return exp.reset_index(drop=True)

    def query_region(self, chrom, start, end, strand=None) -> pd.DataFrame:
        """Query the GTF for regions matching the given parameters.

        Parameters
        ----------
        chrom: str
            The chromosome to query.
        start: int
            The start position to query.
        end: int
            The end position to query.
        strand: str, optional
            The strand to query, defaults to None.
        """
        result = self.gtf.query(
            f'Chromosome == "{chrom}" & Start > {start} & End < {end}'
        )
        if strand is not None:
            result = result.query(f'Strand == "{strand}"')
        return result
