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
        # Filter to transcript features only (keep original Start/End coordinates)
        gtf_df = gtf_df[gtf_df.Feature == "transcript"][
            [
                "Chromosome",
                "Start",
                "End",
                "Strand",
                "gene_name",
                "gene_id",
                "gene_type",
                "transcript_id",
            ]
        ].drop_duplicates().reset_index(drop=True)
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

        # Get access to original GTF data (works for both Gencode and GTF classes)
        # Use original_gtf property (works for Gencode) or _original_gtf attribute (works for GTF)
        try:
            # Try property first (works for Gencode class)
            original_gtf_df = self.original_gtf
        except AttributeError:
            # Fall back to attribute (works for GTF class)
            if hasattr(self, '_original_gtf'):
                original_gtf_df = self._original_gtf
            else:
                raise AttributeError(
                    "Cannot access original GTF data. Need either 'original_gtf' property (Gencode) or '_original_gtf' attribute (GTF)."
                )

        # Get all features for this gene from original GTF
        gene_original_df = original_gtf_df[original_gtf_df.gene_name == gene_name]

        # Calculate full gene body coordinates from original GTF (all features)
        # This represents the true genomic range of the gene (includes TSS and TES regions)
        # Works correctly for both strands:
        #   '+' strand: Start.min() = TSS region, End.max() = TES region
        #   '-' strand: Start.min() = TES region, End.max() = TSS region
        # Both give the full gene body span from first to last feature
        gene_start = gene_original_df.Start.min()
        gene_end = gene_original_df.End.max()

        # Construct tss_list from processed gtf (now has original Start/End)
        # TSS coordinates: '+' strand uses Start, '-' strand uses End
        strand = df.Strand.iloc[0]
        if strand == "+":
            tss_coords = df.Start.values
        elif strand == "-":
            tss_coords = df.End.values
        else:
            tss_coords = df.Start.values  # Fallback

        tss_list = pd.DataFrame({
            'Chromosome': df.Chromosome.values,
            'Start': tss_coords,
            'End': tss_coords,  # TSS is a single coordinate
            'Strand': df.Strand.values,
            'gene_name': df.gene_name.values,
            'gene_id': df.gene_id.values,
            'gene_type': df.gene_type.values if 'gene_type' in df.columns else [None] * len(tss_coords),
            'transcript_id': df.transcript_id.values if 'transcript_id' in df.columns else [None] * len(tss_coords),
        }).reset_index(drop=True)

        # Construct tes_list from processed gtf (now has original Start/End)
        # TES coordinates: '+' strand uses End, '-' strand uses Start
        if strand == "+":
            tes_coords = df.End.values
        elif strand == "-":
            tes_coords = df.Start.values
        else:
            tes_coords = df.End.values  # Fallback

        tes_list = pd.DataFrame({
            'Chromosome': df.Chromosome.values,
            'Start': tes_coords,
            'End': tes_coords,
            'Strand': df.Strand.values,
            'gene_name': df.gene_name.values,
            'gene_id': df.gene_id.values,
            'gene_type': df.gene_type.values if 'gene_type' in df.columns else [None] * len(tes_coords),
        }).reset_index(drop=True)

        return Gene(
            name=df.gene_name.iloc[0],
            id=df.gene_id.iloc[0],
            chrom=df.Chromosome.iloc[0],
            strand=strand,
            tss_list=tss_list,
            tes_list=tes_list,
            gene_start=gene_start,
            gene_end=gene_end,
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

        # Get access to original GTF data (works for both Gencode and GTF classes)
        # Use original_gtf property (works for Gencode) or _original_gtf attribute (works for GTF)
        try:
            # Try property first (works for Gencode class)
            original_gtf_df = self.original_gtf
        except AttributeError:
            # Fall back to attribute (works for GTF class)
            if hasattr(self, '_original_gtf'):
                original_gtf_df = self._original_gtf
            else:
                raise AttributeError(
                    "Cannot access original GTF data. Need either 'original_gtf' property (Gencode) or '_original_gtf' attribute (GTF)."
                )

        # Get all features for this gene from original GTF
        gene_original_df = original_gtf_df[original_gtf_df.gene_id.str.startswith(gene_id)]

        # Calculate full gene body coordinates from original GTF (all features)
        # This represents the true genomic range of the gene (includes TSS and TES regions)
        # Works correctly for both strands:
        #   '+' strand: Start.min() = TSS region, End.max() = TES region
        #   '-' strand: Start.min() = TES region, End.max() = TSS region
        # Both give the full gene body span from first to last feature
        gene_start = gene_original_df.Start.min()
        gene_end = gene_original_df.End.max()

        # Construct tss_list from processed gtf (now has original Start/End)
        # TSS coordinates: '+' strand uses Start, '-' strand uses End
        strand = df.Strand.iloc[0]
        if strand == "+":
            tss_coords = df.Start.values
        elif strand == "-":
            tss_coords = df.End.values
        else:
            tss_coords = df.Start.values  # Fallback

        tss_list = pd.DataFrame({
            'Chromosome': df.Chromosome.values,
            'Start': tss_coords,
            'End': tss_coords,  # TSS is a single coordinate
            'Strand': df.Strand.values,
            'gene_name': df.gene_name.values,
            'gene_id': df.gene_id.values,
            'gene_type': df.gene_type.values if 'gene_type' in df.columns else [None] * len(tss_coords),
            'transcript_id': df.transcript_id.values if 'transcript_id' in df.columns else [None] * len(tss_coords),
        }).reset_index(drop=True)

        # Construct tes_list from processed gtf (now has original Start/End)
        # TES coordinates: '+' strand uses End, '-' strand uses Start
        if strand == "+":
            tes_coords = df.End.values
        elif strand == "-":
            tes_coords = df.Start.values
        else:
            tes_coords = df.End.values  # Fallback

        tes_list = pd.DataFrame({
            'Chromosome': df.Chromosome.values,
            'Start': tes_coords,
            'End': tes_coords,
            'Strand': df.Strand.values,
            'gene_name': df.gene_name.values,
            'gene_id': df.gene_id.values,
            'gene_type': df.gene_type.values if 'gene_type' in df.columns else [None] * len(tes_coords),
        }).reset_index(drop=True)

        return Gene(
            name=df.gene_name.iloc[0],
            id=df.gene_id.iloc[0],
            chrom=df.Chromosome.iloc[0],
            strand=strand,
            tss_list=tss_list,
            tes_list=tes_list,
            gene_start=gene_start,
            gene_end=gene_end,
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
