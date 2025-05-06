"""
Module for handling gene data.

Classes
-------
Gene: A class to represent a gene with TSS information.
TSS: A class to represent a transcription start site (TSS).
GeneExp: A class to represent a gene with expression data.
GeneSets: A class to represent a collection of genes.
"""

from collections.abc import Collection
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd
from tqdm import tqdm

from ..dna.track import Track


class TSS:
    """A class to represent a transcription start site (TSS)."""

    def __init__(self, name, peak_id, chrom, start, strand) -> None:
        """Initialize the TSS class.

        Parameters
        ----------
        name : str
            The name of the TSS.
        peak_id : int
            The ID of the TSS.
        chrom : str
            The chromosome of the TSS.
        start : int
            The start position of the TSS.
        strand : str
            The strand of the TSS.
        """
        self.name = name
        self.peak_id = peak_id
        self.chrom = chrom
        self.start = start
        self.strand = strand

    def __repr__(self) -> str:
        return f"TSS(name={self.name}, peak_id={self.peak_id}, chrom={self.chrom}, strand={self.strand}, start={str(self.start)})"

    def get_sample_from_peak(self, peak_df, focus=100) -> pd.DataFrame:
        """Get the sample from the peak_df, with the peak_id as the center, and focus as the window size

        Parameters
        ----------
        peak_df: pd.DataFrame
            The peak dataframe.
        focus: int
            The window size.

        Returns
        -------
        pd.DataFrame
            The sample from the peak_df.
        """
        return peak_df.iloc[self.peak_id - focus : self.peak_id + focus]


class Gene:
    """A class to represent a gene with TSS information."""

    def __init__(self, name, id, chrom, strand, tss_list) -> None:
        """Initialize the Gene class.

        Parameters
        ----------
        name
            The name of the gene.
        id
            The ID of the gene.
        chrom
            The chromosome of the gene.
        strand
            The strand of the gene.
        tss_list
            The TSS list of the gene.
        """
        self.name = name
        self.id = id
        self.chrom = chrom
        self.strand = strand

        self.tss_list = tss_list

    def __repr__(self) -> str:
        return "Gene(name={}, id={}, chrom={}, strand={}, tss_list={})".format(
            self.name,
            self.id,
            self.chrom,
            self.strand,
            ",".join(self.tss_list.Start.values.astype(str)),
        )

    @property
    def tss(self) -> list[TSS]:
        """Get the TSS list for the gene.

        Returns
        -------
        list[TSS]
            The list of TSS objects.
        """
        return [
            TSS(self.name, self.id, self.chrom, start, self.strand)
            for start in self.tss_list.Start.values
        ]

    @property
    def genomic_range(
        self, upstream=128 * 8192, downstream=128 * 8192
    ) -> tuple[str, int, int, str]:
        return (
            self.chrom,
            self.tss_list.Start.min() - upstream,
            self.tss_list.Start.min() + downstream,
            self.strand,
        )

    def get_track(self, track, upstream=128 * 8192, downstream=128 * 8192, **kwargs):
        return track.get_track(
            chr_name=self.chrom,
            start=self.tss_list.Start.min() - upstream,
            end=self.tss_list.Start.min() + downstream,
            **kwargs,
        )

    def get_tss_track(self, track, upstream=1000, downstream=1000, **kwargs):
        if self.strand == "+":
            return track.get_track(
                chr_name=self.chrom,
                start=self.tss_list.Start.min() - upstream,
                end=self.tss_list.Start.min() + downstream,
                **kwargs,
            )
        else:
            return track.get_track(
                chr_name=self.chrom,
                start=self.tss_list.Start.max() - downstream,
                end=self.tss_list.Start.max() + upstream,
                **kwargs,
            )

    def get_track_obj(
        self, track, upstream=128 * 8192, downstream=128 * 8192, **kwargs
    ) -> Track:
        return track.get_track_obj(
            chr_name=self.chrom,
            start=self.tss_list.Start.min() - upstream,
            end=self.tss_list.Start.min() + downstream,
            **kwargs,
        )

    def get_tss_track_obj(
        self, track, upstream=1000, downstream=1000, **kwargs
    ) -> Track:
        if self.strand == "+":
            return track.get_track_obj(
                chr_name=self.chrom,
                start=self.tss_list.Start.min() - upstream,
                end=self.tss_list.Start.min() + downstream,
                **kwargs,
            )
        else:
            return track.get_track_obj(
                chr_name=self.chrom,
                start=self.tss_list.Start.max() - downstream,
                end=self.tss_list.Start.max() + upstream,
                **kwargs,
            )


class GeneExp(Gene):
    """A class to represent a gene with expression data. Not very useful."""

    def __init__(self, name, id, chrom, strand, tss_list, exp_list) -> None:
        """Initialize the GeneExp class.

        Parameters
        ----------
        name : str
            The name of the gene.
        id : str
            The ID of the gene.
        chrom : str
            The chromosome of the gene.
        strand : str
            The strand of the gene.
        tss_list : pd.DataFrame
            The TSS list of the gene.
        exp_list : pd.DataFrame
            The expression list of the gene.
        """
        super().__init__(name, id, chrom, strand, tss_list)
        self.exp_list = exp_list

    def __repr__(self) -> str:
        return (
            "GeneExp(name={}, id={}, chrom={}, strand={}, tss_list={}, exp={})".format(
                self.name,
                self.id,
                self.chrom,
                self.strand,
                ",".join(self.tss_list.Start.values.astype(str)),
                self.exp_list.mean(),
            )
        )


class GeneSets(Collection):
    """A collection of Genes, initialized from a list of Gene objects or Gene names"""

    def __init__(self, genes: Collection[Gene]) -> None:
        """Initialize the GeneSets class.

        Parameters
        ----------
        genes : Collection[Gene]
            The list of Gene objects.
        """
        self.gene_names = [gene.name for gene in genes]
        self.gene_ids = [gene.id for gene in genes]
        self.data = dict(zip(self.gene_names, genes))
        self.tss_list = pd.concat(
            [gene.tss_list for gene in genes], axis=0
        ).reset_index(drop=True)

    def __contains__(self, x: object) -> bool:
        return x in self.genes

    def __iter__(self):
        return iter(self.genes)

    def __len__(self) -> int:
        return len(self.genes)

    def __repr__(self) -> str:
        return "GeneSets(gene_names={})".format(",".join(self.gene_names))

    def get_tss_track(
        self, track, upstream=1000, downstream=1000, n_jobs=96, **kwargs
    ) -> np.ndarray:
        """Get the TSS track for the gene sets.

        Parameters
        ----------
        track: Track
            The track object.
        upstream: int
            The upstream window size.
        downstream: int
            The downstream window size.
        n_jobs: int
            The number of jobs to run in parallel.
        """
        results = []
        tss_df = []
        with ProcessPoolExecutor(n_jobs) as executor:
            futures = []
            for chr_name in self.tss_list.Chromosome.unique():
                region_list = self.tss_list[self.tss_list.Chromosome == chr_name]
                region_list["Start"] = region_list["Start"] - upstream
                region_list["End"] = region_list["End"] + downstream
                # if strand is +, then get first TSS, else get last TSS
                region_list = region_list.sort_values("Start", ascending=True)
                pos = (
                    region_list.query('Strand == "+"')
                    .groupby("gene_name")
                    .first()
                    .reset_index()
                )
                neg = (
                    region_list.query('Strand == "-"')
                    .groupby("gene_name")
                    .last()
                    .reset_index()
                )
                region_list = pd.concat([pos, neg], axis=0)

                # region_list = pr(region_list).merge().df
                # tss_df.append(region_list)
                # # get center of each region and expand to 2000bp
                # region_list['Start'] = region_list.Start + \
                #     (region_list.End-region_list.Start)//2 - 1000
                # region_list['End'] = region_list.Start + 2000
                tss_df.append(region_list)
                region_list = region_list[["Start", "End"]].values
                # print(region_list.shape)
                # split into 100 region chunks
                if len(region_list) > 50:
                    region_list = np.array_split(region_list, len(region_list) // 4)
                    for chunk in region_list:
                        future = executor.submit(
                            track.get_track_for_regions,
                            chr_name=chr_name,
                            region_list=chunk,
                            **kwargs,
                        )
                        futures.append(future)
                elif len(region_list) == 0:
                    continue
                else:
                    future = executor.submit(
                        track.get_track_for_regions,
                        chr_name=chr_name,
                        region_list=region_list,
                        **kwargs,
                    )
                    futures.append(future)

            for future in tqdm(futures):
                result = future.result()
                if isinstance(result, list):
                    result = np.array(result)
                elif isinstance(result, np.ndarray) and result.ndim == 2:
                    result = result.sum(0)
                results.append(result)

        results = np.vstack(results)
        tss_df = pd.concat(tss_df, axis=0).reset_index(drop=True)
        return results, tss_df
