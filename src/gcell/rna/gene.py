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

    def __init__(self, name, id, chrom, strand, tss_list, tes_list=None, gene_start=None, gene_end=None) -> None:
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
            The TSS list of the gene (DataFrame with Start and End columns).
        tes_list
            Optional TES list of the gene (DataFrame with Start and End columns).
            If not provided, will be constructed from tss_list on first access.
            For backward compatibility, defaults to None.
        gene_start
            Optional gene start coordinate from original GTF (full gene body).
            If not provided, will be calculated from tss_list as fallback.
            Defaults to None.
        gene_end
            Optional gene end coordinate from original GTF (full gene body).
            If not provided, will be calculated from tss_list as fallback.
            Defaults to None.
        """
        self.name = name
        self.id = id
        self.chrom = chrom
        self.strand = strand

        self.tss_list = tss_list
        self._tes_list = tes_list  # Can be provided or constructed on first access
        self._gene_start = gene_start
        self._gene_end = gene_end

    def __repr__(self) -> str:
        tss_coords = ",".join(self.tss_list.Start.values.astype(str))
        try:
            tes_coords = ",".join(self.tes_list.Start.values.astype(str))
        except (ValueError, AttributeError):
            tes_coords = "N/A"
        return "Gene(name={}, id={}, chrom={}, strand={}, tss={}, tes={})".format(
            self.name,
            self.id,
            self.chrom,
            self.strand,
            tss_coords,
            tes_coords,
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
    def tss_coordinate(self) -> int:
        """Get the primary Transcription Start Site (TSS) coordinate for the gene.

        Returns the primary TSS coordinate, symmetric to tes property.
        For '+' strand genes, returns the minimum Start coordinate.
        For '-' strand genes, returns the maximum Start coordinate.

        Returns
        -------
        int
            The TSS coordinate.
        """
        if self.strand == "+":
            return self.tss_list.Start.min()  # Min Start coordinate for '+' strand
        elif self.strand == "-":
            return self.tss_list.Start.max()  # Max Start coordinate for '-' strand
        else:
            return self.tss_list.Start.min()  # Fallback

    @property
    def tes_list(self) -> pd.DataFrame:
        """Get the TES list for the gene, constructed from tss_list if not provided.

        For '+' strand genes, TES coordinates are extracted from End column.
        For '-' strand genes, TES coordinates are extracted from Start column.

        Returns
        -------
        pd.DataFrame
            DataFrame containing TES coordinates, similar structure to tss_list.

        Raises
        ------
        ValueError
            If tss_list doesn't have End column and tes_list was not provided.
        """
        if self._tes_list is not None:
            return self._tes_list

        # Construct from tss_list if it has End column
        if 'End' not in self.tss_list.columns:
            raise ValueError(
                f"Cannot construct tes_list for gene {self.name}: "
                f"tss_list doesn't have End column. Provide tes_list parameter "
                f"when creating Gene object, or ensure tss_list includes End column."
            )

        if self.strand == "+":
            # For '+' strand, TES is at End coordinates
            tes_coords = self.tss_list.End.values
        elif self.strand == "-":
            # For '-' strand, TES is at Start coordinates (5' end)
            tes_coords = self.tss_list.Start.values
        else:
            # Unstranded - use End as fallback
            tes_coords = self.tss_list.End.values

        # Construct tes_list DataFrame similar to tss_list structure
        self._tes_list = pd.DataFrame({
            'Chromosome': self.tss_list.Chromosome.values,
            'Start': tes_coords,
            'End': tes_coords,  # TES is a single coordinate
            'Strand': self.tss_list.Strand.values,
            'gene_name': self.tss_list.gene_name.values,
            'gene_id': self.tss_list.gene_id.values,
            'gene_type': self.tss_list.gene_type.values if 'gene_type' in self.tss_list.columns else [None] * len(tes_coords),
        }).reset_index(drop=True)

        return self._tes_list

    @property
    def tes(self) -> int:
        """Get the Transcription End Site (TES) coordinate for the gene.

        Returns the primary TES coordinate, similar to how TSS works.
        For '+' strand genes, returns the maximum End coordinate.
        For '-' strand genes, returns the minimum Start coordinate.

        Returns
        -------
        int
            The TES coordinate.
        """
        tes_list = self.tes_list
        if self.strand == "+":
            return tes_list.Start.max()  # Max End coordinate for '+' strand
        elif self.strand == "-":
            return tes_list.Start.min()  # Min Start coordinate for '-' strand
        else:
            return tes_list.Start.max()  # Fallback

    @property
    def genomic_range(
        self,
    ) -> tuple[str, int, int, str]:
        """Get the genomic range of the gene (start, end coordinates).

        Uses the full gene body coordinates from the original GTF.
        Requires gene_start and gene_end to be provided during Gene initialization.

        Returns
        -------
        tuple[str, int, int, str]
            (chromosome, start, end, strand)

        Raises
        ------
        ValueError
            If gene_start or gene_end are not available.
        """
        if self._gene_start is not None and self._gene_end is not None:
            return (
                self.chrom,
                self._gene_start,
                self._gene_end,
                self.strand,
            )
        raise ValueError(
            f"genomic_range not available for gene {self.name}: "
            f"gene_start and gene_end must be provided during Gene initialization. "
            f"Use gtf.get_gene() or gencode.get_gene() to create Gene objects with full genomic range."
        )


    def get_track(self, track, upstream=1000, downstream=1000, **kwargs):
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

    def get_track_obj(self, track, upstream=1000, downstream=1000, **kwargs) -> Track:
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

    def get_tes_track(self, track, upstream=1000, downstream=1000, **kwargs):
        """Get track data centered around the gene's Transcription End Site (TES).

        Parameters
        ----------
        track
            The track object to retrieve data from.
        upstream : int, optional
            Number of base pairs upstream of the TES to include.
            For '+' strand genes, this extends before the TES (toward the gene body).
            For '-' strand genes, this extends after the TES (away from the gene body).
            Defaults to 1000.
        downstream : int, optional
            Number of base pairs downstream of the TES to include.
            For '+' strand genes, this extends after the TES (away from the gene body).
            For '-' strand genes, this extends before the TES (toward the gene body).
            Defaults to 1000.
        **kwargs
            Additional arguments passed to track.get_track.

        Returns
        -------
        Track or other return type from track.get_track
        """
        tes_list = self.tes_list
        if self.strand == "+":
            return track.get_track(
                chr_name=self.chrom,
                start=tes_list.Start.max() - upstream,
                end=tes_list.Start.max() + downstream,
                **kwargs,
            )
        else:
            # For '-' strand, upstream/downstream are reversed relative to transcription direction
            return track.get_track(
                chr_name=self.chrom,
                start=tes_list.Start.min() - downstream,
                end=tes_list.Start.min() + upstream,
                **kwargs,
            )

    def get_tes_track_obj(
        self, track, upstream=1000, downstream=1000, **kwargs
    ) -> Track:
        """Get Track object centered around the gene's Transcription End Site (TES).

        Parameters
        ----------
        track
            The track object to retrieve data from.
        upstream : int, optional
            Number of base pairs upstream of the TES to include.
            For '+' strand genes, this extends before the TES (toward the gene body).
            For '-' strand genes, this extends after the TES (away from the gene body).
            Defaults to 1000.
        downstream : int, optional
            Number of base pairs downstream of the TES to include.
            For '+' strand genes, this extends after the TES (away from the gene body).
            For '-' strand genes, this extends before the TES (toward the gene body).
            Defaults to 1000.
        **kwargs
            Additional arguments passed to track.get_track_obj.

        Returns
        -------
        Track
            Track object for the TES region.
        """
        tes_list = self.tes_list
        if self.strand == "+":
            return track.get_track_obj(
                chr_name=self.chrom,
                start=tes_list.Start.max() - upstream,
                end=tes_list.Start.max() + downstream,
                **kwargs,
            )
        else:
            # For '-' strand, upstream/downstream are reversed relative to transcription direction
            return track.get_track_obj(
                chr_name=self.chrom,
                start=tes_list.Start.min() - downstream,
                end=tes_list.Start.min() + upstream,
                **kwargs,
            )


class GeneExp(Gene):
    """A class to represent a gene with expression data. Not very useful."""

    def __init__(self, name, id, chrom, strand, tss_list, exp_list, tes_list=None, gene_start=None, gene_end=None) -> None:
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
        tes_list : pd.DataFrame, optional
            Optional TES list of the gene. If not provided, will be constructed
            from tss_list on first access. Defaults to None.
        gene_start : int, optional
            Optional gene start coordinate from original GTF (full gene body).
            Defaults to None.
        gene_end : int, optional
            Optional gene end coordinate from original GTF (full gene body).
            Defaults to None.
        """
        super().__init__(name, id, chrom, strand, tss_list, tes_list=tes_list, gene_start=gene_start, gene_end=gene_end)
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
