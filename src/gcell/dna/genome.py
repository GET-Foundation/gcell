"""
Genome module for handling genomic sequences and regions.

This module provides classes and functions for working with genomic sequences,
regions, and annotations. It includes functionality for downloading genome files,
accessing sequences, and manipulating genomic regions.

Classes
-------
ChromSize: Handles chromosome size information
ChromGap: Handles chromosome gap information
Genome: Main class for genome sequence access and manipulation
GenomicRegion: Represents a single genomic region
GenomicRegionCollection: Collection of genomic regions with operations

Functions
---------
read_peaks: Read peak files in BED or narrowPeak format
read_blacklist: Read ENCODE blacklist regions
pandas_to_tabix_region: Convert pandas DataFrame to tabix region string
"""

import gzip
from io import StringIO
from pathlib import Path

import numpy as np
import pandas as pd
import pooch
import requests
from pandas.api.types import CategoricalDtype
from pyfaidx import Fasta
from pyranges import PyRanges
from scipy.sparse import csr_matrix

from .. import _settings
from .motif import print_results
from .sequence import DNASequence, DNASequenceCollection


class ChromSize:
    """
    Class for handling chromosome size information.

    Downloads and parses chromosome size files from UCSC genome browser.
    Provides methods to access and manipulate chromosome size data.

    Parameters
    ----------
    assembly : str
        Genome assembly name (e.g. 'hg38', 'mm10')
    annotation_dir : str or Path, optional
        Directory to store annotation files

    Methods
    -------
    get_dict(chr_included=None)
        Get chromosome sizes as dictionary
    save_chrom_sizes()
        Save chromosome sizes to file
    as_pyranges()
        Convert to PyRanges format
    tiling_region(tile_size, tile_overlap=0)
        Create tiled regions across chromosomes
    """

    def __init__(
        self, assembly=None, annotation_dir=_settings.get_setting("annotation_dir")
    ):
        self.assembly = assembly
        self.annotation_dir = Path(annotation_dir) if annotation_dir else None
        if self.assembly is None:
            raise ValueError("assembly is not specified")
        if self.annotation_dir is None:
            raise ValueError("annotation_dir is not specified")

        self.chrom_sizes = self.parse_or_download_chrom_sizes()

    def _download_chrom_sizes(self):
        url = f"http://hgdownload.soe.ucsc.edu/goldenPath/{self.assembly}/bigZips/{self.assembly}.chrom.sizes"
        response = requests.get(url)
        if response.status_code != 200:
            raise ConnectionError("Failed to download chromosome data")
        return self._parse_chrom_data(response.text)

    def _parse_chrom_data(self, data):
        chrom_sizes = {}
        lines = data.strip().split("\n")
        for line in lines:
            parts = line.split("\t")
            if len(parts) == 2:
                chrom, length = parts
                chrom_sizes[chrom] = int(length)
        return chrom_sizes

    def get_dict(self, chr_included=None):
        """
        Get chromosome sizes as dictionary.

        Parameters
        ----------
        chr_included : list, optional
            List of chromosome names to include, defaults to all chromosomes

        Returns
        -------
        dict
            Dictionary of chromosome sizes
        """
        if chr_included is None:
            return self.chrom_sizes
        else:
            return {chr: self.chrom_sizes.get(chr, None) for chr in chr_included}

    @property
    def dict(self):
        """
        Get chromosome sizes as dictionary
        """
        return self.chrom_sizes

    def save_chrom_sizes(self):
        """
        Save chromosome sizes to tab-delimited file
        """
        filepath = self.annotation_dir / f"{self.assembly}_chrom_sizes.txt"
        filepath.write_text(
            "\n".join(
                f"{chrom}\t{length}" for chrom, length in self.chrom_sizes.items()
            )
        )

    def parse_or_download_chrom_sizes(self):
        """
        Parse or download chromosome sizes
        """
        filepath = self.annotation_dir / f"{self.assembly}_chrom_sizes.txt"
        if filepath.exists():
            return self._parse_chrom_data(filepath.read_text())
        else:
            return self._download_chrom_sizes()

    def as_pyranges(self):
        """
        Convert chromosome sizes to PyRanges format
        """
        try:
            import pyranges as pr

            cs = pd.DataFrame(
                {
                    "Chromosome": list(self.chrom_sizes.keys()),
                    "Start": 0,
                    "End": list(self.chrom_sizes.values()),
                }
            ).sort_values(by=["Chromosome", "Start", "End"])
            return pr.PyRanges(cs, int64=True)
        except ImportError:
            raise ImportError("pyranges is not installed")

    def __repr__(self) -> str:
        return (
            f"ChromSize(assembly={self.assembly}, annotation_dir={self.annotation_dir})"
        )

    def tiling_region(self, tile_size: int, tile_overlap: int = 0):
        """
        Create tiled regions across chromosomes

        Parameters
        ----------
        tile_size : int
            The size of the tile
        tile_overlap : int, optional
            The overlap between tiles, defaults to 0

        Returns
        -------
        pandas.DataFrame
            DataFrame of tiled regions
        """
        pr_regions = self.as_pyranges()
        return pr_regions.tile(tile_size=tile_size, overlap=tile_overlap).as_df()


class ChromGap:
    """
    Class for handling chromosome gap information.

    Downloads and parses chromosome gap (AGP) files from UCSC genome browser.
    Provides methods to access gap annotations like telomeres and heterochromatin.

    Parameters
    ----------
    assembly : str
        Genome assembly name (e.g. 'hg38', 'mm10')
    annotation_dir : str or Path, optional
        Directory to store annotation files


    Methods
    -------
    get_telomeres(return_tabix=False)
        Get telomere regions
    get_heterochromatin(return_tabix=False)
        Get heterochromatin regions
    save_agp_data()
        Save AGP data to file
    """

    def __init__(
        self, assembly=None, annotation_dir=_settings.get_setting("annotation_dir")
    ):
        self.assembly = assembly
        self.annotation_dir = Path(annotation_dir) if annotation_dir else None

        if self.assembly is None:
            raise ValueError("assembly is not specified")
        if self.annotation_dir is None:
            raise ValueError("annotation_dir is not specified")

        self.agp_data = self.parse_or_download_agp()

    def _download_agp(self):
        url = f"https://hgdownload.soe.ucsc.edu/goldenPath/{self.assembly}/bigZips/{self.assembly}.agp.gz"
        response = requests.get(url)
        if response.status_code != 200:
            raise ConnectionError("Failed to download AGP data")
        return gzip.decompress(response.content).decode("utf-8")

    def _parse_agp_data(self, data):
        """
        Parse AGP data
        """
        columns = [
            "chrom",
            "start",
            "end",
            "part_number",
            "component_type",
            "component_id",
            "component_start",
            "component_end",
            "orientation",
        ]
        return pd.read_csv(StringIO(data), sep="\t", comment="#", names=columns).rename(
            columns={"chrom": "Chromosome", "start": "Start", "end": "End"}
        )

    def get_telomeres(self, return_tabix=False):
        """
        Get telomere regions

        Parameters
        ----------
        return_tabix : bool, optional
            Whether to return the regions in tabix format string (i.e. "chr1:100-200"), defaults to False
        """
        df = self.agp_data[self.agp_data["component_start"] == "telomere"]
        if return_tabix:
            return pandas_to_tabix_region(df)
        return df

    def get_heterochromatin(self, return_tabix=False):
        """
        Get heterochromatin regions

        Parameters
        ----------
        return_tabix : bool, optional
            Whether to return the regions in tabix format string (i.e. "chr1:100-200"), defaults to False
        """
        df = self.agp_data[
            self.agp_data["component_start"].isin(
                ["heterochromatin", "centromere", "telomere"]
            )
        ]
        if return_tabix:
            return pandas_to_tabix_region(df)
        return df

    def save_agp_data(self):
        """
        Save AGP data to file
        """
        filepath = self.annotation_dir / f"{self.assembly}_agp.txt"
        self.agp_data.to_csv(filepath, sep="\t", index=False)

    def parse_or_download_agp(self):
        """
        Parse or download AGP data
        """
        filepath = self.annotation_dir / f"{self.assembly}_agp.txt"
        if filepath.exists():
            return pd.read_csv(filepath, sep="\t")
        else:
            data = self._download_agp()
            agp_data = self._parse_agp_data(data)
            self.agp_data = agp_data
            self.save_agp_data()
            return agp_data

    def __repr__(self) -> str:
        return (
            f"ChromGap(assembly={self.assembly}, annotation_dir={self.annotation_dir})"
        )


class Genome:
    """
    Main class for accessing and manipulating genome sequences.

    Downloads genome sequence files and provides methods to access sequences
    and genomic regions. Handles chromosome naming conventions and coordinates.

    Parameters
    ----------
    assembly : str
        Genome assembly name (e.g. 'hg38', 'mm10')



    Methods
    -------
    get_sequence(chromosome, start, end, strand='+')
        Get DNA sequence for a genomic region, you can also use the sequence property with `self.genome_seq[chromosome][start:end]`
    random_draw(chromosome, length=4000000, strand='+')
        Get random genomic region
    tiling_region(chromosome, tile_size, step_size, strand='+')
        Create tiled regions across a chromosome
    normalize_chromosome(chromosome)
        Normalize chromosome name format
    """

    def __init__(self, assembly: str, load_genome_seq: bool = True) -> None:
        self.assembly = assembly
        self._download_files_if_not_exist(load_genome_seq)
        if load_genome_seq:
            fasta_file = (
                Path(_settings.get_setting("genome_dir")) / f"{self.assembly}.fa"
            )
            self.genome_seq = Fasta(str(fasta_file))

            if list(self.genome_seq.keys())[0].startswith("chr"):
                self.chr_suffix = "chr"
            else:
                self.chr_suffix = ""
        else:
            self.genome_seq = None
            self.chr_suffix = ""

    @property
    def annotation_dir(self) -> str:
        return str(Path(_settings.get_setting("annotation_dir")))

    @property
    def chrom_sizes(self) -> dict:
        """Get chromosome sizes, downloading if necessary"""
        return ChromSize(self.assembly, self.annotation_dir).dict

    @property
    def chrom_gap(self) -> dict:
        """Get chromosome gap, downloading if necessary"""
        return ChromGap(self.assembly, self.annotation_dir)

    @property
    def blacklist(self) -> dict:
        """Get blacklist, downloading if necessary"""
        if (
            Path(_settings.get_setting("annotation_dir"))
            / f"{self.assembly}-blacklist.v2.bed.gz"
        ).exists():
            return pd.read_csv(
                Path(_settings.get_setting("annotation_dir"))
                / f"{self.assembly}-blacklist.v2.bed.gz",
                sep="\t",
                header=None,
                names=["Chromosome", "Start", "End", "Category"],
            )
        else:
            return pd.DataFrame()

    @property
    def modeling_blacklist(self) -> pd.DataFrame:
        """Concat heterochromatin, blacklist if exists"""
        return pd.concat([self.chrom_gap.get_heterochromatin(), self.blacklist])[
            ["Chromosome", "Start", "End"]
        ]

    def _download_files_if_not_exist(self, download_genome_seq: bool = True) -> None:
        """Download genome files if they don't exist. This will download:
        - genome sequence
        - chromosome sizes
        - chromosome gap
        - blacklist
        """
        if download_genome_seq:
            fasta_file = (
                Path(_settings.get_setting("genome_dir")) / f"{self.assembly}.fa"
            )
            if not fasta_file.exists():
                # Define file name and URL
                fname = f"{self.assembly}.fa.gz"
                url = f"http://hgdownload.cse.ucsc.edu/goldenPath/{self.assembly}/bigZips/{self.assembly}.fa.gz"

                # Register the file with pooch
                _settings.POOCH.path = Path(_settings.get_setting("genome_dir"))
                _settings.POOCH.registry[fname] = (
                    None  # Will be updated with hash later
                )
                _settings.POOCH.urls[fname] = url

                # Download and decompress the file
                downloaded_file = _settings.POOCH.fetch(
                    fname, processor=pooch.Decompress()
                )

                # Move to final location
                Path(downloaded_file).rename(fasta_file)

        chrom_sizes_file = str(
            Path(_settings.get_setting("annotation_dir"))
            / f"{self.assembly}.chrom.sizes"
        )
        if not Path(chrom_sizes_file).exists():
            # Define file name and URL
            fname = f"{self.assembly}.chrom.sizes"
            url = f"http://hgdownload.cse.ucsc.edu/goldenPath/{self.assembly}/bigZips/{fname}"

            # Register the file with pooch
            _settings.POOCH.registry[fname] = None  # Will be updated with hash later
            _settings.POOCH.urls[fname] = url

            # Download chromosome sizes
            downloaded_file = _settings.POOCH.fetch(fname)

            # Move to final location
            Path(downloaded_file).rename(chrom_sizes_file)

        chrom_gap_file = str(
            Path(_settings.get_setting("annotation_dir")) / f"{self.assembly}.agp.gz"
        )
        if not Path(chrom_gap_file).exists():
            # Define file name and URL
            fname = f"{self.assembly}.agp.gz"
            url = f"http://hgdownload.cse.ucsc.edu/goldenPath/{self.assembly}/bigZips/{fname}"

            # Register the file with pooch
            _settings.POOCH.registry[fname] = None  # Will be updated with hash later
            _settings.POOCH.urls[fname] = url

            # Download chromosome gap
            downloaded_file = _settings.POOCH.fetch(fname)

            # Move to final location
            Path(downloaded_file).rename(chrom_gap_file)

        blacklist_file = str(
            Path(_settings.get_setting("annotation_dir"))
            / f"{self.assembly}-blacklist.v2.bed.gz"
        )
        if not Path(blacklist_file).exists():
            try:
                # Define file name and URL
                fname = f"{self.assembly}-blacklist.v2.bed.gz"
                url = f"https://raw.githubusercontent.com/Boyle-Lab/Blacklist/refs/heads/master/lists/{fname}"

                # Register the file with pooch
                _settings.POOCH.registry[fname] = (
                    None  # Will be updated with hash later
                )
                _settings.POOCH.urls[fname] = url

                # Download blacklist
                downloaded_file = _settings.POOCH.fetch(fname)

                # Move to final location
                Path(downloaded_file).rename(blacklist_file)
            except Exception as e:
                print(f"Failed to download blacklist: {e}, using empty blacklist")

    def __repr__(self) -> str:
        return f"Genome: {self.assembly}"

    def normalize_chromosome(self, chromosome):
        """
        Normalize chromosome name to use the specified chromosome suffix {self.chr_suffix}

        Parameters
        ----------
        chromosome : str
            The chromosome name
        """
        if str(chromosome).startswith("chr"):
            chromosome = str(chromosome)[3:]

        return self.chr_suffix + str(chromosome)

    def get_sequence(self, chromosome, start, end, strand="+"):
        """
        Get the sequence of the genomic region

        Parameters
        ----------
        chromosome : str
            The chromosome name
        start : int
            The start position (0-based)
        end : int
            The end position (exclusive)
        strand : str, optional
            The strand ('+' or '-'), defaults to '+'
        """
        if self.genome_seq is None:
            raise ValueError(
                "Genome sequence is not loaded, please set load_genome_seq=True when initializing the Genome object"
            )
        chromosome = self.normalize_chromosome(chromosome)
        if end > self.chrom_sizes[chromosome]:
            end = self.chrom_sizes[chromosome]
            print(
                f"""The end position {end} is larger than the chromosome size {self.chrom_sizes[chromosome]}
                The end position is set to the chromosome size"""
            )
        if start < 0:
            start = 0
            print(
                f"""The start position {start} is smaller than 0
                The start position is set to 0"""
            )
        if strand == "-":
            return DNASequence(
                self.genome_seq[chromosome][start:end].seq.complement(),
                header=f"{chromosome}_{start}_{end}",
            )
        else:
            return DNASequence(
                self.genome_seq[chromosome][start:end].seq,
                header=f"{chromosome}_{start}_{end}",
            )

    def random_draw(self, chromosome, length=4_000_000, strand="+"):
        """
        Randomly draw a genomic region with given length

        Parameters
        ----------
        chromosome : str
            The chromosome name
        length : int, optional
            The length of the genomic region, defaults to 4,000,000
        strand : str, optional
            The strand ('+' or '-'), defaults to '+'
        """
        if self.genome_seq is None:
            raise ValueError(
                "Genome sequence is not loaded, please set load_genome_seq=True when initializing the Genome object"
            )
        chromosome = self.normalize_chromosome(chromosome)
        start = np.random.randint(0, len(self.genome_seq[chromosome]) - length)
        return GenomicRegion(self, chromosome, start, start + length, strand)

    def get_chromosome_size(self, chromosome):
        """
        Get the size of the chromosome

        Parameters
        ----------
        chromosome : str
            The chromosome name
        """
        return self.chrom_sizes[chromosome]

    def tiling_region(self, chromosome, tile_size, step_size, strand="+"):
        """
        Tile the chromosome into smaller regions

        Parameters
        ----------
        chromosome : str
            The chromosome name
        tile_size : int
            The size of the tile
        step_size : int
            The step size between tiles
        strand : str, optional
            The strand ('+' or '-'), defaults to '+'
        """
        chromosome = self.normalize_chromosome(chromosome)
        return GenomicRegionCollection(
            self,
            chromosomes=[chromosome]
            * ((self.chrom_sizes[chromosome]) // step_size + 1),
            starts=range(0, self.chrom_sizes[chromosome], step_size),
            ends=[
                min(s + tile_size, self.chrom_sizes[chromosome])
                for s in range(0, self.chrom_sizes[chromosome], step_size)
            ],
            strands=[strand] * ((self.chrom_sizes[chromosome]) // step_size + 1),
        )


class GenomicRegion:
    """
    Class representing a single genomic region.

    Stores coordinates and provides methods to manipulate and analyze
    a genomic region.

    Parameters
    ----------
    genome : Genome
        Genome object containing sequence data
    chromosome : str
        Chromosome name
    start : int
        Start position (0-based)
    end : int
        End position (exclusive)
    strand : str, optional
        Strand ('+' or '-'), defaults to '+'

    Methods
    -------
    sequence
        Get DNA sequence for region
    get_motif_score(motif)
        Calculate motif scores
    lift_over(target_genome, lo)
        Convert coordinates to different assembly
    get_flanking_region(upstream, downstream)
        Get expanded region
    tiling_region(tile_size, step_size)
        Create tiled sub-regions
    """

    def __init__(
        self, genome: Genome, chromosome: str, start: int, end: int, strand: str = "+"
    ):
        self.genome = genome
        self.chromosome = chromosome
        self.start = start
        self.end = end
        self.strand = strand

    def __repr__(self) -> str:
        return f"[{self.genome.assembly}]{self.chromosome}:{self.start}-{self.end}"

    @property
    def sequence(self):
        """
        Get the sequence of the genomic region
        """
        return self.genome.get_sequence(
            self.chromosome, self.start, self.end, self.strand
        )

    def get_motif_score(self, motif):
        """
        Get the motif score of the genomic region
        """
        return "Not implemented yet"

    def lift_over(self, target_genome, lo):
        """
        Lift over the genomic region to another genome assembly using lo object
        User need to provide the correct lo object

        Parameters
        ----------
        target_genome : Genome
            The target genome assembly
        lo : Liftover
            The lo object
        """
        chromosome, start, end = lo.convert_coordinate(
            self.chromosome, self.start, self.end
        )
        if chromosome:
            return GenomicRegion(target_genome, chromosome, start, end, self.strand)
        else:
            return None

    def get_flanking_region(self, upstream, downstream):
        """
        Get the flanking region of the genomic region

        Parameters
        ----------
        upstream : int
            The number of bases upstream of the genomic region
        downstream : int
            The number of bases downstream of the genomic region
        """
        return GenomicRegion(
            self.genome,
            self.chromosome,
            self.start - upstream,
            self.end + downstream,
            self.strand,
        )

    def tiling_region(self, tile_size, step_size):
        """
        Tile the genomic region into smaller regions

        Parameters
        ----------
        tile_size : int
            The size of the tile
        step_size : int
            The step size between tiles
        """
        return GenomicRegionCollection(
            self.genome,
            chromosomes=[self.chromosome]
            * ((self.end - self.start - tile_size) // step_size + 1),
            starts=range(self.start, self.end, step_size),
            ends=[s + tile_size for s in range(self.start, self.end, step_size)],
            strands=[self.strand]
            * ((self.end - self.start - tile_size) // step_size + 1),
        )

    # def get_hic(self, hic, resolution=25000):
    #     """
    #     Get the Hi-C matrix of the genomic regions
    #     """
    #     # check if hicstraw is imported
    #     try:
    #         import hicstraw
    #     except ImportError:
    #         raise ImportError("hicstraw is not installed")

    #     start = self.start  # // resolution
    #     end = self.end  # // resolution + 1

    #     mzd = hic.getMatrixZoomData(
    #         self.chromosome.replace("chr", ""),
    #         self.chromosome.replace("chr", ""),
    #         "observed",
    #         "KR",
    #         "BP",
    #         resolution,
    #     )
    #     numpy_matrix = mzd.getRecordsAsMatrix(
    #         start,
    #         end,
    #         start,
    #         end,
    #         # start * resolution, end * resolution, start * resolution, end * resolution
    #     )
    #     dst = np.log10(numpy_matrix[1:, 1:] + 1)
    #     return dst


class GenomicRegionCollection(PyRanges):
    """List of GenomicRegion objects"""

    def __init__(
        self,
        genome,
        df=None,
        chromosomes=None,
        starts=None,
        ends=None,
        strands=None,
        int64=False,
        copy_df=True,
    ):
        super().__init__(df, chromosomes, starts, ends, strands, int64, copy_df)
        self.genome = genome

    def __repr__(self) -> str:
        return f"GenomicRegionCollection with {len(self.df)} regions"

    def to_bed(self, file_name):
        """
        Save the genomic region collection to a bed file
        """
        self.df[["Chromosome", "Start", "End", "Strand"]].to_csv(
            file_name, sep="\t", header=False, index=False
        )

    def center_expand(self, target_size):
        """
        Expand the genomic region collection from peak center
        """
        peak_center = (self.df["End"] + self.df["Start"]) // 2
        Start = peak_center - target_size // 2
        End = peak_center + target_size // 2
        index = self.df.index.values
        if "Strand" not in self.df.columns:
            return GenomicRegionCollection(
                genome=self.genome,
                df=pd.DataFrame(
                    {
                        "Chromosome": self.df["Chromosome"],
                        "Start": Start,
                        "End": End,
                        "Index": index,
                    }
                ),
            )
        else:
            return GenomicRegionCollection(
                genome=self.genome,
                df=pd.DataFrame(
                    {
                        "Chromosome": self.df["Chromosome"],
                        "Start": Start,
                        "End": End,
                        "Strand": self.df["Strand"],
                        "Index": index,
                    }
                ),
            )

    # generator of GenomicRegion objects
    def __iter__(self):
        for _, row in self.df.iterrows():
            if "Strand" not in row:
                yield GenomicRegion(
                    row["genome"], row["Chromosome"], row["Start"], row["End"], "+"
                )
            else:
                yield GenomicRegion(
                    row["genome"],
                    row["Chromosome"],
                    row["Start"],
                    row["End"],
                    row["Strand"],
                )

    def __getitem__(self, val):
        if isinstance(val, int):
            row = self.df.iloc[val]
            if "Strand" not in row:
                return GenomicRegion(
                    row["genome"], row["Chromosome"], row["Start"], row["End"], "+"
                )
            else:
                return GenomicRegion(
                    row["genome"],
                    row["Chromosome"],
                    row["Start"],
                    row["End"],
                    row["Strand"],
                )
        else:
            return GenomicRegionCollection(self.genome.iloc[val], self.df.iloc[val])

    # assign a new column to the pyranges
    def __setattr__(self, column_name, column):
        return super().__setattr__(column_name, column)

    def collect_sequence(
        self, mutations=None, upstream=0, downstream=0, target_length=None
    ):
        """
        Collect the sequence of the genomic regions
        """
        if mutations is None:
            return DNASequenceCollection(
                [
                    region.get_flanking_region(upstream, downstream).sequence
                    if target_length is None
                    else region.get_flanking_region(
                        upstream, downstream
                    ).sequence.padding(target_length=target_length)
                    for region in iter(self)
                ]
            )
        else:
            # first ensure the mutations are in the same genome
            if mutations.genome != self.genome:
                raise ValueError("The mutations are not in the same genome")
            # then collect the reference sequence
            ref_seq = self.collect_sequence(mutations=None, upstream=0, downstream=0)
            # overlap the mutations with the genomic regions to determine where the mutations are
            overlap = mutations.overlap(self)
            # calculate the relative position of the mutations in the genomic regions
            relative_pos = overlap.df["Start"].values - overlap.df["Start_b"].values
            # mutate the reference sequence for each mutation
            mutated_seq = [
                seq.mutate(
                    relative_pos[i],
                    mutations.df.iloc[i]["Reference"],
                    mutations.df.iloc[i]["Alternate"],
                )
                for i, seq in enumerate(ref_seq)
            ]

            return DNASequenceCollection(mutated_seq)

    def get_hic(self, hic, resolution=10000):
        """
        Get the Hi-C matrix of the genomic regions
        """
        start = self.df["Start"].min() // resolution
        end = self.df["End"].max() // resolution + 1
        hic_idx = np.array(
            [row.Start // resolution - start + 1 for _, row in self.df.iterrows()]
        )
        mzd = hic.getMatrixZoomData(
            self.df.iloc[0].Chromosome.replace("chr", ""),
            self.df.iloc[0].Chromosome.replace("chr", ""),
            "observed",
            "KR",
            "BP",
            resolution,
        )
        numpy_matrix = mzd.getRecordsAsMatrix(
            start * resolution, end * resolution, start * resolution, end * resolution
        )
        dst = np.log10(numpy_matrix[hic_idx, :][:, hic_idx] + 1)
        return dst

    def scan_motif(
        self,
        motifs,
        mutations=None,
        non_negative=True,
        upstream=0,
        downstream=0,
        raw=False,
    ):
        """
        Scan motif in sequence using MOODS.

        Parameters
        ----------
        seqs: List[Tuple[str, str]]
            A list of tuples containing the header and sequence for each input sequence.
        scanner: MOODS.Scanner
            The MOODS scanner to use for scanning the sequences.
        diff: bool, optional
            Whether to calculate the difference between the scores for the alternate and reference sequences. Defaults to False.

        Returns
        -------
        pd.DataFrame
            A dataframe containing the results of the scan. If `diff` is True, the dataframe will include columns for the difference in scores and the cluster of the motif.
        """
        seqs = self.collect_sequence(mutations, upstream, downstream)
        # initialize the output list
        output = []

        # scan each sequence and add the results to the output list
        headers = []
        lengths = []
        sequences = []
        for s in seqs:
            sequences.append(s.seq)
            headers.append(s.header)
            lengths.append(len(s.seq))

        # concatenate the sequences with 100 Ns between each sequence
        seq_cat = ("N" * 100).join(sequences)
        # get the list of sequence start and end positions in the concatenated sequence
        starts = np.cumsum([0] + lengths[:-1]) + 100 * np.arange(len(seqs))
        ends = starts + np.array(lengths)
        headers = np.array(headers)
        # scan the concatenated sequence
        results = motifs.scanner.scan(seq_cat)
        output = print_results(
            "", seq_cat, motifs.matrices, motifs.matrix_names, results
        )
        # convert the output list to a dataframe
        output = pd.DataFrame(
            output, columns=["header", "motif", "pos", "strand", "score", "seq"]
        )
        output["cluster"] = output.motif.map(motifs.motif_to_cluster)

        # assign header names in output dataframe based on 'pos' and starts/ends
        for i, h in enumerate(headers):
            output.loc[(output.pos >= starts[i]) & (output.pos < ends[i]), "header"] = h

        # remove the rows with multiple Ns
        output = output[~output.seq.str.contains("NN")]
        if raw:
            return output
        output = (
            output.groupby(["header", "pos", "cluster"])
            .score.max()
            .reset_index()
            .groupby(["header", "cluster"])
            .score.sum()
            .reset_index()
        )

        if non_negative:
            output.loc[output.score < 0, "score"] = 0

        motif_c = CategoricalDtype(categories=motifs.cluster_names, ordered=True)
        seq_c = CategoricalDtype(categories=headers, ordered=True)

        row = output.header.astype(seq_c).cat.codes
        col = output.cluster.astype(motif_c).cat.codes

        sparse_matrix = csr_matrix(
            (output["score"], (row, col)),
            shape=(seq_c.categories.size, motif_c.categories.size),
        )

        output = pd.DataFrame.sparse.from_spmatrix(
            sparse_matrix, index=seq_c.categories, columns=motif_c.categories
        )

        return output


def read_peaks(
    peak_file: str,
    return_collection: bool = False,
    genome: Genome | None = None,
    return_columns: list[str] | None = None,
) -> pd.DataFrame | GenomicRegionCollection:
    """
    Read genomic peak regions from BED or narrowPeak files.

    Automatically detects file format and number of columns.
    Can return either a pandas DataFrame or GenomicRegionCollection.

    Parameters
    ----------
    peak_file : str
        Path to peak file in BED or narrowPeak format
    return_collection : bool, optional
        If True, returns GenomicRegionCollection instead of DataFrame
    genome : Genome, optional
        Required if return_collection=True
    return_columns : list of str, optional
        Specific columns to return in DataFrame

    Returns
    -------
    Union[pd.DataFrame, GenomicRegionCollection]
        Peak regions as DataFrame or GenomicRegionCollection

    Raises
    ------
    ValueError
        If file format is invalid or genome not provided when needed
    """
    # Standard column names for different BED formats
    bed_columns = {
        3: ["chromosome", "start", "end"],
        4: ["chromosome", "start", "end", "name"],
        5: ["chromosome", "start", "end", "name", "score"],
        6: ["chromosome", "start", "end", "name", "score", "strand"],
        9: [
            "chromosome",
            "start",
            "end",
            "name",
            "score",
            "strand",
            "thickStart",
            "thickEnd",
            "rgb",
        ],
        12: [
            "chromosome",
            "start",
            "end",
            "name",
            "score",
            "strand",
            "thickStart",
            "thickEnd",
            "rgb",
            "blockCount",
            "blockSizes",
            "blockStarts",
        ],
    }

    narrowpeak_columns = [
        "chromosome",
        "start",
        "end",
        "name",
        "score",
        "strand",
        "signalValue",
        "pValue",
        "qValue",
        "summit",
    ]

    # Read first line to detect number of columns
    with Path(peak_file).open() as f:
        first_line = f.readline().strip()
        num_columns = len(first_line.split("\t"))

    if peak_file.endswith(".narrowPeak"):
        if num_columns != 10:
            raise ValueError(
                f"narrowPeak file should have 10 columns, found {num_columns}"
            )
        peaks = pd.read_csv(peak_file, sep="\t", header=None, names=narrowpeak_columns)
    elif peak_file.endswith(".bed"):
        if num_columns not in bed_columns:
            raise ValueError(f"Unsupported BED format with {num_columns} columns")
        peaks = pd.read_csv(
            peak_file, sep="\t", header=None, names=bed_columns[num_columns]
        )
        # check dtype of the fourth column, if quantitative, rename it to score
        if len(peaks.columns) == 4 and peaks.iloc[:, 3].dtype in [
            "float64",
            "int64",
            "float32",
            "int32",
            "float16",
            "int16",
        ]:
            peaks.rename(columns={3: "score"}, inplace=True)
    else:
        raise ValueError("Peak file must be .bed or .narrowPeak format")

    # Ensure required columns exist
    required_columns = return_columns or ["chromosome", "start", "end"]
    if not all(col in peaks.columns for col in required_columns):
        raise ValueError(f"Peak file must contain {required_columns} columns")

    if return_collection:
        if genome is None:
            raise ValueError("genome parameter is required when return_collection=True")

        # Convert to GenomicRegionCollection format
        df = peaks.rename(
            columns={"chromosome": "Chromosome", "start": "Start", "end": "End"}
        )

        # Add strand column if not present
        if "strand" in peaks.columns:
            df = df.rename(columns={"strand": "Strand"})
        elif "Strand" not in df.columns:
            df["Strand"] = "+"

        return GenomicRegionCollection(genome=genome, df=df)

    if return_columns:
        return peaks[return_columns]
    else:
        return peaks


def read_blacklist(genome: str) -> pd.DataFrame:
    """
    Read the blacklist regions from the ENCODE project.
    """
    supported_genomes = ["mm10", "ce10", "ce11", "dm3", "dm6", "hg19", "hg38"]
    if genome not in supported_genomes:
        raise ValueError(f"Only {supported_genomes} are supported")
    blacklist = pd.read_csv(
        f"https://raw.githubusercontent.com/Boyle-Lab/Blacklist/refs/heads/master/lists/{genome}-blacklist.v2.bed.gz",
        sep="\t",
        header=None,
    )
    blacklist.columns = ["Chromosome", "Start", "End", "name"]
    blacklist["Start"] += 1
    blacklist = PyRanges(blacklist).sort().df[["Chromosome", "Start", "End"]]
    return blacklist


def pandas_to_tabix_region(df, chrom_col="chrom", start_col="start", end_col="end"):
    return " ".join(
        df.apply(
            lambda x: f"{x[chrom_col]}:{x[start_col]}-{x[end_col]}", axis=1
        ).tolist()
    )
