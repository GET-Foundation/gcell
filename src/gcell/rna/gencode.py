"""
Module for handling GENCODE gene annotations.

Classes
-------
Gencode: A class to handle GENCODE gene annotations for different genome assemblies.

"""

from pathlib import Path

import pandas as pd

from .. import _settings
from .gtf import GTF


class Gencode(GTF):
    """A class to handle GENCODE gene annotations for different genome assemblies.

    GENCODE (https://www.gencodegenes.org/) provides high-quality gene annotations
    for human and mouse genomes. This class handles downloading, parsing, and
    accessing GENCODE annotation data.

    The class inherits from GTF and provides convenient methods to access gene
    information like strand, chromosome location, transcription start sites, etc.

    """

    @classmethod
    def from_config(cls, config):
        """Create a Gencode instance from a configuration dictionary.

        Args:
            config (dict): Configuration dictionary containing 'assembly' and
                'gencode_version' keys.

        Returns:
            Gencode: A new Gencode instance initialized with the config values.
        """
        return cls(
            assembly=config.get("assembly"),
            version=config.get("gencode_version"),
        )

    def __init__(
        self, assembly="hg38", version=40, is_basic=True, exclude_chrs=["chrM", "chrY"]
    ):
        """Initialize the Gencode class.

        Parameters
        ----------
        assembly : str
            Genome assembly version. Options: "hg38" (human), "mm10" (mouse), or "hg19" (human). Defaults to "hg38".
        version : int
            GENCODE release version. Defaults to 40.
        exclude_chrs : list
            Chromosomes to exclude from the annotation. Defaults to ["chrM", "chrY"].
        """
        self.assembly = assembly
        self.version = version
        self.is_basic = is_basic
        self.gtf_dir = Path(_settings.get_setting("annotation_dir"))

        # Download annotation if needed
        self._download_files_if_not_exist()

        feather_path = self.gtf_dir / f"gencode.v{str(version)}.{self.assembly}.feather"

        if feather_path.exists():
            self.gtf = pd.read_feather(feather_path)
            self.feather_file = feather_path
        else:
            super().__init__(self.gtf_path, exclude_chrs)
            self.gtf.to_feather(feather_path)
            self.feather_file = feather_path

    def _download_files_if_not_exist(self) -> None:
        """Download GTF file from GENCODE if it doesn't exist locally.

        Downloads the basic annotation GTF file for the specified assembly and version
        from the GENCODE FTP server. The file is cached locally for future use.

        The URLs are constructed based on the assembly type:
        - hg38: Current human genome assembly
        - mm10: Mouse genome assembly
        - hg19: Previous human genome assembly (GRCh37)
        """
        # Set up URLs based on assembly
        if self.assembly == "hg38":
            fname = f"gencode.v{self.version}.basic.annotation.gtf.gz"
            url = f"http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{self.version}/{fname}"
        elif self.assembly == "mm10":
            fname = (
                f"gencode.{self.assembly}.v{str(self.version)}.basic.annotation.gtf.gz"
            )
            url = f"http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_{self.version}/gencode.v{self.version}.basic.annotation.gtf.gz"
        elif self.assembly == "hg19":
            fname = f"gencode.v{str(self.version)}lift37.basic.annotation.gtf.gz"
            url = f"https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{self.version}/GRCh37_mapping/{fname}"

        if not self.is_basic:
            # replace .basic with ''
            fname = fname.replace(".basic", "")
            url = url.replace(".basic", "")

        self.gtf_path = Path(_settings.get_setting("annotation_dir")) / fname
        if not self.gtf_path.exists():
            _settings.download_with_pooch(fname, url, target_dir="annotation_dir")

    @property
    def original_gtf(self):
        if hasattr(self, "_original_gtf"):
            return self._original_gtf
        else:
            from .gtf import read_gtf

            self._original_gtf = read_gtf(str(self.gtf_path)).as_df()
            return self._original_gtf

    @property
    def gene_to_strand(self):
        """Dict mapping gene names to their strand orientation ('+' or '-').

        Returns:
            dict: Dictionary with gene names as keys and strand symbols as values
        """
        return self.gtf.set_index("gene_name")["Strand"].to_dict()

    @property
    def gene_to_strand_numeric(self):
        """Dict mapping gene names to numeric strand values (0 for '+', 1 for '-').

        Returns:
            dict: Dictionary with gene names as keys and numeric strand values
        """
        gtf_copy = self.gtf.copy()
        gtf_copy["Strand"] = gtf_copy["Strand"].map({"+": 0, "-": 1})
        return gtf_copy.set_index("gene_name")["Strand"].to_dict()

    @property
    def gene_to_id(self):
        """Dict mapping gene names to their GENCODE IDs.

        Returns:
            dict: Dictionary with gene names as keys and GENCODE IDs as values
        """
        return self.gtf.set_index("gene_name")["gene_id"].to_dict()

    @property
    def id_to_gene(self):
        """Dict mapping GENCODE IDs to gene names.

        Returns:
            dict: Dictionary with GENCODE IDs as keys and gene names as values
        """
        return self.gtf.set_index("gene_id")["gene_name"].to_dict()

    @property
    def gene_to_chrom(self):
        """Dict mapping gene names to their chromosome locations.

        Returns:
            dict: Dictionary with gene names as keys and chromosome names as values
        """
        return self.gtf.set_index("gene_name")["Chromosome"].to_dict()

    @property
    def gene_to_tss(self):
        """Dict mapping gene names to their transcription start sites (TSS).

        For '+' strand genes, TSS is at Start coordinate.
        For '-' strand genes, TSS is at End coordinate.

        Returns:
            dict: Dictionary with gene names as keys and TSS positions as values
        """
        # Group by gene_name and calculate TSS based on strand
        gene_groups = self.gtf.groupby("gene_name")
        tss_dict = {}

        for gene_name, group in gene_groups:
            strand = group["Strand"].iloc[0]
            if strand == "+":
                # For '+' strand, TSS is at Start coordinate
                tss_dict[gene_name] = group["Start"].min()
            elif strand == "-":
                # For '-' strand, TSS is at End coordinate
                tss_dict[gene_name] = group["End"].max()
            else:
                # Unstranded case - use Start.min() as fallback
                tss_dict[gene_name] = group["Start"].min()

        return tss_dict

    @property
    def gene_to_tes(self):
        """Dict mapping gene names to their transcription end sites (TES).

        For '+' strand genes, TES is at the maximum End coordinate.
        For '-' strand genes, TES is at the minimum Start coordinate.

        Returns:
            dict: Dictionary with gene names as keys and TES positions as values
        """
        # Group by gene_name and calculate TES based on strand
        gene_groups = self.gtf.groupby("gene_name")
        tes_dict = {}

        for gene_name, group in gene_groups:
            strand = group["Strand"].iloc[0]
            if strand == "+":
                # For '+' strand, TES is at the maximum End coordinate
                tes_dict[gene_name] = group["End"].max()
            elif strand == "-":
                # For '-' strand, TES is at the minimum Start coordinate
                tes_dict[gene_name] = group["Start"].min()
            else:
                # Unstranded case - use End.max() as fallback
                tes_dict[gene_name] = group["End"].max()

        return tes_dict

    @property
    def gene_to_type(self):
        """Dict mapping gene names to their gene types (e.g., protein_coding, lincRNA).

        Returns:
            dict: Dictionary with gene names as keys and gene types as values
        """
        return self.gtf.set_index("gene_name")["gene_type"].to_dict()
