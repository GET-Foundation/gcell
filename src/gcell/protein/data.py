"""
Protein Data Management Module
:module: gcell.protein.data

This module provides functionality for managing and accessing protein-related data, including:

* UniProt sequences
* Gene name to UniProt ID mappings
* AlphaFold pLDDT scores
* UniProt XML schema

The module implements a singleton pattern through the ProteinData class to ensure
efficient memory usage when handling large protein datasets.

Key Features:
    * Automatic downloading of required protein data files
    * Lazy loading of data (only loaded when needed)
    * Mapping between different protein identifiers
    * Access to protein sequences and structure confidence scores

Supported Organisms:
    * Human (HUMAN_9606)
    * Mouse (MOUSE_10090)
    * Rat (RAT_10116)
    * Yeast (YEAST_559292)
    * And others (see organism_to_uniprot dictionary)

Example:
    Basic usage for getting protein sequences::

        >>> from gcell.protein.data import get_seq_from_gene_name
        >>> sequence = get_seq_from_gene_name("TP53")
        >>> from gcell.protein.data import get_lddt_from_uniprot_id
        >>> lddt_scores = get_lddt_from_uniprot_id("P04637")

Dependencies:
    * BioPython
    * numpy
    * pandas
    * xmlschema
    * tqdm

.. note::
    All data is loaded lazily to minimize memory usage until needed.
"""

import gzip
import os
import pickle
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
import xmlschema
from Bio import SeqIO
from tqdm import tqdm

from .. import _settings

# Organism mapping
organism_to_uniprot = {
    "human": "HUMAN_9606",
    "mouse": "MOUSE_10090",
    "rat": "RAT_10116",
    "yeast": "YEAST_559292",
    "chicken": "CHICK_9031",
    "drosophila": "DROME_7227",
    "ecoli": "ECOLI_83333",
    "danio": "DANRE_7955",
    "arabidopsis": "ARATH_3702",
    "caenorhabditis": "CAEEL_6239",
    "schizosaccharomyces": "SCHPO_284812",
}

# UCSC reference genome mapping to organism
ucsc_to_organism = {
    "hg38": "human",
    "mm10": "mouse",
    "hg19": "human",
    "mm9": "mouse",
    "rn6": "rat",
    "danRer11": "danio",
    "dm6": "drosophila",
    "ce11": "caenorhabditis",
    "sacCer3": "yeast",
}


def get_protein_data_dir():
    """
    Get the directory path for protein data files.

    Returns:
        Path: Path object pointing to the protein data directory
    """
    return Path(_settings.get_setting("annotation_dir"))


def download_protein_files():
    """
    Download required protein data files if they don't exist locally.

    Downloads:
        - UniProt Swiss-Prot FASTA database
        - Human UniProt ID mapping data
        - AlphaFold pLDDT scores
    """
    data_dir = get_protein_data_dir()
    data_dir.mkdir(parents=True, exist_ok=True)

    # Define files to download
    files = {
        "uniprot_sprot.fasta.gz": "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz",
        "HUMAN_9606_idmapping.dat.gz": "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping.dat.gz",
        "9606.pLDDT.tdt.zip": "https://github.com/normandavey/ProcessedAlphafold/raw/main/9606.pLDDT.tdt.zip",
    }

    for fname, url in files.items():
        _settings.download_with_pooch(fname, url, target_dir="annotation_dir")


class ProteinData:
    """
    Singleton class to manage protein data with lazy loading capabilities.

    This class implements a singleton pattern to ensure only one instance of the
    protein data is loaded in memory. Data is loaded only when first accessed.

    """

    _instance = None

    def __new__(cls):
        if cls._instance is None:
            cls._instance = super().__new__(cls)
            cls._instance._initialized = False
        return cls._instance

    def __init__(self):
        if not self._initialized:
            self._seq = None
            self._genename_to_uniprot = None
            self._lddt = None
            self._schema = None
            self._initialized = True

    def initialize_all(self):
        """Initialize all data if not already loaded"""
        download_protein_files()
        _ = self.schema  # This will trigger _initialize_schema
        _ = self.seq  # This will trigger _initialize_seq
        _ = (
            self.genename_to_uniprot
        )  # This will trigger _initialize_genename_to_uniprot
        _ = self.lddt  # This will trigger _initialize_lddt

    @property
    def seq(self):
        if self._seq is None:
            self._initialize_seq()
        return self._seq

    @property
    def genename_to_uniprot(self):
        if self._genename_to_uniprot is None:
            self._initialize_genename_to_uniprot()
        return self._genename_to_uniprot

    @property
    def lddt(self):
        if self._lddt is None:
            self._initialize_lddt()
        return self._lddt

    @property
    def schema(self):
        if self._schema is None:
            self._initialize_schema()
        return self._schema

    def _initialize_schema(self):
        """Initialize UniProt XML schema"""
        schema_path = os.environ.get("GCELL_UNIPROT_XSD")
        if schema_path is None:
            bundled_schema = Path(__file__).with_name("uniprot.xsd")
            if bundled_schema.exists():
                schema_path = str(bundled_schema)
        self._schema = xmlschema.XMLSchema(
            schema_path or "http://www.uniprot.org/docs/uniprot.xsd"
        )

    def _initialize_seq(self):
        """
        Initialize sequence data from UniProt Swiss-Prot database.

        Loads protein sequences from the downloaded Swiss-Prot FASTA file
        and creates a mapping of UniProt IDs to their sequences.
        """
        self._seq = {}
        fasta_path = _settings.POOCH.fetch("uniprot_sprot.fasta.gz")
        with gzip.open(fasta_path, "rt") as f:
            for record in tqdm(SeqIO.parse(f, "fasta"), desc="Loading sequences"):
                id = record.id.split("|")[1]
                self._seq[id] = record.seq

    def _initialize_genename_to_uniprot(self, organism="human"):
        """
        Initialize gene name to UniProt ID mapping for a specific organism.

        Args:
            organism (str): Organism name (default: "human")
                Must be one of the keys in organism_to_uniprot dictionary
        """
        uniprot_organism_id = organism_to_uniprot[organism]
        mapping_file = _settings.POOCH.fetch(f"{uniprot_organism_id}_idmapping.dat.gz")

        gene_names = []
        with gzip.open(mapping_file, "rt") as f:
            for line in f:
                if "Gene_Name" in line:
                    gene_names.append(line.strip().split("\t"))

        df = pd.DataFrame(gene_names, columns=["uniprot", "type", "gene_name"])
        df = df.groupby("gene_name").first()
        self._genename_to_uniprot = df.to_dict()["uniprot"]

    def _initialize_lddt(self):
        """Initialize pLDDT scores"""
        lddt_path = _settings.POOCH.fetch("9606.pLDDT.tdt.zip")
        pickle_path = get_protein_data_dir() / "9606.pLDDT.pickle"

        if pickle_path.exists():
            with Path.open(pickle_path, "rb") as f:
                self._lddt = pickle.load(f)
        else:
            self._lddt = {}
            with zipfile.ZipFile(lddt_path, "r") as f:
                for line in f.open("9606.pLDDT.tdt"):
                    id, score = line.strip().split(b"\t")
                    # byte to string
                    id = id.decode()
                    self._lddt[id] = np.array(score.split(b",")).astype(float)

            with Path.open(pickle_path, "wb") as f:
                pickle.dump(self._lddt, f)


# Create singleton instance
_protein_data = ProteinData()


# Public API functions
def get_seq_from_uniprot_id(uniprot_id: str):
    """
    Get protein sequence for a given UniProt ID.

    Args:
        uniprot_id (str): UniProt identifier (e.g., "P04637")

    Returns:
        Seq: Biopython Seq object containing the protein sequence

    Raises:
        KeyError: If the UniProt ID is not found in the database
    """
    return _protein_data.seq[uniprot_id]


def get_seq_from_gene_name(gene_name: str):
    """
    Get protein sequence for a given gene name.

    Args:
        gene_name (str): Official gene symbol (e.g., "TP53")

    Returns:
        Seq: Biopython Seq object containing the protein sequence

    Raises:
        KeyError: If the gene name is not found in the database
    """
    return _protein_data.seq[get_uniprot_from_gene_name(gene_name)]


def get_uniprot_from_gene_name(gene_name: str):
    """
    Get UniProt ID for a given gene name.

    Args:
        gene_name (str): Official gene symbol (e.g., "TP53")

    Returns:
        str: UniProt identifier

    Raises:
        KeyError: If the gene name is not found in the database
    """
    return _protein_data.genename_to_uniprot[gene_name]


def get_lddt_from_uniprot_id(uniprot_id: str):
    """
    Get AlphaFold pLDDT scores for a given UniProt ID.

    Args:
        uniprot_id (str): UniProt identifier (e.g., "P04637")

    Returns:
        numpy.ndarray: Array of pLDDT scores for each residue

    Raises:
        KeyError: If the UniProt ID is not found in the database
    """
    return _protein_data.lddt[uniprot_id]


def get_lddt_from_gene_name(gene_name: str):
    """
    Get AlphaFold pLDDT scores for a given gene name.

    Args:
        gene_name (str): Official gene symbol (e.g., "TP53")

    Returns:
        numpy.ndarray: Array of pLDDT scores for each residue

    Raises:
        KeyError: If the gene name is not found in the database
    """
    return _protein_data.lddt[get_uniprot_from_gene_name(gene_name)]


def _get_schema():
    """
    Get the UniProt XML schema.

    Returns:
        xmlschema.XMLSchema: Schema for parsing UniProt XML entries
    """
    return _protein_data.schema


def initialize_all():
    """
    Initialize all protein data.

    This function triggers the download and loading of all protein data files.
    It's useful when you want to preload all data at once instead of using
    lazy loading.
    """
    _protein_data.initialize_all()
