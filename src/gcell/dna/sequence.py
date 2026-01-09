"""
DNA sequence manipulation and analysis module.

This module provides classes for working with DNA sequences, including one-hot encoding,
mutation analysis, and various file format conversions. It supports operations like
reverse complement generation, sequence padding, and motif scanning.

Classes
-------
DNASequence: A class for manipulating individual DNA sequences
DNASequenceCollection: A class for working with collections of DNA sequences
"""

from collections.abc import Iterator
from pathlib import Path

import numcodecs
import numpy as np
import pandas as pd
import zarr
from Bio import SeqIO
from Bio.Seq import Seq
from numpy.typing import NDArray
from scipy.sparse import csr_matrix, save_npz, vstack
from tqdm import tqdm

from .motif import MotifCollection, print_results


class DNASequence(Seq):
    """
    A class for representing and manipulating DNA sequences.

    This class extends Seq with additional functionality for DNA sequence
    manipulation, including one-hot encoding, padding, and mutation operations.

    Parameters
    ----------
    seq : str
        The DNA sequence string
    header : str, optional
        Identifier or description for the sequence (default: "")

    """

    def __init__(self, seq: str, header: str = "") -> None:
        self.header: str = header
        self.seq: str = str(seq).upper()
        self._data: str = str(seq).upper()

        self.one_hot_encoding: dict[str, list[int]] = {
            "A": [1, 0, 0, 0],
            "C": [0, 1, 0, 0],
            "G": [0, 0, 1, 0],
            "T": [0, 0, 0, 1],
            "N": [0, 0, 0, 0],
        }

    def __repr__(self) -> str:
        return self.header

    def get_reverse_complement(self) -> str:
        """
        Generate the reverse complement of the DNA sequence.

        Returns
        -------
        str
            The reverse complement sequence
        """
        complement: dict[str, str] = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}
        return "".join([complement[base] for base in self.seq[::-1]])

    def padding(
        self, left: int = 0, right: int = 0, target_length: int = 0
    ) -> "DNASequence":
        """
        Pad the DNA sequence with N's.

        Parameters
        ----------
        left : int, optional
            Number of N's to add to the left (default: 0)
        right : int, optional
            Number of N's to add to the right (default: 0)
        target_length : int, optional
            Desired total length of sequence (default: 0)

        Returns
        -------
        DNASequence
            A new DNASequence object with padding applied
        """
        if target_length == 0:
            return DNASequence("N" * left + self.seq + "N" * right, self.header)
        elif target_length >= len(self.seq):
            return DNASequence(
                "N" * left + self.seq + "N" * (target_length - len(self.seq) - left),
                self.header,
            )
        elif target_length < len(self.seq):
            return DNASequence(
                self.seq[
                    (len(self.seq) - target_length) // 2 : (
                        len(self.seq) + target_length
                    )
                    // 2
                ],
                self.header,
            )

    def mutate(self, pos: int, alt: str) -> "DNASequence":
        """
        Create a new sequence with a mutation at the specified position.

        Parameters
        ----------
        pos : int
            Position to mutate (0-based indexing)
        alt : str
            Alternative base(s) to insert

        Returns
        -------
        DNASequence
            A new DNASequence object with the mutation applied
        """
        from Bio.Seq import MutableSeq

        if not isinstance(pos, int):
            pos = int(pos)
        if len(alt) == 1:
            seq = MutableSeq(self.seq)
            seq[pos] = alt
        else:
            seq = str(self.seq)
            seq = seq[0:pos] + alt + seq[pos + 1 :]
        return DNASequence(str(seq), self.header)

    # attribute to get one-hot encoding
    @property
    def one_hot(self) -> NDArray[np.int8]:
        """
        Convert the sequence to one-hot encoded format.

        Returns
        -------
        numpy.ndarray
            One-hot encoded representation of the sequence, shape (sequence_length, 4)
        """
        return (
            np.array([self.one_hot_encoding[base] for base in self.seq])
            .astype(np.int8)
            .reshape(-1, 4)
        )

    def save_zarr(
        self,
        zarr_file_path: str | Path,
        included_chromosomes: list[str] = [
            "chr1",
            "chr2",
            "chr3",
            "chr4",
            "chr5",
            "chr6",
            "chr7",
            "chr8",
            "chr9",
            "chr10",
            "chr11",
            "chr12",
            "chr13",
            "chr14",
            "chr15",
            "chr16",
            "chr17",
            "chr18",
            "chr19",
            "chr20",
            "chr21",
            "chr22",
            "chrX",
            "chrY",
        ],
    ) -> None:
        """
        Save the genome sequence data in Zarr format.

        Args:
            zarr_file_path (str): Path to the Zarr file containing genome data.
            included_chromosomes (list): List of chromosomes to be included in the Zarr file.
        """
        zarr_file = zarr.open_group(zarr_file_path, "w")
        for chr in tqdm(included_chromosomes):
            data = self.get_sequence(chr, 0, self.chrom_sizes[chr]).one_hot
            zarr_file.create_array(
                chr,
                data=data,
                chunks=(2000000, 4),
                dtype="i4",
                compressors=numcodecs.Blosc(cname="zstd", clevel=3, shuffle=2),
            )
        return


class DNASequenceCollection:
    """
    A collection of DNA sequences with batch processing capabilities.

    This class provides methods for working with multiple DNA sequences, including
    file I/O operations, motif scanning, and various data format conversions.

    Parameters
    ----------
    sequences : list
        List of Bio.SeqRecord objects or DNASequence objects


    """

    def __init__(self, sequences: list[SeqIO.SeqRecord]) -> None:
        self.sequences = sequences

    def __iter__(self) -> Iterator[DNASequence]:
        """Yield each sequence as a DNASequence object."""
        for seq in self.sequences:
            yield DNASequence(str(seq.seq), seq.id)

    @classmethod
    def from_fasta(cls, filename: str | Path) -> "DNASequenceCollection":
        """
        Create a DNASequenceCollection from a FASTA file.

        Parameters
        ----------
        filename : str or Path
            Path to the FASTA file

        Returns
        -------
        DNASequenceCollection
            New collection containing sequences from the FASTA file
        """
        return cls(list(SeqIO.parse(filename, "fasta")))

    def mutate(
        self, pos_list: list[int], alt_list: list[str]
    ) -> "DNASequenceCollection":
        """Mutate sequences at specified positions."""
        return DNASequenceCollection(
            [
                seq.mutate(pos, alt)
                for seq, pos, alt in zip(self.sequences, pos_list, alt_list)
            ]
        )

    def scan_motif(
        self, motifs: MotifCollection, non_negative: bool = True, raw: bool = False
    ) -> pd.DataFrame:
        """
        Scan sequences for specified motifs.

        Parameters
        ----------
        motifs : MotifCollection
            Collection of motifs to scan for
        non_negative : bool, optional
            If True, set negative scores to 0 (default: True)
        raw : bool, optional
            If True, return raw scanning results (default: False)

        Returns
        -------
        pandas.DataFrame
            Motif scanning results, either as raw data or summarized scores
        """
        seqs = self.sequences
        # initialize the output list
        output = []
        # scan each sequence and add the results to the output list
        headers = []
        lengths = []
        sequences = []
        for s in seqs:
            sequences.append(str(s.seq))
            headers.append(s.header)
            lengths.append(len(str(s.seq)))

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

        if raw is True:
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

        motif_c = pd.CategoricalDtype(categories=motifs.cluster_names, ordered=True)
        seq_c = pd.CategoricalDtype(categories=headers, ordered=True)

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

    def save_npz(self, filename: str | Path) -> None:
        """
        Save a DNASequenceCollection object as one-hot encoding in a sparse matrix in npz format,
        with sequence length information included in the filename
        """
        # create a list to store the sparse matrices
        sparse_matrices = []

        # loop over the sequences and create a sparse matrix for each one-hot encoding
        for seq in tqdm(self.sequences):
            # create the sparse matrix for the one-hot encoding
            sparse_matrix = csr_matrix(seq.one_hot)
            # add the sparse matrix to the list
            sparse_matrices.append(sparse_matrix)

        # concatenate the sparse matrices vertically
        sparse_matrix = vstack(sparse_matrices)

        # save the sparse matrix to a npz file with sequence length information in the filename
        save_npz(filename, sparse_matrix)

    def save_txt(self, filename: str | Path) -> None:
        """Save the DNASequenceCollection object as a text file"""
        with Path(filename).open("w") as f:
            for seq in self.sequences:
                f.write(seq.seq + "\n")

    def save_zarr(
        self,
        filename: str | Path,
        chunks: tuple[int, int, int] = (100, 2000, 4),
        target_length: int = 2000,
    ) -> None:
        """Save the one-hot encoding of a DNASequenceCollection object as a zarr file. Don't use sparse matrix, use compression"""
        # create a list to store the one-hot encoding
        one_hot = []

        # loop over the sequences and create a one-hot encoding for each sequence
        for seq in tqdm(self.sequences):
            # pad the sequence
            if len(seq.seq) != target_length:
                seq = seq.padding(left=0, target_length=target_length)
            # create the one-hot encoding
            one_hot.append(seq.one_hot)

        # concatenate the one-hot encoding vertically
        one_hot = np.stack(one_hot).astype(np.int8)

        # save the one-hot encoding to a zarr file
        zarr.save(filename, one_hot, chunks=chunks)

    def save_zarr_group(
        self,
        zarr_root: str | Path,
        key: str,
        chunks: tuple[int, int, int] = (100, 2000, 4),
        target_length: int = 2000,
    ) -> None:
        """
        Save sequences as one-hot encoded data in a Zarr group.

        Parameters
        ----------
        zarr_root : str
            Root directory for the Zarr storage
        key : str
            Key for storing the dataset within the Zarr group
        chunks : tuple, optional
            Chunk size for Zarr storage (default: (100, 2000, 4))
        target_length : int, optional
            Target sequence length for padding (default: 2000)

        Notes
        -----
        Sequences are padded to target_length if necessary and stored as
        one-hot encoded arrays with zstd compression.
        """
        # create a list to store the one-hot encoding
        one_hot = []

        # loop over the sequences and create a one-hot encoding for each sequence
        for seq in tqdm(self.sequences):
            # pad the sequence
            if len(seq.seq) != target_length:
                seq = seq.padding(left=0, target_length=target_length)
            # create the one-hot encoding
            one_hot.append(seq.one_hot)

        # concatenate the one-hot encoding vertically
        one_hot = np.stack(one_hot).astype(np.int8)

        # Initialize a zarr group/store at the specified root
        zarr_group = zarr.open_group(zarr_root, mode="a")

        # save the one-hot encoding to the specified key in the zarr group
        zarr_group.create_array(
            key,
            data=one_hot,
            chunks=chunks,
            dtype="i1",
            compressors=numcodecs.Blosc(cname="zstd", clevel=3),
        )
