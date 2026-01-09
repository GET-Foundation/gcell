"""
Cell type analysis module for genomic data.

This module provides classes and utilities for analyzing cell type-specific genomic data,
including gene expression, motif analysis, and causal relationships between genes.

Classes
-------
Celltype: Base class for cell type analysis
GETCellType: Extended cell type class with additional functionality
GETHydraCellType: Cell type class optimized for hydra-config-based model analysis
OneTSSJacobian: Class for analyzing Jacobian data for one TSS
OneGeneJacobian: Class for analyzing Jacobian data for one gene
GeneByMotif: Class for analyzing gene by motif relationships

The module supports both local file system and S3 storage backends.
"""

import logging
from concurrent.futures import ProcessPoolExecutor
from io import BytesIO
from pathlib import Path

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import seaborn as sns
import zarr
from plotly.subplots import make_subplots
from scipy.sparse import coo_matrix, csr_matrix
from scipy.stats import zscore
from tqdm import tqdm

from .._logging import get_logger
from ..dna.nr_motif_v1 import NrMotifV1
from ..rna.gene import TSS
from ..utils.causal_lib import get_subnet, plotly_networkx_digraph, preprocess_net
from ..utils.lingam import LiNGAM
from ..utils.s3 import (
    load_np_with_s3,
    load_npz_with_s3,
    load_zarr_with_s3,
    path_exists_with_s3,
)

motif = NrMotifV1.load_from_pickle()
motif_clusters = motif.cluster_names

# Add logger configuration
logger = get_logger(__name__)


def process_chunk(zarr_path, jacobian_group, chunk_idxs):
    """
    Process a subset of gene indices (chunk), computing the motif summary.
    Returns the computed chunk as a NumPy array.
    """
    z_local = zarr.open(zarr_path, mode="r")
    strand_local = z_local["strand"][:]
    jac0_local = z_local["jacobians"][jacobian_group]["0"]["input"]["region_motif"]
    jac1_local = z_local["jacobians"][jacobian_group]["1"]["input"]["region_motif"]
    input_local = z_local["input"]

    results = np.zeros((len(chunk_idxs), jac0_local.shape[2]), dtype=np.float32)

    for idx, i in enumerate(chunk_idxs):
        s = strand_local[i]
        jac = jac0_local[i] if s == 0 else jac1_local[i]
        inp = input_local[i]  # shape: (200, 283)
        results[idx] = (jac * inp).mean(axis=0)

    return chunk_idxs, results


def parallel_get_gene_by_motif(
    zarr_path: str,
    jacobian_group: str = "exp",
    overwrite: bool = True,
    max_workers: int = None,
    chunk_size: int = 64,
):
    """
    Given a GETHydraCellType zarr file, compute a gene-by-motif in parallel,
    storing the result in a 'gene_by_motif' dataset in the same Zarr store after collecting results.
    Each open tss of a gene is included in the gene_by_motif, so the shape
    of the output is (celltype.gene_annot.shape[0], n_features).
    """
    z = zarr.open(zarr_path, mode="a")

    if "gene_by_motif" in z and not overwrite:
        print("gene_by_motif already exists and overwrite=False. Doing nothing.")
        return

    n_genes = z["available_genes"].shape[0]
    n_features = z["input"].shape[2]

    # Prepare the output dataset
    if "gene_by_motif" in z:
        del z["gene_by_motif"]
    out_ds = z.create_array(
        "gene_by_motif", shape=(n_genes, n_features), dtype=np.float32, overwrite=True
    )

    gene_indices = np.arange(n_genes)
    results = np.zeros((n_genes, n_features), dtype=np.float32)

    futures = []
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        for start in range(0, n_genes, chunk_size):
            end = min(start + chunk_size, n_genes)
            chunk_inds = gene_indices[start:end]
            futures.append(
                executor.submit(process_chunk, zarr_path, jacobian_group, chunk_inds)
            )

        for f in tqdm(futures, desc="Computing gene_by_motif"):
            chunk_idxs, chunk_results = f.result()
            results[chunk_idxs, :] = chunk_results

    # Write the results to Zarr after all chunks are processed
    out_ds[:, :] = results
    print("Done. 'gene_by_motif' array has shape:", out_ds.shape)


class OneTSSJacobian:
    """Jacobian data for one TSS.

    Parameters
    ----------
    data : np.ndarray
        The Jacobian data (Peak x Motif x Classes/Strands).
    tss : TSS
        The TSS object.
    region : pd.DataFrame
        The region data.
    features : list
        The features.
    num_cls : int, optional
        The number of classes/strands. Defaults to 2.
    num_region_per_sample : int, optional
        The number of regions per sample. Defaults to 200.
    """

    def __init__(
        self,
        data: np.ndarray,
        tss: TSS,
        region: pd.DataFrame,
        features: list,
        num_cls=2,
        num_region_per_sample=200,
        num_features=283,
    ) -> None:
        # check if the data dimension is correct:
        # (num_cls, num_region_per_sample, num_features)
        assert (
            data.reshape(-1).shape[0] == num_cls * num_region_per_sample * num_features
        )
        self.TSS = tss
        data = data.reshape(num_cls, num_region_per_sample, num_features)[tss.strand]
        data_df = pd.DataFrame(data, columns=features)
        data_df = pd.concat(
            [region.reset_index(drop=True), data_df.reset_index(drop=True)],
            axis=1,
            ignore_index=True,
        )
        data_df.columns = region.columns.tolist() + list(features)
        self.data = data_df
        self.num_cls = num_cls
        self.features = features
        self.num_region_per_sample = num_region_per_sample
        self.num_features = num_features

    def __repr__(self) -> str:
        return f"""TSS: {self.TSS}
        Data shape: {self.data.shape}
        """

    # function for arbitrary transformation of the data

    def transform(self, func):
        """Transform the data."""
        self.data = func(self.data)
        return self.data

    def motif_summary(self, stats="absmean"):
        """Summarize the motif scores."""
        # assert stats in ['mean', 'max']
        motif_data = self.data[self.features]
        if stats == "mean":
            return motif_data.mean(axis=0)
        elif stats == "max":
            return motif_data.max(axis=0)
        elif stats == "absmean":
            return motif_data.abs().mean(axis=0)
        elif stats == "signed_absmean":
            return motif_data.abs().mean(axis=0) * np.sign(motif_data.mean(axis=0))
        # if stats is a function
        elif callable(stats):
            return motif_data.apply(stats, axis=0)

    def region_summary(self, stats="absmean"):
        """Summarize the motif scores."""
        data = self.data[self.features]
        if stats == "mean":
            region_data = data.mean(axis=1)
        elif stats == "max":
            region_data = data.max(axis=1)
        elif stats == "absmean":
            region_data = data.abs().mean(axis=1)
        elif stats == "l2_norm":
            region_data = np.linalg.norm(data, axis=1)
        # if stats is a function
        elif callable(stats):
            region_data = data.apply(stats, axis=1)
        data = self.data.iloc[:, 0:4]
        data["Score"] = region_data
        return data

    def summarize(self, axis="motif", stats="absmean"):
        """Summarize the data."""
        if axis == "motif":
            return self.motif_summary(stats)
        elif axis == "region":
            return self.region_summary(stats)

    def get_pairs_with_l2_cutoff(self, cutoff: float):
        """Get the pairs with L2 Norm cutoff."""
        l2_norm = np.linalg.norm(self.data, axis=1)
        if l2_norm == 0:
            return None
        v = self.data.values
        v[v**2 < cutoff] = 0
        v = csr_matrix(v)
        # Get the row, col, and value arrays from the csr_matrix
        rows, cols = v.nonzero()  # row are region idx, col are motif/feature idx
        values = v.data
        focus = self.num_region_per_sample // 2
        start_idx = self.TSS.peak_id - focus
        gene = self.TSS.gene_name
        # get a dataframe of {'Chromosome', 'Start', 'End', 'Gene', 'Strand', 'Start', 'Pred', 'Accessibility', 'Motif', 'Score'}
        df = self.peak_annot.iloc[
            start_idx : start_idx + self.num_region_per_sample
        ].copy()[["Chromosome", "Start", "End"]]
        df = df.iloc[rows]
        df["Motif"] = [self.data.columns[m] for m in cols]
        df["Score"] = values
        df["Gene"] = gene
        df["Strand"] = self.TSS.strand
        df["TSS"] = self.TSS.start
        df = df[
            ["Chromosome", "Start", "End", "Gene", "Strand", "TSS", "Motif", "Score"]
        ]
        return df


class OneGeneJacobian(OneTSSJacobian):
    """Jacobian data for one gene"""

    def __init__(
        self,
        gene_name: str,
        data: np.ndarray,
        region: pd.DataFrame,
        features: list,
        num_cls=2,
        num_region_per_sample=200,
        num_features=283,
    ) -> None:
        self.gene_name = gene_name
        data_df = pd.DataFrame(data, columns=features)
        data_df = pd.concat(
            [region.reset_index(drop=True), data_df.reset_index(drop=True)],
            axis=1,
            ignore_index=True,
        )
        data_df.columns = region.columns.tolist() + list(features)
        self.data = data_df
        self.num_cls = num_cls
        self.features = features
        self.num_region_per_sample = num_region_per_sample
        self.num_features = num_features

    def __repr__(self) -> str:
        return f"""Gene: {self.gene_name}
        """


class GeneByMotif:
    """
    Class for analyzing gene by motif relationships.

    This class handles the analysis of relationships between genes and motifs,
    including causal graph generation and correlation analysis.

    Parameters
    ----------
    celltype : str, optional
        Cell type name
    interpret_dir : str, optional
        Interpretation directory path
    jacob : pandas.DataFrame, optional
        Jacobian data
    s3_file_sys : S3FileSystem, optional
        S3 filesystem object
    zarr_data_path : str, optional
        Path to zarr data

    """

    def __init__(
        self,
        celltype=None,
        interpret_dir=None,
        jacob=None,
        s3_file_sys=None,
        zarr_data_path=None,
    ) -> None:
        self.celltype = celltype
        self.data = jacob
        self.interpret_dir = interpret_dir
        self.zarr_data_path = zarr_data_path
        self.s3_file_sys = s3_file_sys
        self._corr = None
        self._causal = None

    @staticmethod
    def _process_permutation(args):
        """Process a single permutation. The permutation is required because part of the
        LiNGAM algorithm involves approximate triangularization of the data matrix, which
        will be affected by the initial order of the columns. We try to permute the columns
        to randomize the order. The column name will be carried over after permutation.

        Parameters
        ----------
        args (tuple): Tuple containing (index, data, permute_columns, self)

        Returns
        -------
        tuple: (index, causal_graph)
        """
        i, data, permute_columns, instance = args
        permuted_data = (
            data.iloc[:, np.random.permutation(data.shape[1])]
            if permute_columns
            else data.copy()
        )
        causal_g = instance.create_causal_graph(permuted_data)
        return i, causal_g

    def get_causal(
        self,
        edgelist_file=None,
        permute_columns=True,
        n=3,
        overwrite=False,
        max_workers=None,
    ):
        """
        Generate causal graph from gene by motif data.

        This method creates a causal graph representing relationships between genes
        and motifs using the LiNGAM algorithm. Results can be cached in zarr format.

        Parameters
        ----------
        edgelist_file : str, optional
            Path to save/load edge list
        permute_columns : bool, optional
            Whether to permute columns, by default True
        n : int, optional
            Number of permutations, by default 3
        overwrite : bool, optional
            Whether to overwrite existing data, by default False
        max_workers : int, optional
            Maximum number of parallel workers

        Returns
        -------
        networkx.DiGraph
            Directed graph representing causal relationships
        """
        # Try loading from edgelist if provided
        if (
            edgelist_file is not None
            and path_exists_with_s3(edgelist_file, s3_file_sys=self.s3_file_sys)
            and not overwrite
        ):
            return nx.read_weighted_edgelist(edgelist_file, create_using=nx.DiGraph)

        # Determine which zarr path to use
        zarr_path = self.zarr_data_path
        if zarr_path is None:
            zarr_path = (
                Path(self.interpret_dir)
                / self.celltype
                / "allgenes"
                / f"{self.celltype}.zarr"
            )

        # Check if causal data exists in zarr
        if (
            path_exists_with_s3(zarr_path, s3_file_sys=self.s3_file_sys)
            and not overwrite
        ):
            try:
                return self.load_causal_from_zarr(zarr_path)
            except Exception as e:
                print(f"Failed to load from zarr: {str(e)}")

        # Calculate causal data if not found
        data = zscore(self.data, axis=0)
        data = pd.DataFrame(data, columns=self.data.columns, index=self.data.index)
        zarr_data = load_zarr_with_s3(zarr_path, mode="a", s3_file_sys=self.s3_file_sys)

        # Calculate permutations in parallel
        pending_indices = [
            i for i in range(n) if f"causal_{i}" not in zarr_data or overwrite
        ]
        if pending_indices:
            with ProcessPoolExecutor(max_workers=max_workers) as executor:
                # Create arguments for each permutation - include self instance
                args_list = [(i, data, permute_columns, self) for i in pending_indices]

                # Process permutations in parallel and save results
                for i, causal_g in tqdm(
                    executor.map(self._process_permutation, args_list),
                    total=len(pending_indices),
                    desc="Processing permutations",
                ):
                    self.save_causal_to_zarr(zarr_data, causal_g, i)

        # Compute and save average
        average_causal = self.compute_average_causal(zarr_data, n)
        zarr_data.create_array("causal", data=average_causal, dtype="float32", overwrite=True)

        # Create final graph
        causal_g = nx.from_numpy_array(average_causal, create_using=nx.DiGraph)
        causal_g = nx.relabel_nodes(
            causal_g, dict(zip(range(len(self.data.columns)), self.data.columns))
        )

        # Save edgelist if path provided
        if edgelist_file:
            nx.write_weighted_edgelist(causal_g, edgelist_file)

        return causal_g

    def create_causal_graph(self, data):
        """
        Create a causal graph from the data by running LiNGAM.

        Parameters
        ----------
        data : pd.DataFrame
            The data to create a causal graph from.

        Returns
        -------
        networkx.DiGraph
            The causal graph.
        """
        model = LiNGAM()
        output = model.predict(data)
        causal_g = preprocess_net(
            output.copy(), remove_nodes=False, detect_communities=False
        )
        causal_g_numpy = nx.to_numpy_array(
            causal_g, dtype="float32", nodelist=self.data.columns
        )
        causal_g = nx.from_numpy_array(causal_g_numpy, create_using=nx.DiGraph)
        causal_g = nx.relabel_nodes(
            causal_g, dict(zip(range(len(self.data.columns)), self.data.columns))
        )
        return causal_g

    def load_causal_from_zarr(self, zarr_data_path):
        """Load causal data from zarr file.

        Supports both old and new zarr structures.
        """
        zarr_data = load_zarr_with_s3(
            file_path=zarr_data_path,
            mode="r",  # Changed to read-only mode since we're just loading
            s3_file_sys=self.s3_file_sys,
        )

        # Try to find causal data in zarr structure
        if "causal" in zarr_data:
            causal_array = zarr_data["causal"][:]
        elif "gene_by_motif_causal" in zarr_data:  # Alternative location
            causal_array = zarr_data["gene_by_motif_causal"][:]
        else:
            raise KeyError("No causal data found in zarr file")

        causal_g = nx.from_numpy_array(causal_array, create_using=nx.DiGraph)
        causal_g = nx.relabel_nodes(
            causal_g, dict(zip(range(len(self.data.columns)), self.data.columns))
        )
        return causal_g

    def save_causal_to_zarr(self, zarr_data, causal_g, index=None):
        """Save causal data to zarr file.

        Args:
            zarr_data: Zarr group to save to
            causal_g: Causal graph to save
            index: If provided, saves as causal_{index}, otherwise saves as causal
        """
        causal_g_numpy = nx.to_numpy_array(
            causal_g, dtype="float32", nodelist=self.data.columns
        )

        array_name = f"causal_{index}" if index is not None else "causal"
        zarr_data.create_array(array_name, data=causal_g_numpy, dtype="float32", overwrite=True)

    def compute_average_causal(self, zarr_data, n):
        """Compute the average causal graph across n permutations."""
        causal_arrays = [zarr_data[f"causal_{i}"] for i in range(n)]
        average_causal = np.mean(causal_arrays, axis=0)
        return average_causal

    def set_diagnal_to_zero(self, df: pd.DataFrame):
        for i in range(df.shape[0]):
            df.iloc[i, i] = 0
        return df


class Celltype:
    """
    Base class for cell type analysis.

    This class provides core functionality for analyzing cell type-specific genomic data,
    including gene expression, motif analysis, and regulatory relationships.

    The data-loading logic is now deprecated and only used for the demo data.

    For your own analysis, you should use GETHydraCellType.

    Parameters
    ----------
    features
        Array of feature names/identifiers
    num_region_per_sample
        Number of regions per sample
    celltype
        Name/identifier of the cell type
    data_dir
        Directory containing data files, by default "../pretrain_human_bingren_shendure_apr2023"
    interpret_dir
        Directory for interpretation results, by default "Interpretation"
    assets_dir
        Directory for assets/resources, by default "assets"
    input
        Whether to load input data, by default False
    jacob
        Whether to load Jacobian data, by default False
    embed
        Whether to load embedding data, by default False
    num_cls
        Number of classes, by default 2
    s3_file_sys
        S3 filesystem object for remote storage, by default None
    """

    def __init__(
        self,
        features: np.ndarray,
        num_region_per_sample: int,
        celltype: str,
        data_dir="../pretrain_human_bingren_shendure_apr2023",
        interpret_dir="Interpretation",
        assets_dir="assets/",
        input: bool = False,
        jacob: bool = False,
        embed: bool = False,
        num_cls: int = 2,
        s3_file_sys=None,
    ):
        self._gene_by_motif = None
        self.tss_accessibility = None
        self.input = None
        self.input_all = None
        self.embed = None
        self.jacobs = None
        self.preds = None
        self.obs = None
        self._zarr_data = None
        self.celltype_name = celltype
        self.celltype = celltype
        self.data_dir = data_dir
        self.interpret_dir = interpret_dir
        self.assets_dir = assets_dir
        self.features = features
        self.num_features = features.shape[0]
        self.num_region_per_sample = num_region_per_sample
        self.focus = num_region_per_sample // 2
        self.num_cls = num_cls
        self.s3_file_sys = s3_file_sys
        self.interpret_cell_dir = Path(self.interpret_dir) / celltype / "allgenes"
        print("loading from interpretation dir", self.interpret_cell_dir)
        self.gene_feather_path = f"{self.data_dir}{celltype}.exp.feather"
        # Load gencode_hg38 from feather file.
        # This is only used for backwards compatibility.
        self.gencode_hg38 = pd.read_feather(
            "https://zenodo.org/records/14635090/files/gencode.v40.hg38.feather?download=1"
        )
        self.gencode_hg38["Strand"] = self.gencode_hg38["Strand"].apply(
            lambda x: 0 if x == "+" else 1
        )
        self.gene2strand = self.gencode_hg38.set_index("gene_name").Strand.to_dict()
        if path_exists_with_s3(
            self.interpret_cell_dir / f"{self.celltype}.zarr",
            s3_file_sys=self.s3_file_sys,
        ):
            self._zarr_data = load_zarr_with_s3(
                self.interpret_cell_dir / f"{self.celltype}.zarr",
                mode="a",
                s3_file_sys=self.s3_file_sys,
            )
            self.genelist = (
                self._zarr_data["avaliable_genes"]
                if "avaliable_genes" in self._zarr_data
                else self._zarr_data["gene_names"]
            )
        else:
            self.genelist = load_np_with_s3(
                file_path=self.interpret_cell_dir / "avaliable_genes.npy",
                s3_file_sys=self.s3_file_sys,
            )

        if not path_exists_with_s3(
            self.gene_feather_path, s3_file_sys=self.s3_file_sys
        ):
            self.gene_feather_path = (
                f"{self.interpret_dir}/{celltype}.gene_idx_dict.feather"
            )
        self.peak_annot = pd.read_csv(
            self.data_dir + celltype + ".csv", sep=","
        ).rename(columns={"Unnamed: 0": "index"})
        self.gene_annot = self.load_gene_annot()
        tss_idx = self.gene_annot.level_0.values
        self.tss_idx = tss_idx
        if input:
            self.input = load_npz_with_s3(
                self.data_dir + celltype + ".watac.npz", s3_file_sys=self.s3_file_sys
            )[tss_idx]
            self.input_all = load_npz_with_s3(
                self.data_dir + celltype + ".watac.npz", s3_file_sys=self.s3_file_sys
            )
            self.tss_accessibility = self.input[:, self.num_features - 1]

        self.gene_annot["Strand"] = self.gene_annot["gene_name"].apply(
            lambda x: self.gene2strand.get(x, 0)
        )
        self.tss_strand = self.gene_annot.Strand.astype(int).values
        self.tss_start = self.peak_annot.iloc[tss_idx].Start.values

        if hasattr(self, "_zarr_data"):
            import time

            if jacob:
                start_time = time.time()
                self.jacobs = self._zarr_data["jacobians"]
                # print time with 2 decimals
                print(f"loaded jacobians in {time.time()-start_time:.2f} seconds")
            start_time = time.time()
            self.preds = np.array(self._zarr_data["preds"])
            self.preds = np.array(
                [
                    self.preds[i].reshape(self.num_region_per_sample, self.num_cls)[
                        self.focus, j
                    ]
                    for i, j in enumerate(self.tss_strand)
                ]
            )
            self.obs = np.array(self._zarr_data["obs"])
            self.obs = np.array(
                [
                    self.obs[i].reshape(self.num_region_per_sample, self.num_cls)[
                        self.focus, j
                    ]
                    for i, j in enumerate(self.tss_strand)
                ]
            )
            print(f"loaded preds and obs in {time.time()-start_time:.2f} seconds")
            if embed:
                start_time = time.time()
                self.embed = self._zarr_data["embeds_0"]
                print(f"loaded embeds in {time.time()-start_time:.2f} seconds")

        else:
            if embed:
                if path_exists_with_s3(
                    self.interpret_cell_dir / "embeds_0.npy",
                    s3_file_sys=self.s3_file_sys,
                ):
                    self.embed = load_np_with_s3(
                        self.interpret_cell_dir / "embeds_0.npy",
                        s3_file_sys=self.s3_file_sys,
                    )
                else:
                    raise ValueError("embeds_0.npy not found")

            if jacob:
                # check if os.path.join(self.interpret_cell_dir, "jacobians.zarr") exists, if not save the matrix to zarr file
                if path_exists_with_s3(
                    self.interpret_cell_dir / "jacobians.zarr",
                    s3_file_sys=self.s3_file_sys,
                ):
                    # load from zarr file
                    self.jacobs = load_zarr_with_s3(
                        self.interpret_cell_dir / "jacobians.zarr",
                        mode="r",
                        s3_file_sys=self.s3_file_sys,
                    )
                else:
                    jacob_npz = coo_matrix(
                        load_npz_with_s3(
                            self.interpret_cell_dir / "jacobians.npz",
                            s3_file_sys=self.s3_file_sys,
                        )
                    )
                    z = zarr.zeros(
                        shape=jacob_npz.shape,
                        chunks=(100, jacob_npz.shape[1]),
                        dtype=jacob_npz.dtype,
                    )
                    z.set_coordinate_selection(
                        (jacob_npz.row, jacob_npz.col), jacob_npz.data
                    )
                    # save to zarr file
                    zarr.save(self.interpret_cell_dir / "jacobians.zarr", z)
                    self.jacobs = load_zarr_with_s3(
                        self.interpret_cell_dir / "jacobians.zarr",
                        mode="r",
                        s3_file_sys=self.s3_file_sys,
                    )

            self.preds = load_npz_with_s3(
                self.interpret_cell_dir / "preds.npz",
                s3_file_sys=self.s3_file_sys,
            )
            self.preds = np.array(
                [
                    self.preds[i]
                    .toarray()
                    .reshape(self.num_region_per_sample, self.num_cls)[self.focus, j]
                    for i, j in enumerate(self.tss_strand)
                ]
            )
            self.obs = load_npz_with_s3(
                self.interpret_cell_dir / "obs.npz",
                s3_file_sys=self.s3_file_sys,
            )
            self.obs = np.array(
                [
                    self.obs[i]
                    .toarray()
                    .reshape(self.num_region_per_sample, self.num_cls)[self.focus, j]
                    for i, j in enumerate(self.tss_strand)
                ]
            )

        self.gene_annot["pred"] = self.preds
        self.gene_annot["obs"] = self.obs
        if input:
            self.gene_annot["accessibility"] = self.tss_accessibility.toarray().reshape(
                -1
            )
        else:
            self.gene_annot["accessibility"] = 1
        self.gene_annot["Chromosome"] = self.peak_annot.iloc[tss_idx].Chromosome.values
        self.gene_annot["Start"] = self.tss_start
        self._gene_by_motif = None

    def __repr__(self) -> str:
        return f"""Celltype: {self.celltype}
        Data dir: {self.data_dir}
        Interpretation dir: {self.interpret_dir}
        Number of regions per sample: {self.num_region_per_sample}
        Number of features: {self.num_features}
        Number of genes: {self.gene_annot.gene_name.unique().shape[0]}
        Number of TSS: {self.tss_idx.shape[0]}
        Number of peaks: {self.peak_annot.shape[0]}
        """

    def load_gene_annot(self):
        """Load gene annotations from feather file.

        If the feather file does not exist, construct the gene annotation from gencode.
        Note that this is largely for backwards compatibility.
        """
        if not path_exists_with_s3(
            self.gene_feather_path, s3_file_sys=self.s3_file_sys
        ):
            print(
                "Gene exp feather not found. Constructing gene annotation from gencode..."
            )
            # construct gene annotation from gencode
            from pyranges import PyRanges as pr

            atac = pr(self.peak_annot, int64=True)
            # join the ATAC-seq data with the RNA-seq data
            exp = atac.join(
                pr(self.gencode_hg38, int64=True).extend(300), how="left"
            ).as_df()
            self.gene2strand["-1"] = -1
            exp["Strand"] = exp["gene_name"].apply(lambda x: self.gene2strand[x])
            if self.s3_file_sys is None:
                exp.reset_index(drop=True).to_feather(
                    f"{self.data_dir}{self.celltype}.exp.feather"
                )
            self.gene_feather_path = f"{self.data_dir}{self.celltype}.exp.feather"
            gene_annot = exp
        else:
            print("Gene exp feather found. Loading...")
            gene_annot = pd.read_feather(self.gene_feather_path)
            gene_annot["Strand"] = gene_annot["Strand"].apply(
                lambda x: 0 if x == "+" else 1
            )
        if self.gene_feather_path.endswith(".exp.feather"):
            gene_annot = (
                gene_annot.groupby(["gene_name", "Strand"])["index"]
                .unique()
                .reset_index()
                .dropna()
                .query('gene_name!="-1"')
                .rename(columns={"index": "level_0"})
            )
        gene_annot = gene_annot.explode("level_0").reset_index(drop=True)
        if isinstance(self.genelist, zarr.Array):
            self.genelist = self.genelist[:]
        gene_annot = gene_annot.iloc[self.genelist].reset_index(drop=True)
        return gene_annot

    def get_gene_idx(self, gene_name: str):
        """
        Get the index of a gene in the gene list.

        Parameters
        ----------
        gene_name
            Name of the gene

        Returns
        -------
        numpy.ndarray
            Array of indices where the gene appears
        """
        return self.gene_annot[self.gene_annot["gene_name"] == gene_name].index.values

    def get_tss_idx(self, gene_name: str):
        """Given a gene name, get the TSS index in the peak annotation.

        Parameters
        ----------
        gene_name
            The name of the gene.

        Returns
        -------
        numpy.ndarray
            The index of the TSS in the peak annotation.
        """
        return self.tss_idx[self.get_gene_idx(gene_name)]

    def get_gene_jacobian(
        self, gene_name: str, multiply_input: bool = True
    ) -> list[OneTSSJacobian]:
        """Get the jacobian of a gene. Each gene can have multiple TSSs, so this function returns a list of TSSJacobian objects.
        Each OneTSSJacobian contains a Peak x Motif jacobian score for a given TSS.

        Parameters
        ----------
        gene_name
            The name of the gene.
        multiply_input
            If True, the input is multiplied. Defaults to True.

        Returns
        -------
        list
            A list of TSSJacobian objects.
        """
        gene_idx = self.get_gene_idx(gene_name)
        gene_chr = self.get_gene_chromosome(gene_name)
        jacobs = []
        for i in gene_idx:
            # get a TSS object
            tss = TSS(
                gene_name,
                self.tss_idx[i],
                gene_chr,
                self.tss_start[i],
                self.tss_strand[i],
            )
            jacobs.append(self.get_tss_jacobian(self.jacobs[i], tss, multiply_input))
        return jacobs

    def get_input_data(self, peak_id=None, focus=None, start=None, end=None):
        """Get input data from self.input_all using a slice

        Parameters
        ----------
        peak_id
            The peak id to get the input data from.
        focus
            The focus to get the input data from. Usually the middle of the region.
        start
            The start to get the input data from.
        end
            The end to get the input data from.
        """
        # assert if all are None
        assert not all([peak_id, focus, start, end])
        if start is None:
            start = peak_id - focus
            end = peak_id + focus
        return self.input_all[start:end].toarray().reshape(-1, self.num_features)

    def get_tss_jacobian(self, jacob: np.ndarray, tss: TSS, multiply_input=True):
        """
        Get the jacobian of a TSS.

        Parameters
        ----------
        jacob
            The jacobian to be processed.
        tss
            The TSS object.
        multiply_input
            If True, the input is multiplied. Defaults to True.
        """
        jacob = jacob.reshape(-1, self.num_region_per_sample, self.num_features)
        if multiply_input:
            input = self.get_input_data(peak_id=tss.peak_id, focus=self.focus)
            jacob = jacob * input
        region = tss.get_sample_from_peak(self.peak_annot, self.focus)
        tss_jacob = OneTSSJacobian(jacob, tss, region, self.features, self.num_cls)
        return tss_jacob

    def gene_jacobian_summary(
        self, gene, axis="motif", multiply_input=True, stats="absmean"
    ):
        """
        Summarizes the Jacobian for a given gene.

        This function calculates the Jacobian for a given gene and summarizes it based on the specified axis.
        If the axis is "motif", it concatenates the Jacobian summaries along axis 1 and then sums them.
        If the axis is "region", it concatenates the Jacobian summaries along axis 0, groups them by index, chromosome, and start,
        and then sums the scores.

        Parameters
        ----------
        gene
            The gene for which the Jacobian is to be calculated.
        axis
            The axis along which to summarize the Jacobian. Defaults to "motif".
        multiply_input
            If True, the input is multiplied. Defaults to True.

        Returns
        -------
        pd.DataFrame: A DataFrame containing the summarized Jacobian.
        """
        gene_jacobs = self.get_gene_jacobian(gene, multiply_input)
        if axis == "motif":
            return pd.concat(
                [jac.summarize(axis, stats=stats) for jac in gene_jacobs], axis=1
            ).sum(axis=1)
        elif axis == "region":
            # concat in axis 0 and aggregate overlapping regions by sum the score and divided by number of tss
            return (
                pd.concat([j.summarize(axis, stats=stats) for j in gene_jacobs])
                .groupby(["index", "Chromosome", "Start", "End"])
                .Score.sum()
                .reset_index()
            )

    @property
    def gene_by_motif(self):
        """
        This method retrieves gene data by motif. It first checks if a zarr file exists for the cell type.
        If it does, it opens the zarr file and checks if 'gene_by_motif' is in the keys of the zarr data.
        If 'gene_by_motif' is found, it loads the data into a pandas DataFrame. If not, it computes the
        jacobian for each gene and saves it to the zarr file.

        If a zarr file does not exist, it checks if a feather file exists. If it does, it loads the data
        into a pandas DataFrame. If not, it computes the jacobian for each gene and saves it to a feather file.

        Finally, if the gene_by_motif data is a pandas DataFrame, it is converted to a GeneByMotif object.
        If a zarr file exists, it checks if 'gene_by_motif_corr' is in the keys of the zarr data. If it is,
        it loads the correlation data into the GeneByMotif object. If not, it computes the correlation and
        saves it to the zarr file.

        """
        if self._gene_by_motif is None:
            self._gene_by_motif = self.get_gene_by_motif()
        return self._gene_by_motif

    @gene_by_motif.setter
    def gene_by_motif(self, value):
        self._gene_by_motif = value

    def get_gene_by_motif(self, overwrite: bool = False):
        """
        This method retrieves gene data by motif. It first checks if a zarr file exists for the cell type.
        If it does, it opens the zarr file and checks if 'gene_by_motif' is in the keys of the zarr data.
        If 'gene_by_motif' is found, it loads the data into a pandas DataFrame. If not, it computes the
        jacobian for each gene and saves it to the zarr file.

        If a zarr file does not exist, it checks if a feather file exists. If it does, it loads the data
        into a pandas DataFrame. If not, it computes the jacobian for each gene and saves it to a feather file.

        Finally, if the gene_by_motif data is a pandas DataFrame, it is converted to a GeneByMotif object.
        If a zarr file exists, it checks if 'gene_by_motif_corr' is in the keys of the zarr data. If it is,
        it loads the correlation data into the GeneByMotif object. If not, it computes the correlation and
        saves it to the zarr file.

        Parameters
        ----------
        overwrite
            If True, overwrite the existing data. Defaults to False.

        Returns
        -------
        GeneByMotif: An object that contains the gene data by motif and the correlation data.
        """
        if path_exists_with_s3(
            self.interpret_cell_dir / f"{self.celltype}.zarr",
            s3_file_sys=self.s3_file_sys,
        ):
            self._zarr_data = load_zarr_with_s3(
                self.interpret_cell_dir / f"{self.celltype}.zarr",
                mode="a",
                s3_file_sys=self.s3_file_sys,
            )
            if "gene_by_motif" in self._zarr_data:
                self._gene_by_motif = pd.DataFrame(
                    self._zarr_data["gene_by_motif"], columns=self.features
                )
            else:
                jacobs = []
                for g in tqdm(self.gene_annot.gene_name.unique()):
                    for j in self.get_gene_jacobian(g, multiply_input=True):
                        jacobs.append(j.motif_summary().T)
                jacobs_df = pd.concat(jacobs, axis=1).T
                # save to zarr
                self._zarr_data.create_array(
                    "gene_by_motif", data=jacobs_df.values, dtype="float32"
                )
                self._gene_by_motif = jacobs_df
        elif path_exists_with_s3(
            f"{self.interpret_dir}/{self.celltype}_gene_by_motif.feather",
            s3_file_sys=self.s3_file_sys,
        ):
            self._gene_by_motif = pd.read_feather(
                f"{self.interpret_dir}/{self.celltype}_gene_by_motif.feather"
            ).set_index("index")
        else:
            jacobs = []
            for g in tqdm(self.gene_annot.gene_name.unique()):
                for j in self.get_gene_jacobian(g, multiply_input=True):
                    jacobs.append(j.motif_summary().T)
            jacobs_df = pd.concat(jacobs, axis=1).T
            jacobs_df.reset_index().to_feather(
                f"{self.interpret_dir}/{self.celltype}_gene_by_motif.feather"
            )
            self._gene_by_motif = jacobs_df

        if isinstance(self._gene_by_motif, pd.DataFrame):
            self._gene_by_motif = GeneByMotif(
                self.celltype, self.interpret_dir, self.gene_by_motif, self.s3_file_sys
            )
            if path_exists_with_s3(
                Path(self.interpret_cell_dir) / f"{self.celltype}.zarr",
                s3_file_sys=self.s3_file_sys,
            ):
                if "gene_by_motif_corr" in self._zarr_data:
                    self._gene_by_motif.corr = pd.DataFrame(
                        self._zarr_data["gene_by_motif_corr"],
                        columns=self.features,
                        index=self.features,
                    )

                else:
                    # compute corr and save to zarr also
                    self._zarr_data.create_array(
                        "gene_by_motif_corr",
                        data=self._gene_by_motif.get_corr().values,
                        dtype="float32",
                    )
        self._gene_by_motif.data.set_index(self.gene_annot.gene_name, inplace=True)
        return self._gene_by_motif

    def get_tf_pathway(
        self,
        tf,
        gp=None,
        quantile_cutoff=0.9,
        exp_cutoff=0,
        filter_str="term_size<1000 & term_size>500",
        significance_threshold_method="g_SCS",
    ):
        """
        This function retrieves the pathway for a given transcription factor (tf) using g:Profiler.

        Parameters
        ----------
        tf
            The transcription factor to get the pathway for.
        gp
            An instance of the GProfiler class. If None, a new instance will be created. Defaults to None.
        quantile_cutoff
            The quantile cutoff to use when selecting genes. Defaults to 0.9.
        exp_cutoff
            The expression cutoff to use when querying genes. Defaults to 0.
        filter_str
            The filter string to use when querying the g:Profiler results. Defaults to 'term_size<1000 & term_size>500'.
        significance_threshold_method
            The method to use for determining the significance threshold in g:Profiler. Defaults to 'g_SCS'.

        Returns
        -------
        tuple: A tuple containing the filtered g:Profiler results and the unique genes in the pathways.
        """
        self.get_gene_by_motif()
        if gp is None:
            from gprofiler import GProfiler

            gp = GProfiler(return_dataframe=True)

        selected_index = self.gene_by_motif.data[tf] > self.gene_by_motif.data[
            tf
        ].quantile(quantile_cutoff)
        selected_genes = self.gene_by_motif.data.index[selected_index]
        gene_list = (
            self.gene_annot.query("gene_name.isin(@selected_genes)")
            .query("pred>@exp_cutoff")
            .gene_name.unique()
        )
        go = gp.profile(
            organism="hsapiens",
            query=list(gene_list),
            user_threshold=0.05,
            no_evidences=False,
            significance_threshold_method=significance_threshold_method,
        )
        go_filtered = go.query(filter_str)
        pathway_genes = np.unique(np.concatenate(go_filtered.intersections.values))
        return go_filtered, pathway_genes

    def get_highest_exp_genes(self, genes: list):
        """
        This code takes in a list of genes and returns the gene with the highest expression value.
        """
        return (
            self.gene_annot.query("gene_name.isin(@genes)")
            .sort_values("pred", ascending=False)
            .head(1)
            .gene_name.values[0]
        )

    def get_genes_exp(self, genes: list):
        """
        Get the expression of a list of genes.
        """
        return self.gene_annot.query("gene_name.isin(@genes)").sort_values(
            "pred", ascending=False
        )

    def get_tf_exp_str(self, motif, m):
        """
        This method generates a formatted string of gene names and their corresponding predicted expression values.
        The expression values are averaged and sorted in descending order.

        Parameters
        ----------
        motif
            The motif object containing the cluster gene list.
        m
            The index to access the specific cluster gene list from the motif object.

        Returns
        -------
        str: A string representation of gene names and their corresponding predicted expression values.
            The string is formatted as 'gene_name\tpred', where 'pred' is a 2 digit floating point number.
            Each gene name and its corresponding predicted expression value are separated by a '<br />'.

        Example:
        'gene1\t1.23<br />gene2\t0.56<br />gene3\t0.45'
        """
        if m not in motif.cluster_gene_list:
            return m
        motif_cluster_genes = motif.cluster_gene_list[m]
        motif_cluster_genes_exp = (
            self.get_genes_exp(motif_cluster_genes)
            .groupby("gene_name")
            .pred.mean()
            .sort_values(ascending=False)
        )
        # turn in to a formated table in one string: gene_name\tpred, 2 digit floating point
        return "<br />".join(
            [
                f"{gene_name}\t{pred:.2f}"
                for gene_name, pred in motif_cluster_genes_exp.items()
            ]
        )

    def get_tf_exp_mean(self, motif, m):
        """
        Calculate the mean expression of transcription factors (TFs) for a given motif and cluster.

        Parameters:
        motif (Motif): The motif object containing information about the motif and associated genes.
        m (int): The cluster index for which the mean TF expression is to be calculated.

        Returns:
        float: The mean expression of TFs for the given motif and cluster.
        """
        if m not in motif.cluster_gene_list:
            return np.nan
        motif_cluster_genes = motif.cluster_gene_list[m]
        motif_cluster_genes_exp = (
            self.get_genes_exp(motif_cluster_genes).groupby("gene_name").pred.mean()
        )
        return motif_cluster_genes_exp.mean()

    def get_gene_pred(self, gene_name: str):
        """Get the prediction of a gene."""
        return self.preds[self.get_gene_idx(gene_name)]

    def get_gene_obs(self, gene_name: str):
        """Get the observed value of a gene."""
        return self.obs[self.get_gene_idx(gene_name)]

    def get_gene_annot(self, gene_name: str):
        """Get the gene annotation of a gene."""
        return self.gene_annot[self.gene_annot["gene_name"] == gene_name]

    def get_gene_accessibility(self, gene_name: str):
        """Get the accessibility of a gene."""
        gene_idx = self.get_gene_idx(gene_name)
        return (
            self.tss_accessibility[gene_idx]
            if hasattr(self, "tss_accessibility")
            else None
        )

    def get_gene_strand(self, gene_name: str):
        """Get the strand of a gene."""
        return self.tss_strand[self.get_gene_idx(gene_name)].unique()[0]

    def get_gene_tss_start(self, gene_name: str):
        """Get the start position of a gene."""
        return self.tss_start[self.get_gene_idx(gene_name)]

    def get_gene_tss(self, gene_name: str):
        """Get the TSS objects of a gene."""
        gene_idx = self.get_gene_idx(gene_name)
        gene_chr = self.get_gene_chromosome(gene_name)
        return [
            TSS(
                gene_name,
                self.tss_idx[i],
                gene_chr,
                self.tss_start[i],
                self.tss_strand[i],
            )
            for i in gene_idx
        ]

    def get_gene_chromosome(self, gene_name: str):
        """Get the chromosome of a gene."""
        if gene_name not in self.gene_annot["gene_name"].tolist():
            logger.error(
                "This gene's promoter is not in this cell type's open chromatin."
            )
        return self.gene_annot[self.gene_annot["gene_name"] == gene_name][
            "Chromosome"
        ].values[0]

    def get_gene_jacobian_summary(self, gene_name: str, axis="motif"):
        """Get the jacobian summary of a gene."""
        gene_jacobs = self.get_gene_jacobian(gene_name)
        if axis == "motif":
            return pd.concat([j.summarize(axis) for j in gene_jacobs], axis=1).sum(
                axis=1
            )
        elif axis == "region":
            # concat in axis 0 and aggregate overlapping regions by sum the score and divided by number of tss
            return (
                pd.concat([j.summarize(axis="region") for j in gene_jacobs])
                .groupby(["index", "Chromosome", "Start", "End"])
                .Score.sum()
                .reset_index()
            )

    def plot_gene_motifs(self, gene, motif, overwrite: bool = False):
        r = self.get_gene_jacobian_summary(gene, "motif")
        m = r.sort_values(ascending=False).head(10).index.values
        fig, ax = plt.subplots(2, 5, figsize=(10, 4), sharex=False, sharey=False)
        for i, m_i in enumerate(m):
            if (
                not path_exists_with_s3(
                    file_path=f'{self.assets_dir}{m_i.replace("/", "_")}.png',
                    s3_file_sys=self.s3_file_sys,
                )
                or overwrite
            ):
                motif.get_motif_cluster_by_name(m_i).seed_motif.plot_logo(
                    filename=f'{self.assets_dir}{m_i.replace("/", "_")}.png',
                    logo_title="",
                    size="medium",
                    ic_scale=True,
                )
            # show logo in ax[i] from the png file
            if self.s3_file_sys:
                with self.s3_file_sys.open(
                    f'{self.assets_dir}{m_i.replace("/", "_")}.png', "rb"
                ) as f:
                    img = plt.imread(BytesIO(f.read()))
            else:
                img = plt.imread(f'{self.assets_dir}{m_i.replace("/", "_")}.png')
            ax[i // 5][i % 5].imshow(img)
            ax[i // 5][i % 5].axis("off")
            # add title to highest expressed gene
            if m_i in motif.cluster_gene_list:
                motif_cluster_genes = motif.cluster_gene_list[m_i]
            try:
                ax[i // 5][i % 5].set_title(
                    f"{m_i}:{self.get_highest_exp_genes(motif_cluster_genes)}"
                )
            except Exception as e:
                ax[i // 5][i % 5].set_title(f"{m_i}")

        return fig, ax

    def plotly_motif_subnet(self, motif, m, type="neighbors", threshold="auto"):
        """
        Plots a subnet of motifs.

        This function generates a subnet of motifs based on the given parameters and plots it using plotly.
        The subnet is preprocessed and the TF expression string and mean are calculated for each motif in the cluster gene list.

        Parameters:
        motif (Motif): The motif object to plot.
        m (int): The motif index to plot.
        type (str, optional): The type of subnet to generate. Can be 'neighbors', 'parents', or 'children'. Defaults to 'neighbors'.
        threshold (str or float, optional): The threshold for preprocessing the network. Can be 'auto' or a float. If 'auto', the threshold is set to the 70th percentile of the absolute weight values. Defaults to 'auto'.

        Returns:
        plotly.graph_objs._figure.Figure: The plotly figure object of the plotted subnet.
        """
        causal = self.gene_by_motif.get_causal()
        if threshold == "auto":
            threshold = (
                pd.DataFrame(
                    causal.edges(data="weight"), columns=["From", "To", "Weight"]
                )
                .sort_values("Weight")
                .Weight.abs()
                .quantile(0.7)
            )
        subnet = preprocess_net(causal.copy(), threshold)
        subnet = get_subnet(subnet, m, type)
        tf_exp_str = {m: self.get_tf_exp_str(motif, m) for m in motif.cluster_names}
        tf_exp_str["Accessiblity"] = ["TSS Accessibility"]
        tf_exp_mean = {m: self.get_tf_exp_mean(motif, m) for m in motif.cluster_names}
        tf_exp_mean["Accessiblity"] = [50]
        return plotly_networkx_digraph(subnet, tf_exp_str, tf_exp_mean)

    def plotly_gene_exp(self):
        import pandas as pd
        import plotly.express as px

        df = pd.DataFrame(self.gene_annot)
        fig = px.scatter(
            df.groupby("gene_name")[["obs", "pred", "accessibility"]]
            .mean()
            .reset_index(),
            x="obs",
            y="pred",
            color="accessibility",
            hover_name="gene_name",
            labels={
                "obs": "Observed log10 TPM",
                "pred": "Predicted log10 TPM",
                "accessibility": "TSS Accessibility",
            },
            template="plotly_white",
            width=600,
            height=500,
            opacity=0.5,
            marginal_x="histogram",
            marginal_y="histogram",
        )
        # add a text annotation of pearson correlation
        fig.add_annotation(
            x=0.1,
            y=1.0,
            text=f"Cell type: {self.celltype_name}<br />Pearson correlation: {df.groupby('gene_name')[['obs', 'pred']].mean().corr().values[0,1]:.2f}",
            showarrow=False,
            xref="paper",
            yref="paper",
            font=dict(family="Arial", size=16, color="black"),
            align="left",
            bordercolor="black",
            borderwidth=1,
            borderpad=4,
            bgcolor="white",
            opacity=0.8,
        )
        # all font to Arial
        fig.update_layout(
            font_family="Arial",
            font_color="black",
            title_font_family="Arial",
            title_font_color="black",
            legend_title_font_color="black",
            legend_font_color="black",
            xaxis_title_font_family="Arial",
            yaxis_title_font_family="Arial",
            xaxis_title_font_color="black",
            yaxis_title_font_color="black",
            xaxis_tickfont_family="Arial",
            yaxis_tickfont_family="Arial",
            xaxis_tickfont_color="black",
            yaxis_tickfont_color="black",
            xaxis_tickcolor="black",
            yaxis_tickcolor="black",
        )

        return fig

    def plot_gene_regions(self, gene, plotly=False):
        r = self.get_gene_jacobian_summary(gene, "region")
        js = self.get_gene_jacobian(gene)
        r["End"] = self.peak_annot.iloc[r["index"].values].End.values
        r = r[["index", "Chromosome", "Start", "End", "Score"]]
        r = r.query(f"Chromosome==@r.iloc[{self.focus}].Chromosome")
        r_motif = (
            pd.concat([j.data for j in js], axis=0)
            .drop(["Chromosome", "Start", "End", "Gene"], axis=1, errors="ignore")
            .groupby("index")
            .mean()
        )
        r = r.merge(r_motif, left_on="index", right_index=True)
        if plotly:
            return self.plot_region_plotly(r)
        else:
            return self.plot_region(r)

    def plot_region(df):
        # plot the region using rectangles defined by start and end and height defined by score
        # df: dataframe with columns ['Start', 'End', 'Score']
        # return: plot
        df = df.sort_values("Score", ascending=False)
        df["Height"] = df.Score.abs() / df.Score.abs().max()
        df["Width"] = df.End - df.Start
        df["X"] = df.Start
        df["Y"] = df.Height

        fig, ax = plt.subplots(figsize=(10, 2))
        for i, row in df.iterrows():
            ax.add_patch(plt.Rectangle((row.X, 0), row.Width, row.Height, color="red"))
        ax.set_xlim(df.Start.min(), df.End.max())
        # add y=0
        ax.plot([df.Start.min(), df.End.max()], [0, 0], color="black")
        ax.set_ylim(0, 1)
        ax.set_yticks([])
        ax.set_xlabel(f"Genomic Position on Chromosome {df.Chromosome.iloc[0][3:]}")
        # remove top and right spines
        sns.despine(ax=ax, top=True, right=True, left=True, bottom=False)
        return fig, ax

    def plot_region_plotly(self, df: pd.DataFrame) -> go.Figure:
        # Create a subplot with 2 vertical panels; the second panel will be used for gene annotations
        fig = make_subplots(
            rows=2,
            cols=1,
            shared_xaxes=True,
            vertical_spacing=0.1,
            row_heights=[0.8, 0.2],
        )

        hover_text = df.apply(
            lambda row: f"{row['Chromosome']}:{row['Start']}-{row['End']}<br />Top 5 motifs: {', '.join(row.iloc[5:].sort_values()[-5:].index.values)}<br />Bottom 5 motifs: {', '.join(row.iloc[5:].sort_values()[:5].index.values)}",
            axis=1,
        )
        df["HoverText"] = hover_text

        # Process main DataFrame to sort by 'Score' and normalize height
        sorted_df = df.sort_values("Score", ascending=False)
        max_score = sorted_df["Score"].abs().max()
        sorted_df["NormalizedHeight"] = sorted_df["Score"].abs() / max_score

        # Compute genomic span, normalized positions, and heights
        x_positions = (sorted_df["Start"] + sorted_df["End"]) / 2
        heights = sorted_df["NormalizedHeight"]

        # Add bars from x-axis to scores
        x_bar, y_bar = [], []
        for x_coord, y_coord in zip(x_positions, heights):
            x_bar += [x_coord, x_coord, None]
            y_bar += [0, y_coord, None]
        x_bar = x_bar[:-1]
        y_bar = y_bar[:-1]

        # Prepare hover text for main panel
        # Add bar trace for main genomic data
        fig.add_trace(
            go.Scatter(
                x=x_bar, y=y_bar, orientation="v", text=hover_text, hoverinfo="text"
            ),
            row=1,
            col=1,
        )
        # add a scatter trace for each bar for top 10 % positions
        top10 = sorted_df.sort_values("NormalizedHeight", ascending=False).head(
            int(len(sorted_df) * 0.1)
        )

        fig.add_trace(
            go.Scatter(
                x=top10["Start"],
                y=top10["NormalizedHeight"],
                mode="markers",
                marker=dict(color=top10["NormalizedHeight"], colorscale="plotly3"),
                text=top10["HoverText"],
                hoverinfo="text",
            ),
            row=1,
            col=1,
        )
        # Query gene annotations for the same chromosome and genomic range
        gene_start = sorted_df["Start"].min()
        gene_end = sorted_df["End"].max()
        chromosome = sorted_df.iloc[0].Chromosome
        genes_to_show = self.gene_annot[
            (self.gene_annot["Chromosome"] == chromosome)
            & (self.gene_annot["Start"] >= gene_start)
            & (self.gene_annot["Start"] <= gene_end)
        ]

        # Add scatter trace for gene annotations
        fig.add_trace(
            go.Scatter(
                x=genes_to_show["Start"],
                y=genes_to_show["Strand"],
                mode="markers",
                marker=dict(color=genes_to_show["Strand"], colorscale="spectral"),
                text=genes_to_show["gene_name"],
                hoverinfo="text",
            ),
            row=2,
            col=1,
        )

        # set y-axis ticks for gene annotations (0: '+', 1: '-')
        fig.update_yaxes(
            row=2, col=1, tickmode="array", tickvals=[0, 1], ticktext=["+", "-"]
        )

        # Update layout
        chrom_id = sorted_df.iloc[0]["Chromosome"][3:]
        fig.update_layout(
            xaxis=dict(
                range=[sorted_df["Start"].min(), sorted_df["End"].max()],
                title=f"Genomic Position on Chromosome {chrom_id}",
            ),
            yaxis=dict(range=[0, 1.2], showticklabels=False),
            showlegend=False,
            plot_bgcolor="white",
        )

        return fig


class GETCellType(Celltype):
    def __init__(self, celltype, config, s3_file_sys=None):
        """
        Initialize GETCellType instance from a config.

        Parameters
        ----------
        celltype
            Cell type name
        config
            Configuration object
        s3_file_sys
            S3 file system object

        Returns
        -------
        GETCellType
            New instance

        Note
        -----
        This is deprecated and only used for backwards compatibility for demo website.

        """
        features = config.celltype.features

        if features == "NrMotifV1":
            features = np.array(motif_clusters + ["Accessibility"])
        num_region_per_sample = config.celltype.num_region_per_sample
        data_dir = config.celltype.data_dir
        interpret_dir = config.celltype.interpret_dir
        assets_dir = config.celltype.assets_dir
        input = config.celltype.input
        jacob = config.celltype.jacob
        embed = config.celltype.embed
        num_cls = config.celltype.num_cls
        super().__init__(
            features,
            num_region_per_sample,
            celltype,
            data_dir,
            interpret_dir,
            assets_dir,
            input,
            jacob,
            embed,
            num_cls,
            s3_file_sys,
        )


class GETHydraCellType(Celltype):
    """
    Cell type class optimized for hydra model analysis.

    This class extends the base Celltype class with functionality specific to
    the hydra model architecture, including zarr-based data storage and
    efficient processing of large datasets.

    Parameters
    ----------
    celltype
        Name/identifier of the cell type
    zarr_path
        Path to zarr data

    """

    def __init__(self, celltype: str, zarr_path: str, prediction_target: str = "exp"):
        # Initialize parent class attributes that will be needed
        self._gene_by_motif = None
        self.celltype = celltype
        self.celltype_name = celltype
        self.zarr_path = zarr_path
        self.motif = NrMotifV1.load_from_pickle()
        self.assets_dir = "assets/"  # Default value, could be made configurable
        self.s3_file_sys = None  # Default value, could be made configurable
        self.prediction_target = prediction_target

        self._load_zarr_data()

        # Initialize parent class attributes
        self.features = list(self.motif.cluster_names) + ["Accessibility"]
        self.num_features = len(self.features)
        self.num_region_per_sample = self._zarr_data["input"].shape[1]
        if len(self._zarr_data["preds"][self.prediction_target].shape) == 3:
            self.num_cls = self._zarr_data["preds"][self.prediction_target].shape[2]
        else:
            self.num_cls = 1
        self._process_data()

        # Process jacobians if available
        if "jacobians" in self._zarr_data:
            self.jacobs = self._zarr_data["jacobians"]

    def _load_zarr_data(self):
        """Load zarr data from zarr path."""
        logging.info(f"Loading zarr data from {self.zarr_path}")
        self._zarr_data = zarr.open(self.zarr_path, mode="r")
        self.genelist = [
            x.strip(" ") for x in self._zarr_data["available_genes"][:].reshape(-1)
        ]

    def _process_data(self):
        """Process data using vectorized operations."""
        logging.info(f"Processing data for {self.celltype}")

        # Get all data upfront
        focus_idx = 100  # Since this was hardcoded in model interpretation code
        available_genes = self._zarr_data["available_genes"][:].reshape(-1)
        chromosome = self._zarr_data["chromosome"][:].flatten()
        # strip whitespace from chromosome
        chromosome = np.array([x.strip(" ") for x in chromosome])
        peak_coord = self._zarr_data["peak_coord"][:]
        input_data = self._zarr_data["input"][:]
        if self.prediction_target == "exp":
            strands = self._zarr_data["strand"][:].flatten()
            # Get predictions and observations for all genes at once
            pred_values = np.array(
                [
                    self._zarr_data["preds"][self.prediction_target][
                        i, focus_idx, int(strand)
                    ]
                    for i, strand in enumerate(strands)
                ]
            )
            obs_values = np.array(
                [
                    self._zarr_data["obs"][self.prediction_target][
                        i, focus_idx, int(strand)
                    ]
                    for i, strand in enumerate(strands)
                ]
            )
        else:
            pred_values = self._zarr_data["preds"][self.prediction_target][
                :, focus_idx, :
            ].flatten()
            obs_values = self._zarr_data["obs"][self.prediction_target][
                :, focus_idx, :
            ].flatten()
            strands = np.zeros(len(pred_values))
        # Create DataFrame directly
        self.gene_annot = pd.DataFrame(
            {
                "gene_name": [x.strip(" ") for x in available_genes],
                "Chromosome": [x.strip(" ") for x in chromosome],
                "Start": peak_coord[:, focus_idx, 0].astype(int),
                "End": peak_coord[:, focus_idx, 1].astype(int),
                "Strand": strands,
                "pred": pred_values,
                "obs": obs_values,
                "accessibility": input_data[:, focus_idx, -1],
            }
        )

        # Create peak annotation DataFrame
        self.peak_annot = pd.DataFrame(
            {
                "Chromosome": np.repeat(chromosome, self.num_region_per_sample),
                "Start": peak_coord[:, :, 0].flatten().astype(int),
                "End": peak_coord[:, :, 1].flatten().astype(int),
                "Gene": np.repeat(self.genelist, self.num_region_per_sample),
            }
        ).reset_index()

    def load_gene_annot(self):
        """Override parent method - not needed for GETHydraCellType.

        GETHydraCellType creates gene annotations directly from zarr data in
        `_process_data()`, which is called during initialization. The gene annotation
        is already available in `self.gene_annot`.

        If you need to access gene annotations, use `self.gene_annot` directly.
        If you need strand information for a gene, use `self.get_gene_strand(gene_name)`.

        Raises
        ------
        NotImplementedError
            Always raises to indicate this method is not applicable for GETHydraCellType.
        """
        raise NotImplementedError(
            "load_gene_annot() is not needed for GETHydraCellType. "
            "Gene annotations are automatically created from zarr data during initialization. "
            "Access them via `self.gene_annot` or use `self.get_gene_strand(gene_name)` "
            "to get strand information for a specific gene."
        )

    def get_gene_idx(self, gene_name: str) -> np.ndarray:
        """Get all indices for a given gene name.

        Parameters
        ----------
        gene_name
            The name of the gene.

        Returns
        -------
        np.ndarray
            The indices of the gene.
        """
        return self.gene_annot[self.gene_annot["gene_name"] == gene_name].index.values

    def get_gene_strand(self, gene_name: str) -> int:
        """Get strand information for a gene.

        Parameters
        ----------
        gene_name
            The name of the gene.

        Returns
        -------
        int
            The strand. 0 for '+', 1 for '-'.
        """
        v = self.gene_annot[self.gene_annot["gene_name"] == gene_name]["Strand"].values
        if isinstance(v, np.ndarray):
            return v[0]
        else:
            return v

    def get_powerlaw_at_distance(self, distances, min_distance=5000):
        """Get the powerlaw trained on k562 Hi-C data at a given distance.

        Parameters
        ----------
        distances : array-like
            The distances to the TSS in base pairs.
        min_distance : int, default=5000
            Minimum distance threshold in base pairs. Distances below this value will be clipped to this value.

        Returns
        -------
        numpy.ndarray
            Powerlaw contact frequency values corresponding to the input distances.

        Notes
        -----
        The powerlaw parameters are derived from K562 Hi-C data with:
        - gamma = 1.024238616787792
        - scale = 5.9594510043736655

        For distances below min_distance (5kb), values are clipped as the model
        doesn't accurately represent contact frequencies at short distances.
        """
        gamma = 1.024238616787792
        scale = 5.9594510043736655
        # from https://github.com/broadinstitute/ABC-Enhancer-Gene-Prediction/blob/main/config/config.yaml
        # The powerlaw is computed for distances > 5kb. We don't know what the contact freq looks like at < 5kb.
        # So just assume that everything at < 5kb is equal to 5kb.
        # TO DO: get more accurate powerlaw at < 5kb
        distances = np.clip(distances, min_distance, np.Inf)
        log_dists = np.log(distances + 1)

        powerlaw_contact = np.exp(scale + -1 * gamma * log_dists)
        return powerlaw_contact

    def get_gene_jacobian_region_norm(
        self,
        gene_name: str,
        multiply_atac=True,
        layer="region_embed",
    ) -> list[pd.DataFrame]:
        """Get jacobians L2 norm per region for all TSS of a gene.

        This method calculates the L2 norm of jacobians for each region associated with
        a gene's transcription start sites (TSS), and combines this with distance-based
        powerlaw contact frequency and ATAC-seq accessibility data.

        Parameters
        ----------
        gene_name : str
            The name of the gene to analyze.
        multiply_atac : bool, default=True
            Whether to multiply the ATAC signal by the jacobian L2 norm.
        layer : str, default="region_embed"
            The neural network layer from which to extract jacobians.

        Returns
        -------
        list[pd.DataFrame]
            A list of DataFrames (one per TSS) containing:
            - jacobian: L2 norm of the jacobian (optionally multiplied by ATAC signal)
            - distance: distance from region to gene TSS
            - powerlaw: Hi-C based contact frequency estimate
            - ATAC: accessibility signal
            - powerlaw_atac: product of powerlaw and ATAC signal
            - final_prediction: sum of jacobian and powerlaw_atac
        """
        indices = self.get_gene_idx(gene_name)
        strand = self.get_gene_strand(gene_name)
        if self.prediction_target != "exp":
            strand = 0
        jacobians = []
        for k, i in enumerate(indices):
            jacob = self._zarr_data["jacobians"][self.prediction_target][str(strand)][
                layer
            ][i]
            input_data = self._zarr_data["input"][i]
            atac_data = input_data[:, -1]
            jacob_norm = np.linalg.norm(jacob, axis=1)
            if multiply_atac:
                jacob_norm = jacob_norm * atac_data
            # convert to dataframe
            start_idx = k * self.num_region_per_sample
            end_idx = start_idx + self.num_region_per_sample
            peaks_df = self.peak_annot.query(f"Gene == '{gene_name}'").iloc[
                start_idx:end_idx
            ]
            jacob_norm = pd.DataFrame(
                jacob_norm, index=peaks_df.index.values, columns=["jacobian"]
            ).reset_index()
            gene_tss_start = self.gene_annot.query(f"gene_name == '{gene_name}'")[
                "Start"
            ].values[0]
            jacob_norm["distance"] = np.abs(peaks_df.Start.values - gene_tss_start)
            jacob_norm["powerlaw"] = self.get_powerlaw_at_distance(
                jacob_norm["distance"].values
            )
            jacob_norm["ATAC"] = atac_data
            jacob_norm["powerlaw_atac"] = jacob_norm["powerlaw"] * jacob_norm["ATAC"]
            jacob_norm["final_prediction"] = (
                jacob_norm["jacobian"] + jacob_norm["powerlaw_atac"]
            )
            jacobians.append(jacob_norm)
        return jacobians

    def get_gene_jacobian(
        self,
        gene_name: str,
        multiply_input=True,
        multiply_atac=False,
        layer="input/region_motif",
        return_norm=False,
    ) -> list[OneGeneJacobian]:
        """Get jacobians for all TSS of a gene.

        Parameters
        ----------
        gene_name
            The name of the gene.
        multiply_input
            Whether to multiply the input data by the jacobian, by default True

        Returns
        -------
        list
            The jacobians.
        """
        indices = self.get_gene_idx(gene_name)
        strand = self.get_gene_strand(gene_name)
        if self.prediction_target != "exp":
            strand = 0
        jacobians = []
        for k, i in enumerate(indices):
            jacob = self._zarr_data["jacobians"][self.prediction_target][str(strand)][
                layer
            ][i]
            input_data = self._zarr_data["input"][i]

            if multiply_input:
                jacob = jacob * input_data
            if multiply_atac:
                jacob = jacob * self._zarr_data["input"][i][:, -1]
            if return_norm:
                jacob = np.linalg.norm(jacob, axis=1)
                return jacob
            start_idx = k * self.num_region_per_sample
            end_idx = start_idx + self.num_region_per_sample
            jacobians.append(
                OneGeneJacobian(
                    gene_name=gene_name,
                    data=jacob,
                    region=self.peak_annot.query(f"Gene == '{gene_name}'").iloc[
                        start_idx:end_idx
                    ],
                    features=self.features,
                    num_cls=self.num_cls,
                    num_region_per_sample=self.num_region_per_sample,
                    num_features=self.num_features,
                )
            )

        return jacobians

    def get_gene_chromosome(self, gene_name: str) -> str:
        """Get the chromosome of a gene."""
        return self.gene_annot[self.gene_annot["gene_name"] == gene_name][
            "Chromosome"
        ].values[0]

    def __repr__(self) -> str:
        return f"""GETHydraCelltype: {self.celltype}
        Zarr path: {self.zarr_path}
        Number of regions per sample: {self.num_region_per_sample}
        Number of features: {self.num_features}
        Number of genes: {len(self.genelist)}
        Number of peaks: {len(self.peak_annot)}
        """

    @classmethod
    def from_config(cls, cfg, celltype=None, zarr_path=None, prediction_target=None):
        """
        Create GETHydraCellType instance from configuration.

        Parameters
        ----------
        cfg
            Configuration object
        celltype
            Cell type name
        zarr_path
            Path to zarr data
        prediction_target
            Prediction target. Default is the first loss component in the model.

        Returns
        -------
        GETHydraCellType
            New instance
        """
        if celltype is None:
            celltype = f"{cfg.dataset.leave_out_celltypes}"
        if zarr_path is None:
            zarr_path = f"{cfg.machine.output_dir}/{cfg.run.project_name}/{cfg.run.run_name}/{celltype}.zarr"
        logger.info(
            f"Creating GETHydraCellType instance for {celltype}, loading data from {zarr_path}"
        )
        if prediction_target is None:
            prediction_target = list(cfg.model.cfg.loss.components.keys())[0]
        return cls(
            celltype=celltype,
            zarr_path=zarr_path,
            prediction_target=prediction_target,
        )

    def get_gene_by_motif(self, overwrite: bool = False) -> GeneByMotif:
        """
        Calculate gene by motif data for hydra-based cell-type.

        This method calculates the gene by motif data for a GETHydraCellType.
        It averages the jacobians across all TSS sites for each gene.

        Parameters
        ----------
        overwrite
            Whether to overwrite existing data, by default False

        Returns
        -------
        GeneByMotif
            The gene by motif data.
        """
        if "gene_by_motif" in self._zarr_data and not overwrite:
            gene_by_motif_values = self._zarr_data["gene_by_motif"][:]
            if gene_by_motif_values.shape[0] == len(self.gene_annot):
                gene_by_motif_df = pd.DataFrame(
                    gene_by_motif_values,
                    index=self.gene_annot["gene_name"].values,
                    columns=self.features,
                )
            else:
                gene_by_motif_df = pd.DataFrame(
                    gene_by_motif_values,
                    index=self.gene_annot["gene_name"].unique(),
                    columns=self.features,
                )
        else:
            all_genes = self.gene_annot["gene_name"].unique()
            motif_data = []

            for gene in tqdm(all_genes, desc="Processing genes for motif analysis"):
                try:
                    # Get jacobians for all TSS
                    jacobs = self.get_gene_jacobian(gene, multiply_input=True)

                    # Calculate motif summary for each TSS and average them
                    motif_summaries = []
                    for jacob in jacobs:
                        motif_summary = jacob.motif_summary(stats="mean")
                        motif_summaries.append(motif_summary)

                    # Average across all TSS
                    avg_motif_summary = pd.concat(motif_summaries, axis=1).mean(axis=1)
                    motif_data.append(avg_motif_summary)
                except Exception as e:
                    print(f"Error processing gene {gene}: {str(e)}")
                    continue

            # Create DataFrame with gene by motif data
            gene_by_motif_df = pd.DataFrame(motif_data, index=all_genes)

            # save to zarr
            self._zarr_data = zarr.open(self.zarr_path, mode="a")
            self._zarr_data.create_array(
                "gene_by_motif", data=gene_by_motif_df.values, overwrite=True
            )
        self._gene_by_motif = GeneByMotif(
            celltype=self.celltype,
            interpret_dir=Path(self.zarr_path).parent,
            jacob=gene_by_motif_df,
            s3_file_sys=self.s3_file_sys,
            zarr_data_path=self.zarr_path,
        )
        return self._gene_by_motif

    def plot_gene_motifs(self, gene, motif=None, overwrite: bool = False, n: int = 10):
        """
        Plot top N motifs for a gene. Based on the absolute mean of the motif scores.

        Parameters
        ----------
        gene
            Gene name
        motif
            Motif instance. Not used, self.motif is used instead.
        overwrite
            Whether to overwrite existing plots, by default False
        n
            Number of motifs to plot, by default 10
        """
        r = self.get_gene_jacobian_summary(gene, "motif")
        m = r.sort_values(ascending=False).head(n).index.values
        # remove 'Accessibility' from the list
        m = [m_i for m_i in m if m_i != "Accessibility"]
        fig, ax = plt.subplots(2, 5, figsize=(10, 4), sharex=False, sharey=False)

        for i, m_i in enumerate(m):
            if (
                not path_exists_with_s3(
                    file_path=f'{self.assets_dir}{m_i.replace("/", "_")}.png',
                    s3_file_sys=self.s3_file_sys,
                )
                or overwrite
            ):
                try:
                    self.motif.get_motif_cluster_by_name(m_i).seed_motif.plot_logo(
                        filename=f'{self.assets_dir}{m_i.replace("/", "_")}.png',
                        logo_title="",
                        size="medium",
                        ic_scale=True,
                    )
                except Exception as e:
                    print(f"Error plotting motif {m_i}: {str(e)}")

            if self.s3_file_sys:
                with self.s3_file_sys.open(
                    f'{self.assets_dir}{m_i.replace("/", "_")}.png', "rb"
                ) as f:
                    img = plt.imread(BytesIO(f.read()))
            else:
                img = plt.imread(f'{self.assets_dir}{m_i.replace("/", "_")}.png')

            ax[i // 5][i % 5].imshow(img)
            ax[i // 5][i % 5].axis("off")

            if m_i in self.motif.cluster_gene_list:
                motif_cluster_genes = self.motif.cluster_gene_list[m_i]
                try:
                    ax[i // 5][i % 5].set_title(
                        f"{m_i}:{self.get_highest_exp_genes(motif_cluster_genes)}"
                    )
                except Exception as e:
                    ax[i // 5][i % 5].set_title(f"{m_i}")
            else:
                ax[i // 5][i % 5].set_title(f"{m_i}")

        return fig, ax


def celltype_factory(celltype_class, cell_id, config, zarr_path=None, motif_path=None):
    if celltype_class == "GETCelltype":
        return GETCellType(cell_id, config)
    elif celltype_class == "GETHydraCelltype":
        return GETHydraCellType(cell_id, config, zarr_path, motif_path)
    else:
        raise ValueError(f"Unknown celltype class: {celltype_class}")


class GETDemoLoader:
    def __init__(self):
        import s3fs

        from gcell.config.config import load_config

        print("""
This class will load pre-inferenced data from s3. The checkpoint is a binary atac get model trained on human fetal and adult data.
Gencode v40 was used for gene annotation. The annotation file is stored in a feather file at zenodo: https://zenodo.org/records/14635090/files/gencode.v40.hg38.feather?download=1
The data is stored in the s3 bucket: s3://2023-get-xf2217/get_demo/
Note that not all celltypes have observed expression. In those cases, the observed expression is set to 0.
You can use `available_celltypes` to see which celltypes are available, and `load_celltype` to load a celltype.
              """)
        self.cfg = load_config("s3_interpret")
        self.s3_file_sys = s3fs.S3FileSystem(anon=True)
        self.cfg.celltype.data_dir = (
            f"{self.cfg.s3_uri}/pretrain_human_bingren_shendure_apr2023/fetal_adult/"
        )
        self.cfg.celltype.interpret_dir = (
            f"{self.cfg.s3_uri}/Interpretation_all_hg38_allembed_v4_natac/"
        )
        self.cfg.celltype.motif_dir = (
            f"{self.cfg.s3_uri}/interpret_natac/motif-clustering/"
        )
        self.cfg.celltype.assets_dir = f"{self.cfg.s3_uri}/assets/"
        cell_type_annot = pd.read_csv(
            self.cfg.celltype.data_dir.split("fetal_adult")[0]
            + "data/cell_type_pretrain_human_bingren_shendure_apr2023.txt"
        )
        self.cell_type_id_to_name = dict(
            zip(cell_type_annot["id"], cell_type_annot["celltype"])
        )
        self.cell_type_name_to_id = dict(
            zip(cell_type_annot["celltype"], cell_type_annot["id"])
        )
        self.available_celltypes = sorted(
            [
                self.cell_type_id_to_name[f.split("/")[-1]]
                for f in self.s3_file_sys.glob(self.cfg.celltype.interpret_dir + "*")
            ]
        )

    def load_celltype(self, celltype_name, jacob=False, embed=False):
        """
        Load a celltype from the demo dataset.

        Parameters
        ----------
        celltype_name
            The name of the celltype to load.
        jacob
            Whether to load the Jacobian data.
        embed
            Whether to load the embedding data.

        Returns
        -------
        GETCellType
            The celltype instance.
        """
        celltype_id = self.cell_type_name_to_id[celltype_name]
        if jacob:
            self.cfg.celltype.jacob = True
        if embed:
            self.cfg.celltype.embed = True
        cell = GETCellType(celltype_id, self.cfg, s3_file_sys=self.s3_file_sys)
        return cell
