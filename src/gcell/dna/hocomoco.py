"""
HOCOMOCO motif database interface for Caesar.

This module provides a unified interface for loading and working with HOCOMOCO
transcription factor motif databases, supporting both original PWM format and
pre-aligned tensor format.
"""

import json
import logging
from pathlib import Path
from typing import TYPE_CHECKING, Union

import numpy as np
import torch

from gcell._settings import get_setting

from .motif_utils import (
    compute_pvalue_mapping,
    create_reverse_complement_motifs,
    filter_motifs,
    pad_motifs_to_same_length,
    pvalue_from_logpdf,
    read_annotations,
    read_pwms,
)

if TYPE_CHECKING:
    import pandas as pd

logger = logging.getLogger(__name__)


class HocomocoIO:
    """
    Unified interface for HOCOMOCO motif database access.

    Supports both original PWM format and pre-aligned tensor format,
    providing a consistent interface for motif analysis pipelines.
    """

    def __init__(
        self,
        aligned_motif_path: str | None = None,
        pwm_path: str | None = None,
        annotation_path: str | None = None,
        use_aligned: bool = True,
    ):
        """
        Initialize HOCOMOCO motif database interface.

        Args:
            aligned_motif_path: Path to pre-aligned motif tensor (.pt file)
            pwm_path: Path to original PWM file
            annotation_path: Path to annotation file (required for PWM source)
            use_aligned: Whether to prefer aligned motifs when available
        """
        self.aligned_motif_path = aligned_motif_path
        self.pwm_path = pwm_path
        self.annotation_path = annotation_path
        self.use_aligned = use_aligned

        # Resolve defaults using gcell settings to avoid brittle relative paths
        base_dir = Path(get_setting("annotation_dir")) / "hocomoco"
        base_dir.mkdir(parents=True, exist_ok=True)

        if not aligned_motif_path:
            self.aligned_motif_path = str(base_dir / "motifs_with_rc_aligned.pt")
        if not pwm_path:
            self.pwm_path = str(base_dir / "H13CORE_pwms.txt")
        if not annotation_path:
            self.annotation_path = str(base_dir / "H13CORE_annotation.jsonl")

        # P-value mapping cache directory
        self._pvalue_cache_dir = base_dir / "pvalue_mappings"
        self._pvalue_cache_dir.mkdir(parents=True, exist_ok=True)

        # Cached data
        self._motif_names: list[str] | None = None
        self._motif_kernels: torch.Tensor | None = None
        self._annotations: dict | None = None
        self._similarity_matrix: np.ndarray | None = None
        self._significance_thresholds: dict | None = None
        self._loaded_with_aligned: bool | None = None
        self._tf_name_mapping: dict[str, str] | None = None
        self._pvalue_mappings: dict[str, tuple[int, np.ndarray, float]] | None = None

    @property
    def motif_names(self) -> list[str]:
        """Get list of motif names."""
        if self._motif_names is None:
            self._load_motifs()
        return self._motif_names

    @property
    def motif_kernels(self) -> torch.Tensor:
        """Get motif kernels tensor (N, L, 4)."""
        if self._motif_kernels is None:
            self._load_motifs()
        return self._motif_kernels

    @property
    def similarity_matrix(self) -> np.ndarray | None:
        """Get motif similarity matrix (only available for aligned motifs)."""
        if self._similarity_matrix is None and self._loaded_with_aligned:
            self._load_motifs()
        return self._similarity_matrix

    @property
    def significance_thresholds(self) -> dict | None:
        """Get significance thresholds (available for both aligned and PWM sources)."""
        if self._significance_thresholds is None:
            if self._loaded_with_aligned is None:
                self._load_motifs()
            elif self._loaded_with_aligned:
                # Try to load thresholds from PWM source for aligned motifs
                self._load_significance_thresholds_from_pwm()
            else:
                self._load_motifs()
        return self._significance_thresholds

    @property
    def annotations(self) -> dict | None:
        """Get motif annotations (available for both aligned and PWM sources)."""
        if self._annotations is None:
            if self._loaded_with_aligned is None:
                self._load_motifs()
            elif self._loaded_with_aligned:
                # Try to load annotations from PWM source for aligned motifs
                self._load_significance_thresholds_from_pwm()
            else:
                self._load_motifs()
        return self._annotations

    @property
    def tf_names_list(self) -> list[str]:
        """Get list of unique TF gene names from all motifs."""
        # Load TF name mapping if not already loaded
        if self._tf_name_mapping is None:
            self._load_tf_name_mapping()

        # Get unique TF names from the mapping values
        unique_tf_names = sorted(set(self._tf_name_mapping.values()))
        return unique_tf_names

    def _load_motifs(self):
        """Load motifs from either aligned or PWM source."""
        # Try aligned motifs first if preferred and available
        if self.use_aligned and self._try_load_aligned():
            logger.info("Loaded motifs from pre-aligned tensor")
            self._loaded_with_aligned = True
            return

        # Fall back to PWM source
        if self._try_load_pwm_source():
            logger.info("Loaded motifs from PWM source")
            self._loaded_with_aligned = False
            return

        raise RuntimeError(
            "Could not load motifs from either aligned tensor or PWM source. "
            f"Aligned path: {self.aligned_motif_path}, "
            f"PWM path: {self.pwm_path}, "
            f"Annotation path: {self.annotation_path}. "
            "Place HOCOMOCO files under the configured annotation_dir/hocomoco or pass explicit paths, "
            "and optionally set GCELL_ANNOTATION_DIR to override the base directory."
        )

    def _try_load_aligned(self) -> bool:
        """Try to load pre-aligned motifs."""
        try:
            if not Path(self.aligned_motif_path).exists():
                return False

            logger.info(f"Loading pre-aligned motifs from: {self.aligned_motif_path}")
            motif_data = torch.load(
                self.aligned_motif_path, map_location="cpu", weights_only=False
            )

            self._motif_kernels = motif_data["motif_kernels"]
            self._motif_names = motif_data["motif_names"]
            self._similarity_matrix = motif_data.get("similarity_matrix", None)

            # Try to load significance thresholds from PWM source even when using aligned motifs
            self._load_significance_thresholds_from_pwm()

            logger.info(
                f"Loaded {len(self._motif_names)} aligned motifs with shape {self._motif_kernels.shape}"
            )
            return True

        except Exception as e:
            logger.warning(f"Failed to load aligned motifs: {e}")
            return False

    def _load_significance_thresholds_from_pwm(self):
        """Load significance thresholds from PWM source for aligned motifs."""
        try:
            if (
                not Path(self.pwm_path).exists()
                or not Path(self.annotation_path).exists()
            ):
                logger.warning(
                    "PWM or annotation files not available for threshold loading"
                )
                return

            logger.info("Loading significance thresholds from PWM source")

            # Load raw data
            motif_collection = read_pwms(self.pwm_path)
            annotations = read_annotations(self.annotation_path)

            # Process motifs to get thresholds (no filtering, no RC)
            processed_data = pad_motifs_to_same_length(motif_collection, annotations)
            (
                pwm_motif_names,
                _,
                _,
                significance_thresholds,
                _,
                _,
            ) = processed_data

            # Map thresholds to aligned motif order
            # For RC motifs (ending with _RC), use the threshold from the forward motif
            aligned_thresholds = {}
            for p_value, threshold_array in significance_thresholds.items():
                aligned_threshold_list = []
                for motif_name in self._motif_names:
                    if motif_name.endswith("_RC"):
                        # Remove _RC suffix to get original motif name
                        original_motif = motif_name[:-3]
                    else:
                        original_motif = motif_name

                    # Find the index in PWM motifs
                    try:
                        pwm_index = pwm_motif_names.index(original_motif)
                        aligned_threshold_list.append(threshold_array[pwm_index])
                    except ValueError:
                        logger.warning(
                            f"Motif {original_motif} not found in PWM thresholds"
                        )
                        # Use a default threshold (could be 0 or the mean)
                        aligned_threshold_list.append(0.0)

                aligned_thresholds[p_value] = np.array(aligned_threshold_list)

            self._significance_thresholds = aligned_thresholds
            self._annotations = annotations

            logger.info(
                f"Loaded significance thresholds for {len(pwm_motif_names)} motifs"
            )

        except Exception as e:
            logger.warning(
                f"Failed to load significance thresholds from PWM source: {e}"
            )

    def _try_load_pwm_source(self) -> bool:
        """Try to load motifs from PWM source."""
        try:
            if (
                not Path(self.pwm_path).exists()
                or not Path(self.annotation_path).exists()
            ):
                return False

            logger.info(f"Loading motifs from PWM source: {self.pwm_path}")

            # Load raw data
            motif_collection = read_pwms(self.pwm_path)
            self._annotations = read_annotations(self.annotation_path)

            # Process motifs with default settings (no filtering, no RC)
            processed_data = pad_motifs_to_same_length(
                motif_collection, self._annotations
            )
            (
                self._motif_names,
                padded_kernels,
                _,
                self._significance_thresholds,
                _,
                _,
            ) = processed_data

            # Convert to tensor
            self._motif_kernels = torch.from_numpy(padded_kernels).float()

            logger.info(
                f"Loaded {len(self._motif_names)} motifs from PWM source with shape {self._motif_kernels.shape}"
            )
            return True

        except Exception as e:
            logger.warning(f"Failed to load PWM source: {e}")
            return False

    def filter_motifs(self, motif_selection: str) -> "HocomocoIO":
        """
        Create a new HocomocoIO instance with filtered motifs.

        Args:
            motif_selection: Motif selection criteria

        Returns:
            New HocomocoIO instance with filtered motifs
        """
        if self._loaded_with_aligned or self._loaded_with_aligned is None:
            raise ValueError(
                "Motif filtering is only supported when loading from PWM source"
            )

        logger.info(f"Filtering motifs with criteria: {motif_selection}")

        # Reload from PWM source with filtering
        motif_collection = read_pwms(self.pwm_path)
        annotations = read_annotations(self.annotation_path)

        # Apply filtering
        filtered_collection, filtered_annotations = filter_motifs(
            motif_collection, annotations, motif_selection
        )

        # Process filtered motifs
        processed_data = pad_motifs_to_same_length(
            filtered_collection, filtered_annotations
        )
        motif_names, padded_kernels, _, significance_thresholds, _, _ = processed_data

        # Create new instance
        filtered_io = HocomocoIO(use_aligned=False)
        filtered_io._motif_names = motif_names
        filtered_io._motif_kernels = torch.from_numpy(padded_kernels).float()
        filtered_io._annotations = filtered_annotations
        filtered_io._significance_thresholds = significance_thresholds
        filtered_io._loaded_with_aligned = False

        logger.info(f"Filtered to {len(motif_names)} motifs")
        return filtered_io

    def create_reverse_complements(self) -> "HocomocoIO":
        """
        Create a new HocomocoIO instance with reverse complement motifs added.

        Returns:
            New HocomocoIO instance with reverse complement motifs
        """
        if self._loaded_with_aligned or self._loaded_with_aligned is None:
            logger.info("Reverse complement creation is already done in aligned motifs")
            return self

        logger.info("Creating reverse complement motifs")

        # Reload from PWM source
        motif_collection = read_pwms(self.pwm_path)
        annotations = read_annotations(self.annotation_path)

        # Create reverse complements
        rc_collection, rc_annotations = create_reverse_complement_motifs(
            motif_collection, annotations
        )

        # Process motifs with RC
        processed_data = pad_motifs_to_same_length(rc_collection, rc_annotations)
        motif_names, padded_kernels, _, significance_thresholds, _, _ = processed_data

        # Create new instance
        rc_io = HocomocoIO(use_aligned=False)
        rc_io._motif_names = motif_names
        rc_io._motif_kernels = torch.from_numpy(padded_kernels).float()
        rc_io._annotations = rc_annotations
        rc_io._significance_thresholds = significance_thresholds
        rc_io._loaded_with_aligned = False

        logger.info(f"Created {len(motif_names)} motifs with reverse complements")
        return rc_io

    def get_threshold(self, p_value: str = "p0.0001") -> torch.Tensor | None:
        """
        Get significance thresholds for a given p-value.

        Args:
            p_value: P-value threshold string (e.g., "p0.0001")

        Returns:
            Tensor of thresholds or None if not available
        """
        if self.significance_thresholds and p_value in self.significance_thresholds:
            return torch.from_numpy(self.significance_thresholds[p_value]).float()
        return None

    def prepare_for_gpu(
        self,
        device: str = "cuda",
        dtype: torch.dtype = torch.float32,
        p_value: str = "p0.0001",
    ) -> tuple[torch.Tensor, torch.Tensor | None]:
        """
        Prepare motif data for GPU computation.

        Args:
            device: Target device (cuda/cpu)
            dtype: Numerical precision
            p_value: P-value threshold to prepare

        Returns:
            Tuple of (motif_kernels_gpu, thresholds_gpu)
        """
        kernels_gpu = self.motif_kernels.to(device=device, dtype=dtype)

        thresholds_gpu = None
        threshold_tensor = self.get_threshold(p_value)
        if threshold_tensor is not None:
            thresholds_gpu = threshold_tensor.to(device=device, dtype=dtype)

        return kernels_gpu, thresholds_gpu

    def get_motif_index(self, motif_name: str) -> int | None:
        """
        Get the index of a motif by name.

        Args:
            motif_name: Name of the motif to find

        Returns:
            Index of the motif or None if not found
        """
        try:
            return self.motif_names.index(motif_name)
        except ValueError:
            return None

    def get_motif_by_index(self, index: int) -> tuple[str, torch.Tensor]:
        """
        Get motif name and kernel by index.

        Args:
            index: Index of the motif

        Returns:
            Tuple of (motif_name, motif_kernel)
        """
        return self.motif_names[index], self.motif_kernels[index]

    def _load_tf_name_mapping(self):
        """Load mapping from motif names to TF gene names from annotation file."""
        if not self.annotation_path:
            logger.warning(
                "Annotation path not available. TF name mapping will be empty."
            )
            self._tf_name_mapping = {}
            return

        annotation_path = Path(self.annotation_path)
        if not annotation_path.exists():
            logger.warning(
                f"Annotation file not found: {self.annotation_path}. TF name mapping will be empty."
            )
            self._tf_name_mapping = {}
            return

        self._tf_name_mapping = {}
        try:
            with annotation_path.open() as f:
                for line in f:
                    annotation = json.loads(line.strip())
                    motif_name = annotation.get("name")
                    if not motif_name:
                        continue

                    # Try to get the gene symbol from masterlist_info first (preferred)
                    masterlist_info = annotation.get("masterlist_info", {})
                    species_info = masterlist_info.get("species", {})
                    human_info = species_info.get("HUMAN", {})
                    gene_symbol = human_info.get("gene_symbol")

                    if gene_symbol:
                        self._tf_name_mapping[motif_name] = gene_symbol
                    else:
                        # Fall back to tf field if gene_symbol not available
                        tf_name = annotation.get("tf")
                        if tf_name:
                            self._tf_name_mapping[motif_name] = tf_name

            logger.info(
                f"Loaded TF name mapping for {len(self._tf_name_mapping)} motifs"
            )

        except Exception as e:
            logger.warning(f"Failed to load TF name mapping: {e}")
            self._tf_name_mapping = {}

    def get_tf_name_from_motif(self, motif_name: str) -> str:
        """
        Extract transcription factor gene name from motif using the HOCOMOCO annotation file.

        Args:
            motif_name: Name of the motif (e.g., "ZN239.H13CORE.0.P.C" or "ZN239.H13CORE.0.P.C_RC")

        Returns:
            Transcription factor gene name (e.g., "ZNF239")

        Raises:
            ValueError: If motif not found in annotation file or annotation path not available
        """
        # Load TF name mapping if not already loaded
        if self._tf_name_mapping is None:
            self._load_tf_name_mapping()

        # Handle reverse complement motifs by removing _RC suffix
        lookup_name = motif_name
        if motif_name.endswith("_RC"):
            lookup_name = motif_name[:-3]  # Remove "_RC" suffix

        # Look up motif in cached mapping
        if lookup_name in self._tf_name_mapping:
            return self._tf_name_mapping[lookup_name]

        raise ValueError(
            f"Motif {motif_name} (lookup: {lookup_name}) not found in annotation file"
        )

    def resolve_motif_index(self, motif_identifier: str | int) -> int:
        """
        Resolves a motif name or index to a numeric index.

        Args:
            motif_identifier: Motif name or index

        Returns:
            Resolved motif index

        Raises:
            ValueError: If motif cannot be resolved
        """
        motif_names = self.motif_names

        if isinstance(motif_identifier, int):
            if 0 <= motif_identifier < len(motif_names):
                return motif_identifier
            else:
                raise ValueError(
                    f"Motif index {motif_identifier} out of range [0, {len(motif_names) - 1}]"
                )

        # Try exact name match first
        try:
            return motif_names.index(motif_identifier)
        except ValueError:
            pass

        # Try prefix matching with preference for forward strand (non-RC)
        prefix_matches = []
        for i, name in enumerate(motif_names):
            if name.startswith(motif_identifier + ".") or name.startswith(
                motif_identifier + "_"
            ):
                prefix_matches.append((i, name))

        if len(prefix_matches) == 1:
            return prefix_matches[0][0]
        elif len(prefix_matches) > 1:
            # Prefer non-reverse complement matches
            non_rc_matches = [
                (i, name) for i, name in prefix_matches if not name.endswith("_RC")
            ]
            if len(non_rc_matches) == 1:
                return non_rc_matches[0][0]
            elif len(non_rc_matches) > 1:
                # If multiple non-RC matches, take the first one
                return non_rc_matches[0][0]
            else:
                # Only RC matches available, take the first one
                return prefix_matches[0][0]

        # Try partial name match (case-insensitive)
        partial_matches = []
        motif_identifier_upper = motif_identifier.upper()
        for i, name in enumerate(motif_names):
            if motif_identifier_upper in name.upper():
                partial_matches.append((i, name))

        if len(partial_matches) == 1:
            return partial_matches[0][0]
        elif len(partial_matches) > 1:
            # Prefer non-reverse complement matches
            non_rc_matches = [
                (i, name) for i, name in partial_matches if not name.endswith("_RC")
            ]
            if len(non_rc_matches) == 1:
                return non_rc_matches[0][0]
            elif len(non_rc_matches) > 1:
                # Multiple matches - provide helpful error message
                match_names = [name for _, name in non_rc_matches]
                raise ValueError(
                    f"Ambiguous motif name '{motif_identifier}'. Non-RC matches: {match_names[:5]}"
                )
            else:
                # Only RC matches
                match_names = [name for _, name in partial_matches]
                raise ValueError(
                    f"Ambiguous motif name '{motif_identifier}'. RC matches: {match_names[:5]}"
                )

        raise ValueError(
            f"Motif '{motif_identifier}' not found in {len(motif_names)} available motifs"
        )

    def get_tf_tss_scores(
        self,
        zarr_io,
        gencode=None,
        window_size: int = 1000,
        conv_size: int = 100,
        aggregator: str = "by_key",
        output_format: str = "raw_dict",
        normalize: bool = True,
        verbose: bool = True,
        return_dataframe: bool = False,
    ) -> Union[dict[str, np.ndarray], "pd.DataFrame"]:
        """
        Extract TSS expression scores for all TFs from a zarr track.

        Args:
            zarr_io: CelltypeDenseZarrIO instance with expression data
            gencode: Gencode instance for gene annotation (auto-initialized from zarr assembly if None)
            window_size: Size of window around TSS (default: 1000, so ±500 bp)
            conv_size: Convolution size for signal processing (default: 100)
            aggregator: Aggregation method (default: 'by_key')
            output_format: Output format (default: 'raw_dict')
            normalize: Whether to normalize by library size (default: True)
            verbose: Whether to show progress and errors (default: True)
            return_dataframe: Whether to return pandas DataFrame instead of dict (default: False)

        Returns:
            Dictionary mapping TF names to max TSS scores across all TSSs, or pandas DataFrame if return_dataframe=True
        """
        try:
            from tqdm import tqdm
        except ImportError:
            # Fallback if tqdm is not available
            def tqdm(iterable, **kwargs):
                return iterable

        # Auto-initialize gencode if not provided
        if gencode is None:
            try:
                from gcell.rna.gencode import Gencode

                # Try to get assembly from zarr_io attributes
                assembly = getattr(zarr_io, "assembly", None)
                if assembly is None:
                    # Try common attribute names
                    for attr in ["genome_assembly", "genome", "reference_genome"]:
                        assembly = getattr(zarr_io, attr, None)
                        if assembly:
                            break

                if assembly is None:
                    raise ValueError(
                        "Could not determine assembly from zarr_io. "
                        "Please provide gencode instance or ensure zarr_io has assembly attribute."
                    )

                gencode = Gencode(assembly=assembly)
                if verbose:
                    logger.info(f"Auto-initialized Gencode with assembly: {assembly}")

            except ImportError:
                raise ImportError(
                    "gcell.rna.gencode.Gencode is required for auto-initialization. "
                    "Please install gcell or provide gencode instance."
                )

        def flatten_dict_mean(d):
            """Take a dict of arrays, return the mean of each array concatenated in key order."""
            return np.array([float(np.mean(v)) for k, v in d.items()])

        tf_tss_dict = {}
        half_window = window_size // 2

        iterator = (
            tqdm(self.tf_names_list, desc="Processing TFs")
            if verbose
            else self.tf_names_list
        )

        for tf in iterator:
            try:
                gene = gencode.get_gene(tf)
                gene_tss_exp = []

                for tss in gene.tss_list.Start.values:
                    # Extract track data around TSS
                    track_kwargs = {
                        "aggregator": aggregator,
                        "conv_size": conv_size,
                        "output_format": output_format,
                    }

                    if normalize and hasattr(zarr_io, "libsize"):
                        track_kwargs["normalize_factor"] = zarr_io.libsize

                    d = zarr_io.get_track(
                        gene.chrom, tss - half_window, tss + half_window, **track_kwargs
                    )
                    gene_tss_exp.append(flatten_dict_mean(d))

                # Take maximum across all TSSs for this gene
                if gene_tss_exp:
                    tf_tss_dict[tf] = np.stack(gene_tss_exp).max(axis=0)

            except Exception as e:
                if verbose:
                    logger.warning(f"Error processing {tf}: {e}")
                continue

        if verbose:
            logger.info(
                f"Successfully processed {len(tf_tss_dict)} out of {len(self.tf_names_list)} TFs"
            )

        # Return DataFrame if requested
        if return_dataframe:
            try:
                import pandas as pd

                # Get column names from zarr_io (cell type IDs)
                columns = getattr(zarr_io, "ids", None)
                if columns is None:
                    # Try other common attribute names
                    for attr in ["cell_types", "celltype_ids", "labels"]:
                        columns = getattr(zarr_io, attr, None)
                        if columns is not None:
                            break

                if columns is None:
                    # Use numeric indices if no labels available
                    if tf_tss_dict:
                        first_array = next(iter(tf_tss_dict.values()))
                        columns = [f"celltype_{i}" for i in range(len(first_array))]
                    else:
                        columns = []

                # Create DataFrame with TF names as index and cell types as columns
                df = pd.DataFrame(tf_tss_dict, index=columns).T
                return df

            except ImportError:
                raise ImportError(
                    "pandas is required for return_dataframe=True. "
                    "Please install pandas or set return_dataframe=False."
                )

        return tf_tss_dict

    def parse_motif_pairs(self, pairs_string: str) -> list[tuple[int, int]]:
        """
        Parses a string specification of motif pairs into a list of index tuples.

        Supports formats:
        - "motif1,motif2;motif3,motif4" (semicolon-separated pairs)
        - "0,1;2,3" (index-based pairs)
        - "CTCF,TP53;NFKB1,STAT1" (name-based pairs)

        Args:
            pairs_string: String specification of motif pairs

        Returns:
            List of (motif_i_index, motif_j_index) tuples
        """
        if not pairs_string:
            return []

        pairs = []
        pair_specifications = pairs_string.split(";")

        for pair_spec in pair_specifications:
            try:
                motif_i, motif_j = pair_spec.strip().split(",")
                motif_i, motif_j = motif_i.strip(), motif_j.strip()

                # Try to convert to integers if possible, otherwise resolve names
                try:
                    motif_i_idx = int(motif_i)
                except ValueError:
                    motif_i_idx = self.resolve_motif_index(motif_i)

                try:
                    motif_j_idx = int(motif_j)
                except ValueError:
                    motif_j_idx = self.resolve_motif_index(motif_j)

                pairs.append((motif_i_idx, motif_j_idx))

            except ValueError as e:
                logger.warning(
                    f"Invalid motif pair specification: '{pair_spec}'. Error: {e}"
                )
                continue

        return pairs

    def __len__(self) -> int:
        """Get number of motifs."""
        return len(self.motif_names)

    def _get_pvalue_mapping_cache_path(self, motif_name: str) -> Path:
        """Get cache path for a motif's p-value mapping."""
        # Sanitize motif name for filesystem
        safe_name = motif_name.replace("/", "_").replace("\\", "_").replace(":", "_")
        return self._pvalue_cache_dir / f"{safe_name}.npz"

    def _load_pvalue_mappings(self):
        """Load or compute p-value mappings for all motifs."""
        if self._pvalue_mappings is not None:
            return

        logger.info("Loading p-value mappings for motifs")
        motif_names = self.motif_names
        self._pvalue_mappings = {}

        # Load raw PWMs
        motif_collection = read_pwms(self.pwm_path)

        for motif_name in motif_names:
            # Handle RC motifs - use the original motif for mapping
            original_name = motif_name
            if motif_name.endswith("_RC"):
                original_name = motif_name[:-3]

            cache_path = self._get_pvalue_mapping_cache_path(original_name)

            # Try to load from cache
            if cache_path.exists():
                try:
                    cached_data = np.load(cache_path)
                    smallest = int(cached_data["smallest"])
                    logpdf = cached_data["logpdf"]
                    bin_size = float(cached_data["bin_size"])
                    self._pvalue_mappings[motif_name] = (smallest, logpdf, bin_size)
                    continue
                except Exception as e:
                    logger.warning(
                        f"Failed to load cached p-value mapping for {motif_name}: {e}"
                    )

            # Compute mapping if not cached or cache failed
            if original_name in motif_collection:
                pwm = motif_collection[original_name]
                try:
                    smallest, logpdf, bin_size = compute_pvalue_mapping(pwm)
                    self._pvalue_mappings[motif_name] = (smallest, logpdf, bin_size)

                    # Cache the result
                    np.savez(
                        cache_path,
                        smallest=smallest,
                        logpdf=logpdf,
                        bin_size=bin_size,
                    )
                except Exception as e:
                    logger.warning(
                        f"Failed to compute p-value mapping for {motif_name}: {e}"
                    )
            else:
                logger.warning(
                    f"Motif {original_name} not found in PWM collection for p-value mapping"
                )

        logger.info(f"Loaded p-value mappings for {len(self._pvalue_mappings)} motifs")

    @property
    def pvalue_mappings(self) -> dict[str, tuple[int, np.ndarray, float]]:
        """Get p-value mappings for all motifs."""
        if self._pvalue_mappings is None:
            self._load_pvalue_mappings()
        return self._pvalue_mappings

    def get_score_threshold_for_pvalue(
        self, p_value: float, motif_indices: list[int] | None = None
    ) -> torch.Tensor:
        """
        Get score thresholds for a given p-value.

        Args:
            p_value: Desired p-value threshold (e.g., 0.0001)
            motif_indices: Optional list of motif indices to get thresholds for.
                          If None, returns thresholds for all motifs.

        Returns:
            Tensor of score thresholds with shape (n_motifs,)
        """
        mappings = self.pvalue_mappings
        motif_names = self.motif_names

        if motif_indices is None:
            motif_indices = list(range(len(motif_names)))

        thresholds = []
        log_pvalue = np.log2(p_value)

        for idx in motif_indices:
            if idx >= len(motif_names):
                thresholds.append(float("-inf"))
                continue

            motif_name = motif_names[idx]
            if motif_name not in mappings:
                # Fallback to old threshold system if available
                if (
                    self.significance_thresholds
                    and "p0.0001" in self.significance_thresholds
                ):
                    thresholds.append(
                        float(self.significance_thresholds["p0.0001"][idx])
                    )
                else:
                    thresholds.append(float("-inf"))
                continue

            smallest, logpdf, bin_size = mappings[motif_name]

            # Find the score threshold: find where logpdf <= log_pvalue
            # logpdf is in descending order (cumulative from right)
            threshold_idx = np.searchsorted(-logpdf, -log_pvalue, side="right")
            threshold_idx = max(0, min(threshold_idx, len(logpdf) - 1))

            # Convert index back to score
            score_threshold = (threshold_idx + smallest) * bin_size
            thresholds.append(score_threshold)

        return torch.tensor(thresholds, dtype=torch.float32)

    def scores_to_pvalues(
        self, scores: torch.Tensor, motif_indices: torch.Tensor | None = None
    ) -> torch.Tensor:
        """
        Convert motif scores to p-values.

        Args:
            scores: Tensor of scores with shape (n_positions,) or (n_motifs, n_positions)
                   or (batch, n_motifs, n_positions)
            motif_indices: Optional tensor of motif indices. If scores is 2D or 3D,
                          should have shape (n_motifs,) indicating which motif each
                          dimension corresponds to.

        Returns:
            Tensor of p-values with same shape as scores
        """
        mappings = self.pvalue_mappings
        motif_names = self.motif_names

        scores_np = scores.cpu().numpy()
        original_shape = scores_np.shape

        # Handle different input shapes
        if scores_np.ndim == 1:
            # Single motif, single position or single motif multiple positions
            if motif_indices is None:
                motif_idx = 0
            else:
                motif_idx = int(
                    motif_indices.item()
                    if motif_indices.numel() == 1
                    else motif_indices[0]
                )

            motif_name = motif_names[motif_idx]
            if motif_name in mappings:
                smallest, logpdf, bin_size = mappings[motif_name]
                pvalues_np = pvalue_from_logpdf(scores_np, smallest, logpdf, bin_size)
            else:
                pvalues_np = np.full_like(scores_np, 1.0)  # Default to non-significant
        else:
            # Multi-dimensional case
            if scores_np.ndim == 2:
                n_motifs, n_positions = scores_np.shape
                pvalues_np = np.zeros_like(scores_np)

                for i in range(n_motifs):
                    if motif_indices is not None:
                        motif_idx = int(motif_indices[i].item())
                    else:
                        motif_idx = i

                    motif_name = (
                        motif_names[motif_idx] if motif_idx < len(motif_names) else None
                    )
                    if motif_name and motif_name in mappings:
                        smallest, logpdf, bin_size = mappings[motif_name]
                        pvalues_np[i] = pvalue_from_logpdf(
                            scores_np[i], smallest, logpdf, bin_size
                        )
                    else:
                        pvalues_np[i] = 1.0
            else:  # ndim == 3
                batch_size, n_motifs, n_positions = scores_np.shape
                pvalues_np = np.zeros_like(scores_np)

                for b in range(batch_size):
                    for i in range(n_motifs):
                        if motif_indices is not None:
                            motif_idx = int(motif_indices[i].item())
                        else:
                            motif_idx = i

                        motif_name = (
                            motif_names[motif_idx]
                            if motif_idx < len(motif_names)
                            else None
                        )
                        if motif_name and motif_name in mappings:
                            smallest, logpdf, bin_size = mappings[motif_name]
                            pvalues_np[b, i] = pvalue_from_logpdf(
                                scores_np[b, i], smallest, logpdf, bin_size
                            )
                        else:
                            pvalues_np[b, i] = 1.0

        return torch.from_numpy(pvalues_np).to(scores.device).to(scores.dtype)

    def __repr__(self) -> str:
        """String representation."""
        source = "aligned" if self._loaded_with_aligned else "PWM source"
        return f"HocomocoIO({len(self)} motifs from {source})"
