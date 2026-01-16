"""ChIP-Atlas query interface.

This module provides the main ChipAtlas class for querying and downloading
ChIP-Atlas data.
"""

from __future__ import annotations

import io
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Literal

import pandas as pd

from gcell.epigenome.chipatlas.metadata import (
    ANTIGEN_CLASSES,
    BASE_URL,
    SUPPORTED_ASSEMBLIES,
    ChipAtlasMetadata,
    get_chipatlas_dir,
)

if TYPE_CHECKING:
    from collections.abc import Sequence

# Threshold values for peak calling (Q-value thresholds)
THRESHOLD_MAP = {
    5: "05",  # Q < 1E-05
    10: "10",  # Q < 1E-10
    20: "20",  # Q < 1E-20
}


@dataclass
class ChipAtlasExperiment:
    """Represents a ChIP-Atlas experiment.

    Attributes
    ----------
    experiment_id : str
        Experiment identifier (SRX, ERX, or DRX format)
    assembly : str
        Genome assembly (hg38, hg19, mm10, etc.)
    antigen : str
        Target antigen/protein name
    antigen_class : str
        Antigen class (TFs and others, Histone, etc.)
    cell_type : str
        Cell type name
    cell_type_class : str
        Cell type class (Blood, Breast, etc.)
    title : str
        Experiment title from GEO/SRA
    cell_type_description : str
        Detailed cell type description
    """

    experiment_id: str
    assembly: str
    antigen: str = ""
    antigen_class: str = ""
    cell_type: str = ""
    cell_type_class: str = ""
    title: str = ""
    cell_type_description: str = ""
    _chipatlas: ChipAtlas | None = field(default=None, repr=False)

    @property
    def bigwig_url(self) -> str:
        """URL for the BigWig coverage file."""
        return f"{BASE_URL}/{self.assembly}/eachData/bw/{self.experiment_id}.bw"

    def peak_url(self, threshold: int = 5) -> str:
        """Get URL for peak BED file.

        Parameters
        ----------
        threshold : int
            Peak calling threshold (5, 10, or 20)

        Returns
        -------
        str
            URL for the peak BED file
        """
        thresh_str = THRESHOLD_MAP.get(threshold, "05")
        return (
            f"{BASE_URL}/{self.assembly}/eachData/bed{thresh_str}/"
            f"{self.experiment_id}.{thresh_str}.bed"
        )

    def bigbed_url(self, threshold: int = 5) -> str:
        """Get URL for BigBed peak file.

        Parameters
        ----------
        threshold : int
            Peak calling threshold (5, 10, or 20)

        Returns
        -------
        str
            URL for the BigBed file
        """
        thresh_str = THRESHOLD_MAP.get(threshold, "05")
        return (
            f"{BASE_URL}/{self.assembly}/eachData/bb{thresh_str}/"
            f"{self.experiment_id}.{thresh_str}.bb"
        )

    def get_peaks(
        self,
        threshold: int = 5,
        as_pyranges: bool = False,
    ) -> pd.DataFrame:
        """Download and return peak data.

        Parameters
        ----------
        threshold : int
            Peak calling threshold (5, 10, or 20)
        as_pyranges : bool
            If True, return as PyRanges object

        Returns
        -------
        pd.DataFrame or PyRanges
            Peak data with columns: Chromosome, Start, End, Score
        """
        chipatlas = ChipAtlas() if self._chipatlas is None else self._chipatlas

        return chipatlas.get_peaks(
            self.experiment_id,
            assembly=self.assembly,
            threshold=threshold,
            as_pyranges=as_pyranges,
        )

    def get_bigwig_path(self, cache: bool = True) -> Path:
        """Download and return path to BigWig file.

        Parameters
        ----------
        cache : bool
            If True, cache the downloaded file

        Returns
        -------
        Path
            Path to the downloaded BigWig file
        """
        chipatlas = ChipAtlas() if self._chipatlas is None else self._chipatlas

        return chipatlas.download_bigwig(
            self.experiment_id,
            assembly=self.assembly,
            cache=cache,
        )


class ChipAtlas:
    """Main interface for querying ChIP-Atlas data.

    ChIP-Atlas is a comprehensive database integrating publicly available
    ChIP-seq, ATAC-seq, DNase-seq, and Bisulfite-seq data. This class provides
    efficient querying and downloading capabilities.

    Parameters
    ----------
    force_refresh : bool, optional
        If True, force re-download of metadata. Default is False.
    max_workers : int, optional
        Maximum number of parallel workers for downloads. Default is 4.
    cache_dir : str | Path, optional
        Directory for caching downloaded files. Default is the gcell cache.
    metadata_mode : {"full", "lite"}, optional
        Mode for metadata caching:
        - "full": Download complete metadata (~500MB). Slower first-time setup
          but comprehensive coverage of all experiments. Default.
        - "lite": Download partial metadata (~1MB). Quick setup for exploration
          but limited to ~5000 experiments per table.

    Examples
    --------
    >>> from gcell.epigenome import ChipAtlas

    Quick exploration with lite mode:

    >>> # Fast setup for exploration (~1MB download)
    >>> ca = ChipAtlas(metadata_mode="lite")
    >>> antigens = ca.get_antigens(assembly="hg38")

    >>> # Upgrade to full when ready
    >>> ca.upgrade_metadata()

    Full mode (default) for comprehensive queries:

    >>> # Complete metadata (~500MB, one-time download)
    >>> ca = ChipAtlas()  # or ChipAtlas(metadata_mode="full")

    Search for experiments:

    >>> # Find TP53 ChIP-seq experiments in human
    >>> experiments = ca.search(antigen="TP53", assembly="hg38")
    >>> print(f"Found {len(experiments)} experiments")

    >>> # Find all H3K27ac experiments in blood cells
    >>> experiments = ca.search(
    ...     antigen="H3K27ac",
    ...     assembly="hg38",
    ...     cell_type_class="Blood",
    ... )

    Download peak data:

    >>> # Get peaks for a specific experiment
    >>> peaks = ca.get_peaks("SRX123456", threshold=10)

    >>> # Download multiple experiments in parallel
    >>> peaks_dict = ca.get_peaks_batch(
    ...     ["SRX123456", "SRX789012"],
    ...     threshold=10,
    ... )

    Browse available data:

    >>> # Get all available antigens
    >>> antigens = ca.get_antigens(assembly="hg38")

    >>> # Get all cell types
    >>> celltypes = ca.get_celltypes(assembly="hg38")
    """

    def __init__(
        self,
        force_refresh: bool = False,
        max_workers: int = 4,
        cache_dir: str | Path | None = None,
        metadata_mode: str = "full",
        max_rows: int | None = None,
    ):
        self.metadata = ChipAtlasMetadata(
            force_refresh=force_refresh,
            max_workers=max_workers,
            metadata_mode=metadata_mode,
            max_rows=max_rows,
        )
        self.max_workers = max_workers
        self._cache_dir = Path(cache_dir) if cache_dir else get_chipatlas_dir()

    @property
    def cache_dir(self) -> Path:
        """Directory for caching downloaded files."""
        return self._cache_dir

    @property
    def supported_assemblies(self) -> list[str]:
        """List of supported genome assemblies."""
        return SUPPORTED_ASSEMBLIES.copy()

    @property
    def antigen_classes(self) -> list[str]:
        """List of antigen classes."""
        return ANTIGEN_CLASSES.copy()

    def search(
        self,
        assembly: str | None = None,
        antigen: str | None = None,
        antigen_class: str | None = None,
        cell_type: str | None = None,
        cell_type_class: str | None = None,
        title_contains: str | None = None,
        limit: int | None = None,
        as_experiments: bool = False,
    ) -> pd.DataFrame | list[ChipAtlasExperiment]:
        """Search for experiments matching criteria.

        Parameters
        ----------
        assembly : str, optional
            Genome assembly (hg38, hg19, mm10, etc.)
        antigen : str, optional
            Antigen/target name (e.g., TP53, H3K27ac). Partial match supported.
        antigen_class : str, optional
            Antigen class (TFs and others, Histone, etc.)
        cell_type : str, optional
            Cell type name. Partial match supported.
        cell_type_class : str, optional
            Cell type class (Blood, Breast, etc.)
        title_contains : str, optional
            Search for experiments with title containing this string
        limit : int, optional
            Maximum number of results to return
        as_experiments : bool, optional
            If True, return list of ChipAtlasExperiment objects

        Returns
        -------
        pd.DataFrame or list[ChipAtlasExperiment]
            Matching experiments

        Examples
        --------
        >>> ca = ChipAtlas()
        >>> # Find TP53 experiments
        >>> df = ca.search(antigen="TP53", assembly="hg38")
        >>> # Get as experiment objects for easy data access
        >>> experiments = ca.search(
        ...     antigen="TP53", assembly="hg38", as_experiments=True
        ... )
        >>> peaks = experiments[0].get_peaks()
        """
        df = self.metadata.search_experiments(
            assembly=assembly,
            antigen=antigen,
            antigen_class=antigen_class,
            cell_type=cell_type,
            cell_type_class=cell_type_class,
            title_contains=title_contains,
            limit=limit,
        )

        if as_experiments:
            return self._df_to_experiments(df)

        return df

    def _df_to_experiments(self, df: pd.DataFrame) -> list[ChipAtlasExperiment]:
        """Convert DataFrame to list of ChipAtlasExperiment objects."""
        experiments = []
        for _, row in df.iterrows():
            exp = ChipAtlasExperiment(
                experiment_id=row.get("experiment_id", ""),
                assembly=row.get("assembly", ""),
                antigen=row.get("antigen", ""),
                antigen_class=row.get("antigen_class", ""),
                cell_type=row.get("cell_type", ""),
                cell_type_class=row.get("cell_type_class", ""),
                title=row.get("title", ""),
                cell_type_description=row.get("cell_type_description", ""),
                _chipatlas=self,
            )
            experiments.append(exp)
        return experiments

    def get_experiment(self, experiment_id: str) -> ChipAtlasExperiment | None:
        """Get a specific experiment by ID.

        Parameters
        ----------
        experiment_id : str
            Experiment identifier (SRX, ERX, or DRX format)

        Returns
        -------
        ChipAtlasExperiment or None
            Experiment object or None if not found
        """
        df = self.metadata.search_experiments(limit=1)
        # Need to query directly
        import sqlite3

        from gcell.epigenome.chipatlas.metadata import get_db_path

        with sqlite3.connect(get_db_path()) as conn:
            conn.row_factory = sqlite3.Row
            cursor = conn.execute(
                "SELECT * FROM experiments WHERE experiment_id = ?",
                (experiment_id,),
            )
            row = cursor.fetchone()

        if row is None:
            return None

        return ChipAtlasExperiment(
            experiment_id=row["experiment_id"],
            assembly=row["assembly"],
            antigen=row["antigen"] or "",
            antigen_class=row["antigen_class"] or "",
            cell_type=row["cell_type"] or "",
            cell_type_class=row["cell_type_class"] or "",
            title=row["title"] or "",
            cell_type_description=row["cell_type_description"] or "",
            _chipatlas=self,
        )

    def get_peaks(
        self,
        experiment_id: str,
        assembly: str | None = None,
        threshold: int = 5,
        as_pyranges: bool = False,
        cache: bool = True,
    ) -> pd.DataFrame:
        """Download and return peak data for an experiment.

        Parameters
        ----------
        experiment_id : str
            Experiment identifier
        assembly : str, optional
            Genome assembly. If not provided, will be looked up.
        threshold : int
            Peak calling threshold (5, 10, or 20)
        as_pyranges : bool
            If True, return as PyRanges object
        cache : bool
            If True, cache the downloaded file

        Returns
        -------
        pd.DataFrame or PyRanges
            Peak data with columns: Chromosome, Start, End, Score

        Raises
        ------
        ValueError
            If experiment not found or invalid threshold
        """
        if threshold not in THRESHOLD_MAP:
            raise ValueError(
                f"Invalid threshold. Must be one of: {list(THRESHOLD_MAP)}"
            )

        # Get assembly if not provided
        if assembly is None:
            exp = self.get_experiment(experiment_id)
            if exp is None:
                raise ValueError(f"Experiment not found: {experiment_id}")
            assembly = exp.assembly

        # Build URL
        thresh_str = THRESHOLD_MAP[threshold]
        url = (
            f"{BASE_URL}/{assembly}/eachData/bed{thresh_str}/"
            f"{experiment_id}.{thresh_str}.bed"
        )

        # Check cache
        cache_path = self._cache_dir / "peaks" / f"{experiment_id}.{thresh_str}.bed"

        if cache and cache_path.exists():
            df = pd.read_csv(
                cache_path,
                sep="\t",
                header=None,
                names=["Chromosome", "Start", "End", "Score"],
            )
        else:
            # Download
            import urllib.request

            try:
                with urllib.request.urlopen(url, timeout=30) as response:
                    content = response.read().decode("utf-8")
            except Exception as e:
                raise RuntimeError(
                    f"Failed to download peaks for {experiment_id}: {e}"
                ) from e

            df = pd.read_csv(
                io.StringIO(content),
                sep="\t",
                header=None,
                names=["Chromosome", "Start", "End", "Score"],
            )

            # Cache if requested
            if cache:
                cache_path.parent.mkdir(parents=True, exist_ok=True)
                df.to_csv(cache_path, sep="\t", header=False, index=False)

        if as_pyranges:
            import pyranges as pr

            return pr.PyRanges(df)

        return df

    def get_peaks_batch(
        self,
        experiment_ids: Sequence[str],
        assembly: str | None = None,
        threshold: int = 5,
        max_workers: int | None = None,
        on_error: Literal["raise", "skip", "warn"] = "warn",
    ) -> dict[str, pd.DataFrame]:
        """Download peak data for multiple experiments in parallel.

        Parameters
        ----------
        experiment_ids : Sequence[str]
            List of experiment identifiers
        assembly : str, optional
            Genome assembly. If not provided, will be looked up for each.
        threshold : int
            Peak calling threshold (5, 10, or 20)
        max_workers : int, optional
            Maximum number of parallel workers. Default uses instance setting.
        on_error : str
            How to handle errors: "raise", "skip", or "warn"

        Returns
        -------
        dict[str, pd.DataFrame]
            Dictionary mapping experiment IDs to peak DataFrames

        Examples
        --------
        >>> ca = ChipAtlas()
        >>> peaks = ca.get_peaks_batch(["SRX123456", "SRX789012"])
        >>> for exp_id, df in peaks.items():
        ...     print(f"{exp_id}: {len(df)} peaks")
        """
        workers = max_workers or self.max_workers
        results = {}

        def download_one(
            exp_id: str,
        ) -> tuple[str, pd.DataFrame | None, Exception | None]:
            try:
                df = self.get_peaks(exp_id, assembly=assembly, threshold=threshold)
                return exp_id, df, None
            except Exception as e:
                return exp_id, None, e

        with ThreadPoolExecutor(max_workers=workers) as executor:
            futures = [
                executor.submit(download_one, exp_id) for exp_id in experiment_ids
            ]

            for future in as_completed(futures):
                exp_id, df, error = future.result()
                if error is not None:
                    if on_error == "raise":
                        raise error
                    elif on_error == "warn":
                        import warnings

                        warnings.warn(
                            f"Failed to download peaks for {exp_id}: {error}",
                            stacklevel=2,
                        )
                elif df is not None:
                    results[exp_id] = df

        return results

    def download_bigwig(
        self,
        experiment_id: str,
        assembly: str | None = None,
        cache: bool = True,
        output_path: str | Path | None = None,
    ) -> Path:
        """Download BigWig coverage file for an experiment.

        Parameters
        ----------
        experiment_id : str
            Experiment identifier
        assembly : str, optional
            Genome assembly. If not provided, will be looked up.
        cache : bool
            If True, cache the downloaded file
        output_path : str | Path, optional
            Custom output path. If not provided, uses cache directory.

        Returns
        -------
        Path
            Path to the downloaded BigWig file
        """
        # Get assembly if not provided
        if assembly is None:
            exp = self.get_experiment(experiment_id)
            if exp is None:
                raise ValueError(f"Experiment not found: {experiment_id}")
            assembly = exp.assembly

        # Build URL
        url = f"{BASE_URL}/{assembly}/eachData/bw/{experiment_id}.bw"

        # Determine output path
        if output_path:
            out_path = Path(output_path)
        else:
            out_path = self._cache_dir / "bigwig" / f"{experiment_id}.bw"

        # Check cache
        if cache and out_path.exists():
            return out_path

        # Download
        out_path.parent.mkdir(parents=True, exist_ok=True)

        import urllib.request

        try:
            urllib.request.urlretrieve(url, out_path)
        except Exception as e:
            raise RuntimeError(
                f"Failed to download BigWig for {experiment_id}: {e}"
            ) from e

        return out_path

    def download_bigwig_batch(
        self,
        experiment_ids: Sequence[str],
        assembly: str | None = None,
        max_workers: int | None = None,
        on_error: Literal["raise", "skip", "warn"] = "warn",
    ) -> dict[str, Path]:
        """Download BigWig files for multiple experiments in parallel.

        Parameters
        ----------
        experiment_ids : Sequence[str]
            List of experiment identifiers
        assembly : str, optional
            Genome assembly. If not provided, will be looked up for each.
        max_workers : int, optional
            Maximum number of parallel workers
        on_error : str
            How to handle errors: "raise", "skip", or "warn"

        Returns
        -------
        dict[str, Path]
            Dictionary mapping experiment IDs to BigWig file paths
        """
        workers = max_workers or self.max_workers
        results = {}

        def download_one(exp_id: str) -> tuple[str, Path | None, Exception | None]:
            try:
                path = self.download_bigwig(exp_id, assembly=assembly)
                return exp_id, path, None
            except Exception as e:
                return exp_id, None, e

        with ThreadPoolExecutor(max_workers=workers) as executor:
            futures = [
                executor.submit(download_one, exp_id) for exp_id in experiment_ids
            ]

            for future in as_completed(futures):
                exp_id, path, error = future.result()
                if error is not None:
                    if on_error == "raise":
                        raise error
                    elif on_error == "warn":
                        import warnings

                        warnings.warn(
                            f"Failed to download BigWig for {exp_id}: {error}",
                            stacklevel=2,
                        )
                elif path is not None:
                    results[exp_id] = path

        return results

    def get_antigens(
        self,
        assembly: str | None = None,
        antigen_class: str | None = None,
        min_experiments: int = 0,
    ) -> pd.DataFrame:
        """Get available antigens/targets.

        Parameters
        ----------
        assembly : str, optional
            Filter by genome assembly
        antigen_class : str, optional
            Filter by antigen class
        min_experiments : int, optional
            Minimum number of experiments required

        Returns
        -------
        pd.DataFrame
            DataFrame with antigen information
        """
        return self.metadata.get_antigens(
            assembly=assembly,
            antigen_class=antigen_class,
            min_experiments=min_experiments,
        )

    def get_celltypes(
        self,
        assembly: str | None = None,
        cell_type_class: str | None = None,
        min_experiments: int = 0,
    ) -> pd.DataFrame:
        """Get available cell types.

        Parameters
        ----------
        assembly : str, optional
            Filter by genome assembly
        cell_type_class : str, optional
            Filter by cell type class
        min_experiments : int, optional
            Minimum number of experiments required

        Returns
        -------
        pd.DataFrame
            DataFrame with cell type information
        """
        return self.metadata.get_celltypes(
            assembly=assembly,
            cell_type_class=cell_type_class,
            min_experiments=min_experiments,
        )

    def get_celltype_classes(self, assembly: str | None = None) -> list[str]:
        """Get available cell type classes.

        Parameters
        ----------
        assembly : str, optional
            Filter by genome assembly

        Returns
        -------
        list[str]
            List of cell type class names
        """
        return self.metadata.get_celltype_classes(assembly=assembly)

    def get_antigen_classes(self, assembly: str | None = None) -> list[str]:
        """Get available antigen classes.

        Parameters
        ----------
        assembly : str, optional
            Filter by genome assembly

        Returns
        -------
        list[str]
            List of antigen class names
        """
        return self.metadata.get_antigen_classes(assembly=assembly)

    def get_assemblies(self) -> list[str]:
        """Get available genome assemblies.

        Returns
        -------
        list[str]
            List of genome assembly names
        """
        return self.metadata.get_assemblies()

    def summary(self, assembly: str | None = None) -> pd.DataFrame:
        """Get summary statistics of available data.

        Parameters
        ----------
        assembly : str, optional
            Filter by genome assembly

        Returns
        -------
        pd.DataFrame
            Summary statistics by assembly and antigen class
        """
        import sqlite3

        from gcell.epigenome.chipatlas.metadata import get_db_path

        query = """
            SELECT
                assembly,
                antigen_class,
                COUNT(*) as experiment_count,
                COUNT(DISTINCT antigen) as unique_antigens,
                COUNT(DISTINCT cell_type) as unique_celltypes
            FROM experiments
        """
        params = []

        if assembly:
            query += " WHERE assembly = ?"
            params.append(assembly)

        query += " GROUP BY assembly, antigen_class ORDER BY assembly, antigen_class"

        with sqlite3.connect(get_db_path()) as conn:
            df = pd.read_sql_query(query, conn, params=params)

        return df

    def refresh_metadata(self) -> None:
        """Force refresh of the metadata database."""
        self.metadata.refresh()

    def upgrade_metadata(self) -> None:
        """Upgrade from lite mode to full metadata.

        Downloads the complete metadata if currently in lite mode.
        This is useful when you started with lite mode for quick exploration
        and now need comprehensive coverage.

        Examples
        --------
        >>> ca = ChipAtlas(metadata_mode="lite")  # Quick start
        >>> # ... explore data ...
        >>> ca.upgrade_metadata()  # Download full metadata when ready
        """
        self.metadata.upgrade_to_full()

    @property
    def metadata_mode(self) -> str:
        """Current metadata mode ('full', 'lite', or 'unknown')."""
        return self.metadata.get_mode()

    @property
    def is_lite_mode(self) -> bool:
        """Check if using lite mode (partial metadata)."""
        return self.metadata.is_lite_mode()
