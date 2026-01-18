"""Utilities for integrating ChIP-Atlas with gcell Track visualization.

This module provides functions to stream BigWig data from ChIP-Atlas
and create Track objects for visualization.
"""

from __future__ import annotations

import warnings
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from pathlib import Path

    from gcell.dna.track import Track
    from gcell.epigenome.chipatlas.query import ChipAtlasExperiment

# Threshold for auto-selecting streaming vs caching (100kb)
STREAMING_THRESHOLD_BP = 100_000


def _check_remote_support() -> bool:
    """Check if pyBigWig has remote URL support compiled in.

    Returns
    -------
    bool
        True if remote streaming is supported, False otherwise.
    """
    import pyBigWig

    return bool(getattr(pyBigWig, "remote", 0))


def stream_bigwig_region(
    url_or_path: str | Path,
    chrom: str,
    start: int,
    end: int,
) -> np.ndarray:
    """Stream a region from a BigWig file (local or remote URL).

    Uses pyBigWig's HTTP Range request support to fetch only the needed
    region from remote BigWig files, avoiding full file downloads.

    Parameters
    ----------
    url_or_path : str | Path
        Local path or remote URL to BigWig file
    chrom : str
        Chromosome name
    start : int
        Start position
    end : int
        End position

    Returns
    -------
    np.ndarray
        Signal values for the region (1bp resolution).
        Returns zeros if the file cannot be opened.

    Raises
    ------
    RuntimeError
        If URL streaming is requested but pyBigWig lacks remote support.

    Examples
    --------
    >>> # Stream from ChIP-Atlas URL
    >>> url = "https://chip-atlas.dbcls.jp/data/hg38/eachData/bw/SRX123456.bw"
    >>> signal = stream_bigwig_region(url, "chr1", 1000000, 1001000)
    """
    import pyBigWig

    path_str = str(url_or_path)
    is_remote = path_str.startswith(("http://", "https://"))

    # Check remote support for URLs
    if is_remote and not _check_remote_support():
        raise RuntimeError(
            "pyBigWig was compiled without remote URL support (pyBigWig.remote=0). "
            "To enable streaming, reinstall pyBigWig with curl development headers: "
            "1) Install libcurl-dev: apt-get install libcurl4-openssl-dev "
            "2) Reinstall pyBigWig: pip install --force-reinstall --no-binary pyBigWig pyBigWig"
        )

    try:
        bw = pyBigWig.open(path_str)
        if bw is None:
            # Return zeros if file cannot be opened
            return np.zeros(end - start, dtype=np.float32)
        try:
            values = bw.values(chrom, start, end)
            # Replace None/NaN with 0
            values = np.array(values, dtype=np.float32)
            values = np.nan_to_num(values, nan=0.0)
            return values
        finally:
            bw.close()
    except (RuntimeError, OSError):
        # Return zeros if file cannot be opened (e.g., ERX experiments)
        return np.zeros(end - start, dtype=np.float32)


def stream_bigwig_regions_parallel(
    urls_or_paths: list[str | Path],
    chrom: str,
    start: int,
    end: int,
    max_workers: int = 4,
) -> list[np.ndarray]:
    """Stream regions from multiple BigWig files in parallel.

    Parameters
    ----------
    urls_or_paths : list[str | Path]
        List of local paths or remote URLs to BigWig files
    chrom : str
        Chromosome name
    start : int
        Start position
    end : int
        End position
    max_workers : int
        Maximum number of parallel workers. Default is 4.

    Returns
    -------
    list[np.ndarray]
        List of signal arrays in the same order as input URLs
    """
    results = [None] * len(urls_or_paths)

    def fetch_one(idx: int, url_or_path: str | Path) -> tuple[int, np.ndarray]:
        signal = stream_bigwig_region(url_or_path, chrom, start, end)
        return idx, signal

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = [
            executor.submit(fetch_one, i, url)
            for i, url in enumerate(urls_or_paths)
        ]

        for future in as_completed(futures):
            idx, signal = future.result()
            results[idx] = signal

    return results


def experiments_to_track(
    experiments: list[ChipAtlasExperiment],
    chrom: str,
    start: int,
    end: int,
    assembly: str | None = None,
    use_cache: bool | None = None,
    labels: list[str] | None = None,
    normalizer=None,
    conv_size: int | None = 50,
    max_workers: int = 4,
) -> Track:
    """Create a Track object from ChIP-Atlas experiments.

    Downloads/streams BigWig data for multiple experiments and creates
    a Track object for visualization. Uses parallel fetching for speed.

    Parameters
    ----------
    experiments : list[ChipAtlasExperiment]
        List of ChipAtlasExperiment objects
    chrom : str
        Chromosome name
    start : int
        Start position
    end : int
        End position
    assembly : str, optional
        Genome assembly. If not provided, uses first experiment's assembly.
    use_cache : bool, optional
        If True, download and cache full BigWig files. If False, stream
        directly using HTTP Range requests. If None (default), automatically
        chooses: streaming for regions < 100kb (if pyBigWig.remote=1),
        otherwise downloads full files.
    labels : list[str], optional
        Custom labels for each experiment. If not provided, uses
        "{antigen}_{cell_type}" format.
    normalizer : Normalizer, optional
        Normalizer object to apply to tracks
    conv_size : int, optional
        Convolution window size. Default is 50.
    max_workers : int
        Maximum number of parallel workers for fetching. Default is 4.

    Returns
    -------
    Track
        Track object ready for visualization

    Examples
    --------
    >>> from gcell.epigenome import ChipAtlas
    >>> from gcell.epigenome.chipatlas.track_utils import experiments_to_track
    >>>
    >>> ca = ChipAtlas()
    >>> experiments = ca.search(
    ...     antigen="H3K27ac",
    ...     cell_type="K562",
    ...     assembly="hg38",
    ...     as_experiments=True,
    ...     limit=3,
    ... )
    >>>
    >>> track = experiments_to_track(
    ...     experiments,
    ...     chrom="chr1",
    ...     start=1000000,
    ...     end=1100000,
    ... )
    >>> track.plot_tracks()
    """
    from gcell.dna.track import Track

    if not experiments:
        raise ValueError("No experiments provided")

    if assembly is None:
        assembly = experiments[0].assembly

    # Check if remote streaming is supported
    has_remote = _check_remote_support()

    # Auto-select streaming vs caching based on region size and remote support
    region_size = end - start
    if use_cache is None:
        if not has_remote:
            # No remote support, must download
            use_cache = True
            warnings.warn(
                "pyBigWig.remote=0: Remote streaming not available. "
                "Downloading full BigWig files instead. "
                "To enable streaming, reinstall pyBigWig with curl support.",
                stacklevel=2,
            )
        else:
            # Use streaming for small regions (faster), caching for large regions
            use_cache = region_size > STREAMING_THRESHOLD_BP
    elif not use_cache and not has_remote:
        # User explicitly requested streaming but it's not available
        warnings.warn(
            "pyBigWig.remote=0: Cannot stream from URLs. "
            "Falling back to downloading full BigWig files. "
            "To enable streaming, reinstall pyBigWig with curl support.",
            stacklevel=2,
        )
        use_cache = True

    # Generate labels
    exp_labels = []
    for i, exp in enumerate(experiments):
        if labels and i < len(labels):
            label = labels[i]
        else:
            label = f"{exp.antigen}_{exp.cell_type}"
            if not label.strip("_"):
                label = exp.experiment_id
        exp_labels.append(label)

    # Get URLs or paths for all experiments
    if use_cache:
        # Download and cache files first (uses existing batch download if available)
        urls_or_paths = [exp.get_bigwig_path(cache=True) for exp in experiments]
    else:
        # Use remote URLs for streaming
        urls_or_paths = [exp.bigwig_url for exp in experiments]

    # Fetch all signals in parallel
    signals = stream_bigwig_regions_parallel(
        urls_or_paths, chrom, start, end, max_workers=max_workers
    )

    # Build tracks dict
    tracks = dict(zip(exp_labels, signals))

    return Track(
        chrom=chrom,
        start=start,
        end=end,
        assembly=assembly,
        tracks=tracks,
        normalizer=normalizer,
        conv_size=conv_size,
    )


def search_and_plot(
    antigen: str | None = None,
    cell_type: str | None = None,
    assembly: str = "hg38",
    chrom: str = "chr1",
    start: int = 1000000,
    end: int = 1100000,
    limit: int = 5,
    gene_annot=None,
    genes_to_highlight: list[str] | None = None,
    **search_kwargs,
):
    """Search ChIP-Atlas and plot tracks in one step.

    Convenience function that searches for experiments, downloads data,
    and creates a visualization.

    Parameters
    ----------
    antigen : str, optional
        Antigen to search for (e.g., "H3K27ac", "TP53")
    cell_type : str, optional
        Cell type to search for
    assembly : str
        Genome assembly. Default is "hg38".
    chrom : str
        Chromosome to plot. Default is "chr1".
    start : int
        Start position. Default is 1000000.
    end : int
        End position. Default is 1100000.
    limit : int
        Maximum number of experiments. Default is 5.
    gene_annot : Gencode, optional
        Gene annotation object for gene track
    genes_to_highlight : list[str], optional
        Genes to highlight in the plot
    **search_kwargs
        Additional arguments passed to ChipAtlas.search()

    Returns
    -------
    tuple
        (Track object, matplotlib figure, list of axes)

    Examples
    --------
    >>> from gcell.epigenome.chipatlas.track_utils import search_and_plot
    >>>
    >>> # Plot H3K27ac signal around MYC gene
    >>> track, fig, axes = search_and_plot(
    ...     antigen="H3K27ac",
    ...     cell_type="K562",
    ...     assembly="hg38",
    ...     chrom="chr8",
    ...     start=127735000,
    ...     end=127745000,
    ...     limit=3,
    ... )
    """
    from gcell.epigenome import ChipAtlas

    ca = ChipAtlas()

    # Search for experiments
    experiments = ca.search(
        antigen=antigen,
        cell_type=cell_type,
        assembly=assembly,
        limit=limit,
        as_experiments=True,
        **search_kwargs,
    )

    if not experiments:
        raise ValueError(
            f"No experiments found for antigen={antigen}, cell_type={cell_type}"
        )

    # Create track
    track = experiments_to_track(
        experiments,
        chrom=chrom,
        start=start,
        end=end,
        assembly=assembly,
    )

    # Plot
    fig, axes = track.plot_tracks(
        gene_annot=gene_annot,
        genes_to_highlight=genes_to_highlight,
    )

    return track, fig, axes


def compare_celltypes(
    antigen: str,
    cell_types: list[str],
    assembly: str = "hg38",
    chrom: str = "chr1",
    start: int = 1000000,
    end: int = 1100000,
    experiments_per_celltype: int = 1,
    gene_annot=None,
    genes_to_highlight: list[str] | None = None,
):
    """Compare the same antigen across different cell types.

    Parameters
    ----------
    antigen : str
        Antigen to compare (e.g., "H3K27ac")
    cell_types : list[str]
        List of cell types to compare
    assembly : str
        Genome assembly. Default is "hg38".
    chrom : str
        Chromosome to plot
    start : int
        Start position
    end : int
        End position
    experiments_per_celltype : int
        Number of experiments per cell type. Default is 1.
    gene_annot : Gencode, optional
        Gene annotation object
    genes_to_highlight : list[str], optional
        Genes to highlight

    Returns
    -------
    tuple
        (Track object, matplotlib figure, list of axes)

    Examples
    --------
    >>> from gcell.epigenome.chipatlas.track_utils import compare_celltypes
    >>>
    >>> track, fig, axes = compare_celltypes(
    ...     antigen="H3K27ac",
    ...     cell_types=["K562", "HepG2", "GM12878"],
    ...     chrom="chr8",
    ...     start=127735000,
    ...     end=127745000,
    ... )
    """
    from gcell.epigenome import ChipAtlas

    ca = ChipAtlas()
    all_experiments = []
    labels = []

    for ct in cell_types:
        experiments = ca.search(
            antigen=antigen,
            cell_type=ct,
            assembly=assembly,
            limit=experiments_per_celltype,
            as_experiments=True,
        )

        for exp in experiments:
            all_experiments.append(exp)
            labels.append(f"{ct}_{antigen}")

    if not all_experiments:
        raise ValueError(f"No experiments found for {antigen} in {cell_types}")

    # Create track with custom labels
    track = experiments_to_track(
        all_experiments,
        chrom=chrom,
        start=start,
        end=end,
        assembly=assembly,
        labels=labels,
    )

    # Plot
    fig, axes = track.plot_tracks(
        gene_annot=gene_annot,
        genes_to_highlight=genes_to_highlight,
    )

    return track, fig, axes


def peaks_to_bed_track(
    experiments: list[ChipAtlasExperiment],
    chrom: str,
    start: int,
    end: int,
    threshold: int = 5,
) -> np.ndarray:
    """Convert peak data from experiments to a binary BED track.

    Parameters
    ----------
    experiments : list[ChipAtlasExperiment]
        List of experiments
    chrom : str
        Chromosome name
    start : int
        Start position
    end : int
        End position
    threshold : int
        Peak threshold (5, 10, or 20)

    Returns
    -------
    np.ndarray
        Binary array with shape (end-start,) indicating peak regions
    """

    track = np.zeros(end - start, dtype=np.float32)

    for exp in experiments:
        try:
            peaks = exp.get_peaks(threshold=threshold)
            # Filter to region
            peaks = peaks[
                (peaks["Chromosome"] == chrom)
                & (peaks["End"] > start)
                & (peaks["Start"] < end)
            ]

            for _, peak in peaks.iterrows():
                peak_start = max(0, peak["Start"] - start)
                peak_end = min(end - start, peak["End"] - start)
                track[peak_start:peak_end] += 1

        except Exception:
            # Skip if peaks not available
            continue

    return track
