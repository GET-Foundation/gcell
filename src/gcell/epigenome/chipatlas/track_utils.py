"""Utilities for integrating ChIP-Atlas with gcell Track visualization.

This module provides functions to stream BigWig data from ChIP-Atlas
and create Track objects for visualization.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from pathlib import Path

    from gcell.dna.track import Track
    from gcell.epigenome.chipatlas.query import ChipAtlasExperiment


def stream_bigwig_region(
    url_or_path: str | Path,
    chrom: str,
    start: int,
    end: int,
) -> np.ndarray:
    """Stream a region from a BigWig file (local or remote URL).

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
        Signal values for the region (1bp resolution)

    Examples
    --------
    >>> # Stream from ChIP-Atlas URL
    >>> url = "https://chip-atlas.dbcls.jp/data/hg38/eachData/bw/SRX123456.bw"
    >>> signal = stream_bigwig_region(url, "chr1", 1000000, 1001000)
    """
    import pyBigWig

    bw = pyBigWig.open(str(url_or_path))
    try:
        values = bw.values(chrom, start, end)
        # Replace None/NaN with 0
        values = np.array(values, dtype=np.float32)
        values = np.nan_to_num(values, nan=0.0)
        return values
    finally:
        bw.close()


def experiments_to_track(
    experiments: list[ChipAtlasExperiment],
    chrom: str,
    start: int,
    end: int,
    assembly: str | None = None,
    use_cache: bool = True,
    labels: list[str] | None = None,
    normalizer=None,
    conv_size: int | None = 50,
) -> Track:
    """Create a Track object from ChIP-Atlas experiments.

    Downloads/streams BigWig data for multiple experiments and creates
    a Track object for visualization.

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
    use_cache : bool
        If True, download and cache BigWig files. If False, stream directly.
    labels : list[str], optional
        Custom labels for each experiment. If not provided, uses
        "{antigen}_{cell_type}" format.
    normalizer : Normalizer, optional
        Normalizer object to apply to tracks
    conv_size : int, optional
        Convolution window size. Default is 50.

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

    # Build track data
    tracks = {}
    for i, exp in enumerate(experiments):
        # Generate label
        if labels and i < len(labels):
            label = labels[i]
        else:
            label = f"{exp.antigen}_{exp.cell_type}"
            if not label.strip("_"):
                label = exp.experiment_id

        # Get BigWig data
        if use_cache:
            # Download and cache
            bw_path = exp.get_bigwig_path(cache=True)
            signal = stream_bigwig_region(bw_path, chrom, start, end)
        else:
            # Stream directly from URL
            signal = stream_bigwig_region(exp.bigwig_url, chrom, start, end)

        tracks[label] = signal

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
