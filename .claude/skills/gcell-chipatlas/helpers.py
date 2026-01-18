"""Helper functions for ChIP-Atlas queries.

These functions wrap common analysis patterns to reduce boilerplate code.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    import pandas as pd
    from gcell.dna.track import Track
    from gcell.rna.gene import Gene


def get_promoter_region(
    gene_name: str,
    assembly: str = "hg38",
    upstream: int = 2000,
    downstream: int = 500,
) -> tuple[str, int, int, str]:
    """Get promoter region coordinates for a gene (strand-aware).

    Parameters
    ----------
    gene_name : str
        Gene symbol (e.g., "TP53", "MYC")
    assembly : str
        Genome assembly ("hg38", "hg19", "mm10")
    upstream : int
        Base pairs upstream of TSS (default: 2000)
    downstream : int
        Base pairs downstream of TSS (default: 500)

    Returns
    -------
    tuple[str, int, int, str]
        (chrom, start, end, strand)

    Example
    -------
    >>> chrom, start, end, strand = get_promoter_region("TP53")
    >>> print(f"{chrom}:{start}-{end} ({strand})")
    chr17:7685038-7690038 (-)
    """
    from gcell.rna.gencode import Gencode

    gencode = Gencode(assembly=assembly)
    gene = gencode.get_gene(gene_name)

    tss = gene.tss_coordinate
    if gene.strand == "-":
        # Minus strand: upstream is higher coordinates
        start = tss - downstream
        end = tss + upstream
    else:
        # Plus strand: upstream is lower coordinates
        start = tss - upstream
        end = tss + downstream

    return gene.chrom, start, end, gene.strand


def check_binding_at_gene(
    gene_name: str,
    antigen: str,
    cell_type: str | None = None,
    assembly: str = "hg38",
    upstream: int = 2000,
    downstream: int = 500,
    threshold: float = 5.0,
    max_experiments: int = 10,
) -> "pd.DataFrame":
    """Check if an antigen binds at a gene's promoter.

    Parameters
    ----------
    gene_name : str
        Gene symbol (e.g., "TP53", "MYC")
    antigen : str
        Antigen to search (e.g., "CTCF", "H3K27ac")
    cell_type : str, optional
        Cell type filter (e.g., "K-562")
    assembly : str
        Genome assembly
    upstream : int
        BP upstream of TSS to check
    downstream : int
        BP downstream of TSS to check
    threshold : float
        Signal threshold for "binding" call
    max_experiments : int
        Max experiments to query

    Returns
    -------
    pd.DataFrame
        Results with columns: experiment_id, max_signal, has_binding

    Example
    -------
    >>> df = check_binding_at_gene("TP53", "CTCF", cell_type="K-562")
    >>> print(df[df['has_binding']])
    """
    import pandas as pd

    from gcell.epigenome import ChipAtlas

    # Get promoter region
    chrom, start, end, strand = get_promoter_region(
        gene_name, assembly, upstream, downstream
    )

    # Search experiments
    ca = ChipAtlas(metadata_mode="full")
    search_results = ca.search(
        antigen=antigen, cell_type=cell_type, assembly=assembly
    )

    if len(search_results) == 0:
        return pd.DataFrame(columns=["experiment_id", "max_signal", "has_binding"])

    exp_ids = search_results["experiment_id"].head(max_experiments).tolist()

    # Stream signals
    signals = ca.get_signal_at_region(exp_ids, chrom, start, end, assembly=assembly)

    # Analyze binding
    results = []
    for exp_id, signal in signals.items():
        max_sig = float(np.max(signal)) if signal is not None else 0.0
        results.append({
            "experiment_id": exp_id,
            "max_signal": max_sig,
            "has_binding": max_sig > threshold,
        })

    return pd.DataFrame(results)


def find_celltype(
    query: str,
    assembly: str = "hg38",
) -> "pd.DataFrame":
    """Find cell types matching a query string.

    Handles common variations (e.g., "K562" matches "K-562").

    Parameters
    ----------
    query : str
        Partial cell type name (e.g., "K562", "HepG2")
    assembly : str
        Genome assembly

    Returns
    -------
    pd.DataFrame
        Matching cell types with experiment counts

    Example
    -------
    >>> find_celltype("K562")
       cell_type  experiment_count
    0      K-562              6900
    """
    import re

    from gcell.epigenome import ChipAtlas

    ca = ChipAtlas(metadata_mode="full")
    celltypes = ca.get_celltypes(assembly=assembly)

    # Build flexible regex: insert optional hyphen/space between any characters
    # "K562" -> "K[-\s]?562", "HepG2" -> "Hep[-\s]?G[-\s]?2"
    # This handles: K-562, K 562, K562, Hep G2, HepG2, etc.
    chars = list(query)
    pattern_parts = [chars[0]] if chars else []
    for c in chars[1:]:
        pattern_parts.append(r"[-\s]?")
        pattern_parts.append(c)
    pattern = "".join(pattern_parts)

    return celltypes[
        celltypes["cell_type"].str.contains(pattern, case=False, na=False, regex=True)
    ].reset_index(drop=True)


def plot_signal_at_region(
    chrom: str,
    start: int,
    end: int,
    antigen: str,
    cell_type: str | None = None,
    assembly: str = "hg38",
    max_experiments: int = 5,
    out_file: str | None = None,
    show_genes: bool = True,
    max_workers: int = 5,
) -> "Track":
    """Plot ChIP-seq signal tracks at a genomic region.

    Streams BigWig data in parallel and creates a multi-track visualization
    with optional gene annotations.

    Parameters
    ----------
    chrom : str
        Chromosome (e.g., "chr17")
    start : int
        Start coordinate
    end : int
        End coordinate
    antigen : str
        Antigen to search (e.g., "CTCF", "H3K27ac")
    cell_type : str, optional
        Cell type filter (e.g., "K-562")
    assembly : str
        Genome assembly ("hg38", "hg19", "mm10")
    max_experiments : int
        Maximum number of experiments to plot
    out_file : str, optional
        Output file path for the plot (e.g., "output.png")
    show_genes : bool
        Whether to show gene annotations (default: True)
    max_workers : int
        Number of parallel workers for BigWig streaming

    Returns
    -------
    Track
        Track object that can be further customized

    Example
    -------
    >>> # Plot CTCF at TERT locus +/- 50kb
    >>> chrom, start, end, _ = get_promoter_region("TERT", upstream=50000, downstream=50000)
    >>> track = plot_signal_at_region(chrom, start, end, "CTCF", "K-562", out_file="ctcf_tert.png")
    """
    from gcell.dna.track import Track
    from gcell.epigenome import ChipAtlas
    from gcell.epigenome.chipatlas.track_utils import stream_bigwig_regions_parallel

    ca = ChipAtlas(metadata_mode="full")

    # Search for experiments
    experiments = ca.search(
        antigen=antigen,
        cell_type=cell_type,
        assembly=assembly,
        as_experiments=True,
        limit=max_experiments,
    )

    if len(experiments) == 0:
        raise ValueError(f"No experiments found for {antigen} in {cell_type or 'any cell type'}")

    # Get BigWig URLs and labels
    urls = [exp.bigwig_url for exp in experiments]
    labels = [exp.experiment_id for exp in experiments]

    # Stream signals in parallel
    signals = stream_bigwig_regions_parallel(urls, chrom, start, end, max_workers=max_workers)
    tracks_dict = dict(zip(labels, signals))

    # Create Track object
    track = Track(chrom=chrom, start=start, end=end, assembly=assembly, tracks=tracks_dict)

    # Get gene annotations if requested
    gene_annot = None
    if show_genes:
        try:
            from gcell.rna.gencode import Gencode
            gene_annot = Gencode(assembly=assembly)
        except Exception:
            pass  # Skip gene annotations if unavailable

    # Plot with gene body visualization
    if out_file:
        track.plot_tracks_with_genebody(
            out_file=out_file,
            gene_annot=gene_annot,
            isoforms=True,
            gene_track_height=2,
        )

    return track


def plot_signal_at_gene(
    gene_name: str,
    antigen: str,
    cell_type: str | None = None,
    assembly: str = "hg38",
    upstream: int = 50000,
    downstream: int = 50000,
    max_experiments: int = 5,
    out_file: str | None = None,
    show_genes: bool = True,
    max_workers: int = 5,
) -> "Track":
    """Plot ChIP-seq signal tracks around a gene.

    Convenience wrapper around plot_signal_at_region that centers on a gene.

    Parameters
    ----------
    gene_name : str
        Gene symbol (e.g., "TP53", "MYC", "TERT")
    antigen : str
        Antigen to search (e.g., "CTCF", "H3K27ac")
    cell_type : str, optional
        Cell type filter (e.g., "K-562")
    assembly : str
        Genome assembly
    upstream : int
        Base pairs upstream of TSS (default: 50000)
    downstream : int
        Base pairs downstream of TSS (default: 50000)
    max_experiments : int
        Maximum number of experiments to plot
    out_file : str, optional
        Output file path for the plot
    show_genes : bool
        Whether to show gene annotations
    max_workers : int
        Number of parallel workers for BigWig streaming

    Returns
    -------
    Track
        Track object that can be further customized

    Example
    -------
    >>> track = plot_signal_at_gene("TERT", "CTCF", "K-562", out_file="ctcf_tert.png")
    """
    chrom, start, end, strand = get_promoter_region(
        gene_name, assembly, upstream, downstream
    )

    return plot_signal_at_region(
        chrom=chrom,
        start=start,
        end=end,
        antigen=antigen,
        cell_type=cell_type,
        assembly=assembly,
        max_experiments=max_experiments,
        out_file=out_file,
        show_genes=show_genes,
        max_workers=max_workers,
    )


def signal_at_genes(
    gene_names: list[str],
    antigen: str,
    cell_type: str | None = None,
    assembly: str = "hg38",
    upstream: int = 2000,
    downstream: int = 500,
    max_experiments: int = 5,
) -> "pd.DataFrame":
    """Get signal summary at promoters of multiple genes.

    Parameters
    ----------
    gene_names : list[str]
        List of gene symbols
    antigen : str
        Antigen to query
    cell_type : str, optional
        Cell type filter
    assembly : str
        Genome assembly
    upstream : int
        BP upstream of TSS
    downstream : int
        BP downstream of TSS
    max_experiments : int
        Max experiments to average over

    Returns
    -------
    pd.DataFrame
        Summary with columns: gene, chrom, tss, mean_signal, max_signal

    Example
    -------
    >>> df = signal_at_genes(["TP53", "MYC", "BRCA1"], "H3K27ac", "K-562")
    >>> print(df.sort_values("max_signal", ascending=False))
    """
    import pandas as pd

    from gcell.epigenome import ChipAtlas
    from gcell.rna.gencode import Gencode

    ca = ChipAtlas(metadata_mode="full")
    gencode = Gencode(assembly=assembly)

    # Get experiment IDs
    search_results = ca.search(
        antigen=antigen, cell_type=cell_type, assembly=assembly
    )
    if len(search_results) == 0:
        return pd.DataFrame()

    exp_ids = search_results["experiment_id"].head(max_experiments).tolist()

    results = []
    for gene_name in gene_names:
        try:
            gene = gencode.get_gene(gene_name)
            chrom, start, end, strand = get_promoter_region(
                gene_name, assembly, upstream, downstream
            )

            signals = ca.get_signal_at_region(
                exp_ids, chrom, start, end, assembly=assembly
            )

            # Average across experiments
            all_signals = [s for s in signals.values() if s is not None]
            if all_signals:
                combined = np.mean(all_signals, axis=0)
                results.append({
                    "gene": gene_name,
                    "chrom": chrom,
                    "tss": gene.tss_coordinate,
                    "strand": strand,
                    "mean_signal": float(np.mean(combined)),
                    "max_signal": float(np.max(combined)),
                })
        except Exception:
            continue

    return pd.DataFrame(results)
