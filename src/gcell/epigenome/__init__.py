"""Epigenome analysis module for gcell.

This module provides tools for working with epigenomic data including
ChIP-seq, ATAC-seq, DNase-seq, and bisulfite sequencing data from
public databases like ChIP-Atlas.
"""

from gcell.epigenome.chipatlas import (
    ChipAtlas,
    ChipAtlasEnrichment,
    ChipAtlasExperiment,
    EnrichmentResult,
)

__all__ = [
    "ChipAtlas",
    "ChipAtlasExperiment",
    "ChipAtlasEnrichment",
    "EnrichmentResult",
]
