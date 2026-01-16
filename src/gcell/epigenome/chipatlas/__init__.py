"""ChIP-Atlas query module for gcell.

ChIP-Atlas is a comprehensive database of public ChIP-seq, ATAC-seq,
DNase-seq, and Bisulfite-seq data. This module provides efficient
querying and downloading capabilities.

Example
-------
>>> from gcell.epigenome import ChipAtlas
>>> ca = ChipAtlas()
>>> # Search for TP53 ChIP-seq experiments in human
>>> experiments = ca.search(antigen="TP53", assembly="hg38")
>>> # Download peak data for specific experiment
>>> peaks = ca.get_peaks("SRX123456", threshold=10)

References
----------
ChIP-Atlas: https://chip-atlas.org/
Publication: Oki et al., EMBO Reports (2018)
"""

from gcell.epigenome.chipatlas.enrichment import (
    ChipAtlasEnrichment,
    EnrichmentResult,
)
from gcell.epigenome.chipatlas.query import ChipAtlas, ChipAtlasExperiment

__all__ = [
    "ChipAtlas",
    "ChipAtlasExperiment",
    "ChipAtlasEnrichment",
    "EnrichmentResult",
]
