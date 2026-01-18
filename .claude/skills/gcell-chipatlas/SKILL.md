---
name: gcell-chipatlas
description: |
  ChIP-Atlas epigenome data query using gcell. Use this skill when users ask about:
  - Finding ChIP-seq, ATAC-seq, or DNase-seq experiments
  - Searching for transcription factor binding data
  - Finding histone modification data (H3K27ac, H3K4me3, etc.)
  - Querying epigenome data by cell type or tissue
  - Downloading peak files or BigWig coverage data
  - Enrichment analysis for genomic regions or gene lists
  Triggers: ChIP-seq, ChIP-Atlas, ATAC-seq, DNase-seq, histone modification, H3K27ac, H3K4me3, transcription factor binding, epigenome, peak data, BigWig, enrichment analysis
refs:
  - _refs/gencode-gene-api.md
---

# ChIP-Atlas Epigenome Data Query

ChIP-Atlas integrates publicly available ChIP-seq, ATAC-seq, DNase-seq data from NCBI SRA.

## Helper Functions (Recommended)

Use helpers from this skill for common tasks:

```python
import sys
sys.path.insert(0, ".claude/skills/gcell-chipatlas")
from helpers import (
    check_binding_at_gene,
    get_promoter_region,
    find_celltype,
    signal_at_genes,
    plot_signal_at_gene,
    plot_signal_at_region,
)

# Plot ChIP-seq signal at a gene locus (easiest way to visualize)
track = plot_signal_at_gene("TERT", "CTCF", "K-562", out_file="ctcf_tert.png")
track = plot_signal_at_gene("ESR1", "ESR1", "MCF-7", upstream=20000, downstream=20000, out_file="esr1.png")

# Plot at arbitrary coordinates
track = plot_signal_at_region("chr17", 7670000, 7690000, "H3K27ac", "K-562", out_file="tp53_h3k27ac.png")

# Check if CTCF binds at TP53 promoter in K-562 cells
df = check_binding_at_gene("TP53", "CTCF", cell_type="K-562")
print(df[df['has_binding']])

# Get promoter coordinates (strand-aware)
chrom, start, end, strand = get_promoter_region("MYC", upstream=2000, downstream=500)

# Find correct cell type name
find_celltype("K562")  # Returns "K-562"

# Compare signal across multiple genes
df = signal_at_genes(["TP53", "MYC", "BRCA1"], "H3K27ac", cell_type="K-562")
```

### Function Signatures

```python
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
) -> Track:
    """Plot ChIP-seq signal tracks around a gene."""

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
) -> Track:
    """Plot ChIP-seq signal tracks at a genomic region."""

def check_binding_at_gene(
    gene_name: str,
    antigen: str,
    cell_type: str | None = None,
    assembly: str = "hg38",
    upstream: int = 2000,
    downstream: int = 500,
    threshold: float = 5.0,
    max_experiments: int = 10,
) -> pd.DataFrame:
    """Check if an antigen binds at a gene's promoter.
    Returns DataFrame with columns: experiment_id, max_signal, has_binding"""

def get_promoter_region(
    gene_name: str,
    assembly: str = "hg38",
    upstream: int = 2000,
    downstream: int = 500,
) -> tuple[str, int, int, str]:
    """Get promoter region coordinates (strand-aware).
    Returns (chrom, start, end, strand)"""

def find_celltype(
    query: str,
    assembly: str = "hg38",
) -> pd.DataFrame:
    """Find cell types matching a query string (handles variations like K562 -> K-562)."""

def signal_at_genes(
    gene_names: list[str],
    antigen: str,
    cell_type: str | None = None,
    assembly: str = "hg38",
    upstream: int = 2000,
    downstream: int = 500,
    max_experiments: int = 5,
) -> pd.DataFrame:
    """Get signal summary at promoters of multiple genes.
    Returns DataFrame with columns: gene, chrom, tss, strand, mean_signal, max_signal"""
```

## Core API

### Initialize and Search

```python
from gcell.epigenome import ChipAtlas

ca = ChipAtlas(metadata_mode="full")  # Cached SQLite DB

# Search experiments (sorted by library size by default - largest first)
df = ca.search(antigen="CTCF", cell_type="K-562", assembly="hg38")
df = ca.search(antigen="H3K27ac", assembly="hg38")  # Active enhancers

# Disable library size sorting if needed
df = ca.search(antigen="CTCF", assembly="hg38", sort_by_library_size=False)

# Browse available data
ca.get_antigens(assembly="hg38", antigen_class="TFs and others")
ca.get_celltypes(assembly="hg38")
```

### Stream Signal at Region

```python
# Stream BigWig data (no full download, ~3-5s per file for index fetch)
signals = ca.get_signal_at_region(
    experiment_ids=["SRX190161", "SRX190162"],
    chrom="chr17", start=7675594, end=7681594,
    assembly="hg38"
)
# Returns dict[exp_id, np.ndarray] at 1bp resolution
```

### Download Peaks

```python
# Get peaks (threshold: 5=Q<1E-05, 10=Q<1E-10, 20=Q<1E-20)
peaks = ca.get_peaks("SRX190161", threshold=5, as_pyranges=True)
```

### Track Visualization

**Preferred: Use helper functions** (see above) for one-liner plotting with gene annotations.

For manual control:

```python
from gcell.dna.track import Track
from gcell.epigenome.chipatlas.track_utils import stream_bigwig_regions_parallel

experiments = ca.search(antigen="H3K27ac", cell_type="K-562",
                        assembly="hg38", as_experiments=True, limit=5)

# Stream BigWig data in parallel
urls = [exp.bigwig_url for exp in experiments]
labels = [exp.experiment_id for exp in experiments]
signals = stream_bigwig_regions_parallel(urls, "chr8", 127735000, 127745000)

# Create and plot track
track = Track(chrom="chr8", start=127735000, end=127745000, assembly="hg38",
              tracks=dict(zip(labels, signals)))
track.plot_tracks_with_genebody(out_file="h3k27ac_myc.png", gene_annot=gencode, isoforms=True)
```

## Key Parameters

| Parameter | Values |
|-----------|--------|
| assembly | hg38, hg19, mm10, mm9, rn6, dm6, ce11, sacCer3 |
| antigen_class | "TFs and others", "Histone", "ATAC-Seq", "DNase-seq" |
| threshold | 5 (Q<1E-05), 10 (Q<1E-10), 20 (Q<1E-20) |

## Performance Notes

- **SQLite metadata**: Loads instantly after first download (~500MB)
- **BigWig streaming**: ~3-5s per file for index fetch, then ~0.5s per region
- **Cell type names**: Often hyphenated (e.g., "K-562" not "K562") - use `find_celltype()`
- **Library size sorting**: By default, searches return experiments sorted by read count (largest first) for better signal quality
