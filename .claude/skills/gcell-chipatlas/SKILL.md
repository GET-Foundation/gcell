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
---

# ChIP-Atlas Epigenome Data Query

ChIP-Atlas is a comprehensive database integrating publicly available ChIP-seq, ATAC-seq, DNase-seq, and Bisulfite-seq data from NCBI SRA.

## Quick Start

```python
from gcell.epigenome import ChipAtlas

# Initialize ChipAtlas (downloads metadata on first use)
ca = ChipAtlas()

# Search for experiments
experiments = ca.search(antigen="TP53", assembly="hg38")
print(f"Found {len(experiments)} TP53 experiments")
```

## Metadata Modes

On first use, ChipAtlas downloads metadata to enable fast local queries. Choose the mode that fits your needs:

### Lite Mode (Quick Setup)

```python
# Fast setup (~1MB download) for exploration
ca = ChipAtlas(metadata_mode="lite")

# Limited to ~5000 experiments per table, but instant setup
antigens = ca.get_antigens(assembly="hg38")

# Check current mode
print(f"Mode: {ca.metadata_mode}")  # "lite"

# Upgrade to full when ready for comprehensive queries
ca.upgrade_metadata()
```

### Full Mode (Default, Comprehensive)

```python
# Complete metadata (~500MB, one-time download)
ca = ChipAtlas()  # or ChipAtlas(metadata_mode="full")

# Covers all experiments in ChIP-Atlas
# Takes a few minutes on first run, but only needs to download once
```

### When to Use Each Mode

| Mode | Download Size | Use Case |
|------|--------------|----------|
| `lite` | ~1MB | Quick exploration, testing, demos |
| `full` | ~500MB | Production queries, comprehensive searches |

### Refresh Metadata

```python
# Force re-download (e.g., to get latest experiments)
ca = ChipAtlas(force_refresh=True)

# Or refresh existing instance
ca.refresh_metadata()
```

## Searching for Experiments

### By Transcription Factor / Antigen

```python
# Find TP53 ChIP-seq experiments
df = ca.search(antigen="TP53", assembly="hg38")

# Find CTCF experiments in mouse
df = ca.search(antigen="CTCF", assembly="mm10")

# Find all experiments for a TF
df = ca.search(antigen="STAT3", assembly="hg38")
```

### By Histone Modification

```python
# H3K27ac (active enhancers)
df = ca.search(antigen="H3K27ac", assembly="hg38")

# H3K4me3 (active promoters)
df = ca.search(antigen="H3K4me3", assembly="hg38")

# H3K27me3 (repressive mark)
df = ca.search(antigen="H3K27me3", assembly="hg38")

# H3K4me1 (poised enhancers)
df = ca.search(antigen="H3K4me1", assembly="hg38")
```

### By Cell Type

```python
# Find experiments in a specific cell type
df = ca.search(cell_type="HeLa", assembly="hg38")

# Find experiments in a cell type class
df = ca.search(cell_type_class="Blood", assembly="hg38")
df = ca.search(cell_type_class="Pluripotent stem cell", assembly="hg38")
df = ca.search(cell_type_class="Breast", assembly="hg38")
```

### Combined Searches

```python
# H3K27ac in blood cells
df = ca.search(
    antigen="H3K27ac",
    cell_type_class="Blood",
    assembly="hg38"
)

# TP53 in cancer cell lines
df = ca.search(
    antigen="TP53",
    cell_type="cancer",  # partial match
    assembly="hg38"
)
```

### By Experiment Type

```python
# ATAC-seq experiments
df = ca.search(antigen_class="ATAC-Seq", assembly="hg38")

# DNase-seq experiments
df = ca.search(antigen_class="DNase-seq", assembly="hg38")

# All TF ChIP-seq
df = ca.search(antigen_class="TFs and others", assembly="hg38")

# All histone ChIP-seq
df = ca.search(antigen_class="Histone", assembly="hg38")
```

## Downloading Peak Data

```python
# Get peaks for a specific experiment
peaks = ca.get_peaks("SRX190161", threshold=5)
# threshold: 5 (Q<1E-05), 10 (Q<1E-10), 20 (Q<1E-20)

# Returns DataFrame with columns: Chromosome, Start, End, Score
print(peaks.head())

# Get as PyRanges object for interval operations
peaks_pr = ca.get_peaks("SRX190161", as_pyranges=True)
```

### Batch Download

```python
# Download peaks for multiple experiments in parallel
exp_ids = ["SRX190161", "SRX190162", "SRX190163"]
peaks_dict = ca.get_peaks_batch(exp_ids, threshold=10)

for exp_id, df in peaks_dict.items():
    print(f"{exp_id}: {len(df)} peaks")
```

## Downloading BigWig Files

```python
# Download BigWig coverage file
bw_path = ca.download_bigwig("SRX190161")
print(f"Downloaded to: {bw_path}")

# Batch download
paths = ca.download_bigwig_batch(exp_ids)
```

## Using Experiment Objects

```python
# Get experiments as objects for easy data access
experiments = ca.search(
    antigen="TP53",
    assembly="hg38",
    as_experiments=True,
    limit=10
)

# Access data directly from experiment
exp = experiments[0]
print(f"Experiment: {exp.experiment_id}")
print(f"Cell type: {exp.cell_type}")
print(f"Title: {exp.title}")

# Download data
peaks = exp.get_peaks(threshold=10)
bw_path = exp.get_bigwig_path()
```

## Browsing Available Data

```python
# List available genome assemblies
assemblies = ca.get_assemblies()
print(assemblies)  # ['hg38', 'hg19', 'mm10', ...]

# List available antigens/targets
antigens = ca.get_antigens(assembly="hg38")
print(antigens.head(20))

# List antigens by class
tf_antigens = ca.get_antigens(
    assembly="hg38",
    antigen_class="TFs and others",
    min_experiments=10  # At least 10 experiments
)

# List available cell types
celltypes = ca.get_celltypes(assembly="hg38")
print(celltypes.head(20))

# List cell type classes
classes = ca.get_celltype_classes(assembly="hg38")
print(classes)

# Get summary statistics
summary = ca.summary(assembly="hg38")
print(summary)
```

## Enrichment Analysis

Find TFs or histone marks enriched in your regions of interest.

```python
from gcell.epigenome import ChipAtlasEnrichment

enrichment = ChipAtlasEnrichment()

# Analyze genomic regions (BED format)
bed_data = """chr1\t1000\t2000
chr1\t3000\t4000
chr2\t5000\t6000"""

result = enrichment.analyze_regions(
    bed_data,
    assembly="hg38",
    antigen_class="TFs and others",
)

if result.is_successful:
    print(result.results.head())
```

### Analyze Gene List

```python
# Find TFs enriched near gene promoters
genes = ["TP53", "MYC", "BRCA1", "PTEN", "RB1"]

result = enrichment.analyze_genes(
    genes,
    assembly="hg38",
    antigen_class="TFs and others",
    distance_up=5000,    # 5kb upstream of TSS
    distance_down=1000,  # 1kb downstream of TSS
)

if result.is_successful:
    # Results show enriched TFs
    print(result.results)
```

## Supported Assemblies

| Assembly | Species |
|----------|---------|
| hg38 | Human (GRCh38) |
| hg19 | Human (GRCh37) |
| mm10 | Mouse (GRCm38) |
| mm9 | Mouse (NCBI37) |
| rn6 | Rat |
| dm6 | Drosophila |
| dm3 | Drosophila (legacy) |
| ce11 | C. elegans |
| ce10 | C. elegans (legacy) |
| sacCer3 | Yeast |

## Antigen Classes

| Class | Description |
|-------|-------------|
| TFs and others | Transcription factors |
| Histone | Histone modifications (H3K27ac, H3K4me3, etc.) |
| RNA polymerase | RNA Pol II, etc. |
| ATAC-Seq | Chromatin accessibility |
| DNase-seq | DNase hypersensitivity |
| Bisulfite-Seq | DNA methylation |
| Input control | Control experiments |

## Data Location

- Metadata: `~/.gcell_data/annotations/chipatlas/chipatlas.db`
- Cached peaks: `~/.gcell_data/annotations/chipatlas/peaks/`
- Cached BigWig: `~/.gcell_data/annotations/chipatlas/bigwig/`

## Common Workflows

### Find Regulatory Elements for a Gene

```python
from gcell.epigenome import ChipAtlas

ca = ChipAtlas()

# 1. Find TF binding experiments
tf_exps = ca.search(
    antigen_class="TFs and others",
    cell_type="K562",
    assembly="hg38",
    limit=50
)

# 2. Download peaks
exp_ids = tf_exps['experiment_id'].tolist()
all_peaks = ca.get_peaks_batch(exp_ids)

# 3. Check which TFs have peaks near your gene of interest
# (Combine with gcell-gene skill for TSS coordinates)
```

### Compare Enhancer Activity Across Cell Types

```python
# Get H3K27ac data for different cell types
cell_types = ["K562", "HepG2", "GM12878"]

for ct in cell_types:
    exps = ca.search(
        antigen="H3K27ac",
        cell_type=ct,
        assembly="hg38",
        limit=5
    )
    print(f"{ct}: {len(exps)} experiments")
```

## Track Visualization

ChIP-Atlas integrates with gcell's Track class for BigWig visualization.

```python
from gcell.dna.track import Track
```

The Track class provides three plotting methods:
- `plot_tracks()` - Basic track visualization
- `plot_tracks_with_motif_density()` - Tracks with motif annotations
- `plot_tracks_with_genebody()` - Tracks with detailed gene structure visualization

### Manual Track Construction from BigWig

The most flexible approach is to manually load BigWig data and construct a Track object:

```python
from gcell.epigenome import ChipAtlas
from gcell.dna.track import Track  # Basic Track class
import pyBigWig
import numpy as np

ca = ChipAtlas()

# Download BigWig file
exp_id = "SRX3322103"
bw_path = ca.download_bigwig(exp_id)

# Load data from BigWig
bw = pyBigWig.open(str(bw_path))
values = bw.values("chr5", 1293068, 1297068)
values = np.nan_to_num(np.array(values), nan=0)
bw.close()

# Create Track object - tracks parameter is a dict of {name: array}
track = Track(
    chrom="chr5",
    start=1293068,
    end=1297068,
    assembly="hg38",
    tracks={"CTCF_Epithelial": values},  # Can add multiple tracks
    conv_size=50  # Smoothing window size
)

# Plot
track.plot_tracks(out_file="output.png")
```

### Multiple Tracks from Different Experiments

```python
from gcell.epigenome import ChipAtlas
from gcell.dna.track import Track
import pyBigWig
import numpy as np

ca = ChipAtlas()

# Define region
chrom, start, end = "chr5", 1293068, 1297068

# Load multiple experiments
experiments = {
    "SRX3322103": "Prostate_Epithelial",
    "SRX10142196": "Mammary_Epithelial_1",
    "SRX10142197": "Mammary_Epithelial_2",
}

tracks_dict = {}
for exp_id, name in experiments.items():
    bw_path = ca.download_bigwig(exp_id)
    bw = pyBigWig.open(str(bw_path))
    values = bw.values(chrom, start, end)
    values = np.nan_to_num(np.array(values), nan=0)
    bw.close()
    tracks_dict[name] = values

# Create multi-track visualization
track = Track(
    chrom=chrom,
    start=start,
    end=end,
    assembly="hg38",
    tracks=tracks_dict,
    conv_size=50
)

track.plot_tracks(out_file="multi_tracks.png")
```

### Stream BigWig Data

```python
from gcell.epigenome.chipatlas import stream_bigwig_region

# Stream directly from ChIP-Atlas URL (no download needed)
url = "https://chip-atlas.dbcls.jp/data/hg38/eachData/bw/SRX190161.bw"
signal = stream_bigwig_region(url, "chr1", 1000000, 1100000)
print(f"Signal shape: {signal.shape}")
```

### Create Track from Experiments

```python
from gcell.epigenome import ChipAtlas
from gcell.epigenome.chipatlas import experiments_to_track

ca = ChipAtlas()

# Search for experiments
experiments = ca.search(
    antigen="H3K27ac",
    cell_type="K562",
    assembly="hg38",
    as_experiments=True,
    limit=3
)

# Create Track object for visualization
track = experiments_to_track(
    experiments,
    chrom="chr8",
    start=127735000,
    end=127745000,  # MYC locus
)

# Plot tracks
track.plot_tracks()
```

### Quick Search and Plot

```python
from gcell.epigenome.chipatlas import search_and_plot

# Search and plot in one step
track, fig, axes = search_and_plot(
    antigen="H3K27ac",
    cell_type="K562",
    assembly="hg38",
    chrom="chr8",
    start=127735000,
    end=127745000,
    limit=3,
)
```

### Compare Cell Types with Track Visualization

```python
from gcell.epigenome.chipatlas import compare_celltypes

# Compare H3K27ac across cell types
track, fig, axes = compare_celltypes(
    antigen="H3K27ac",
    cell_types=["K562", "HepG2", "GM12878"],
    chrom="chr8",
    start=127735000,
    end=127745000,
)
```

### Add Gene Annotations

```python
from gcell.rna.gencode import Gencode
from gcell.dna.track import Track
from gcell.epigenome import ChipAtlas
import pyBigWig
import numpy as np

# Load gene annotations
gencode = Gencode(assembly="hg38")

# Get gene coordinates - note: tss returns a list of TSS objects
gene = gencode.get_gene("MYC")
tss_position = gene.tss[0].start  # Get TSS position from first TSS object

# Define region around gene promoter
chrom = gene.chrom
start = tss_position - 5000
end = tss_position + 5000

# Load BigWig data
ca = ChipAtlas()
exp_id = "SRX190161"
bw_path = ca.download_bigwig(exp_id)
bw = pyBigWig.open(str(bw_path))
values = bw.values(chrom, start, end)
values = np.nan_to_num(np.array(values), nan=0)
bw.close()

# Create Track object
track = Track(
    chrom=chrom,
    start=start,
    end=end,
    assembly="hg38",
    tracks={"H3K27ac": values},
    conv_size=50
)

# Plot with gene body annotations
track.plot_tracks_with_genebody(
    gene_annot=gencode,
    genes_to_highlight=["MYC"],
    out_file="tracks_with_genes.png"
)
```

### Combine with Peak Data

```python
from gcell.epigenome.chipatlas import peaks_to_bed_track

# Get binary peak track
peak_track = peaks_to_bed_track(
    experiments,
    chrom="chr8",
    start=127735000,
    end=127745000,
    threshold=10,
)

# Add to visualization as BED regions
import numpy as np
peaks = np.where(peak_track > 0)[0]
# Use with track.plot_tracks(beds=[peak_regions])
```

### Stream vs Cache

```python
# Option 1: Stream directly (no disk usage, slower for repeated access)
track = experiments_to_track(
    experiments,
    chrom="chr8",
    start=127735000,
    end=127745000,
    use_cache=False,  # Stream from URL
)

# Option 2: Download and cache (faster for repeated access)
track = experiments_to_track(
    experiments,
    chrom="chr8",
    start=127735000,
    end=127745000,
    use_cache=True,  # Download first, then read locally
)
```
