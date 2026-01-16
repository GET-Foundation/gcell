---
name: gcell-dna
description: |
  DNA sequence and motif analysis using gcell. Use this skill when users ask about:
  - Loading and manipulating DNA sequences
  - Genome assembly access (hg38, hg19, mm10)
  - Motif scanning with HOCOMOCO motifs
  - Working with genomic regions (BED, narrowPeak files)
  - Transcription factor binding site prediction
  Triggers: DNA sequence, motif scanning, HOCOMOCO, genomic regions, BED file, narrowPeak, genome, transcription factor binding
---

# DNA Sequence and Motif Analysis

## Genome Access

```python
from gcell.dna.genome import Genome

# Load genome assembly
genome = Genome('hg38')  # Human GRCh38
genome = Genome('hg19')  # Human GRCh37
genome = Genome('mm10')  # Mouse mm10

# Get sequence from coordinates
seq = genome.get_sequence('chr1', 1000000, 1001000)
seq = genome.get_sequence('chr17', 41196312, 41277500)  # BRCA1 region
```

## DNA Sequences

```python
from gcell.dna.sequence import DNASequence

# Create sequence object
seq = DNASequence('ATCGATCGATCG')

# Sequence operations
rc = seq.reverse_complement()
gc_content = seq.gc_content()
```

## Genomic Regions

```python
from gcell.dna.region import GenomicRegionCollection

# Read BED/narrowPeak files
regions = GenomicRegionCollection.read_bed('peaks.bed')
regions = GenomicRegionCollection.read_narrowpeak('peaks.narrowPeak')

# Get sequences for all regions
sequences = regions.get_sequences(genome)

# Filter regions
filtered = regions.filter(lambda r: r.length > 200)

# Merge overlapping regions
merged = regions.merge()
```

## Motif Scanning with HOCOMOCO

```python
from gcell.dna.motif import MotifCollection, Motif

# Load HOCOMOCO motifs (comprehensive TF database)
motifs = MotifCollection.from_hocomoco(version='v11')

# Load specific motif
stat3 = motifs['STAT3']
p53 = motifs['TP53']

# Scan sequence for motif hits
hits = motifs.scan(seq, pvalue=1e-4)

# Scan multiple sequences
for region_seq in sequences:
    hits = motifs.scan(region_seq, pvalue=1e-4)
    for hit in hits:
        print(hit.motif_name, hit.start, hit.end, hit.strand, hit.score)
```

## Single Motif Operations

```python
# Get PWM (position weight matrix)
pwm = motif.pwm

# Score a sequence
score = motif.score(seq)

# Get consensus sequence
consensus = motif.consensus

# Scan with single motif
hits = motif.scan(seq, pvalue=1e-4)
```

## Key Classes

| Class | Purpose |
|-------|---------|
| `Genome` | Reference genome access |
| `DNASequence` | DNA sequence operations |
| `GenomicRegionCollection` | BED/peak file handling |
| `MotifCollection` | Collection of TF motifs |
| `Motif` | Single motif PWM and scanning |

## Data Location

- Genomes: `~/.gcell_data/genomes/`
- Override: `GCELL_GENOME_DIR` environment variable
