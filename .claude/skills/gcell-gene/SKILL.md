---
name: gcell-gene
description: |
  Gene annotations and TSS analysis using gcell. Use this skill when users ask about:
  - GENCODE gene annotations
  - Transcription start sites (TSS)
  - Gene coordinates and metadata
  - Transcript information
  - Querying genes by genomic region
  Triggers: gene annotation, GENCODE, TSS, transcription start site, gene coordinates, transcript, GTF
---

# Gene Annotations

## Loading GENCODE Annotations

```python
from gcell.rna.gencode import Gencode

# Load annotations for specific genome
gc = Gencode('hg38')  # Human GRCh38
gc = Gencode('hg19')  # Human GRCh37
gc = Gencode('mm10')  # Mouse mm10
```

## Accessing Gene Information

```python
# Get gene by symbol
gene = gc.genes['TP53']
gene = gc.genes['BRCA1']
gene = gc.genes['MYC']

# Gene attributes
print(gene.gene_id)      # Ensembl ID
print(gene.gene_name)    # Symbol
print(gene.chrom)        # Chromosome
print(gene.start)        # Start coordinate
print(gene.end)          # End coordinate
print(gene.strand)       # + or -
print(gene.gene_type)    # protein_coding, lncRNA, etc.
```

## Transcription Start Sites (TSS)

```python
# Get primary TSS
tss = gene.tss

# Get all TSS for a gene (from all transcripts)
tss_list = gene.get_all_tss()

# TSS object attributes
for tss in tss_list:
    print(tss.chrom, tss.position, tss.strand)
```

## Query Genes by Region

```python
# Find all genes in a genomic region
genes_in_region = gc.query('chr17', 41196312, 41277500)

# Returns list of Gene objects
for gene in genes_in_region:
    print(gene.gene_name, gene.start, gene.end)
```

## Transcript Information

```python
# Get transcripts for a gene
transcripts = gene.transcripts

for tx in transcripts:
    print(tx.transcript_id)
    print(tx.start, tx.end)
    print(tx.is_canonical)
```

## Key Classes

| Class | Purpose |
|-------|---------|
| `Gencode` | GENCODE annotation database |
| `Gene` | Gene with coordinates and metadata |
| `TSS` | Transcription start site |
| `Transcript` | Transcript information |

## Data Location

- Annotations: `~/.gcell_data/annotations/`
- Override: `GCELL_ANNOTATION_DIR` environment variable
