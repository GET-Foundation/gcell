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
refs:
  - _refs/gencode-gene-api.md
---

# Gene Annotations

## Loading GENCODE Annotations

```python
from gcell.rna.gencode import Gencode

# Load annotations for specific genome
gencode = Gencode(assembly="hg38")  # Human GRCh38
gencode = Gencode(assembly="hg19")  # Human GRCh37
gencode = Gencode(assembly="mm10")  # Mouse mm10
```

## Accessing Gene Information

```python
# Get gene by symbol
gene = gencode.get_gene("TP53")
gene = gencode.get_gene("BRCA1")
gene = gencode.get_gene("MYC")

# Gene attributes
print(gene.id)         # Ensembl ID: "ENSG00000141510"
print(gene.name)       # Symbol: "TP53"
print(gene.chrom)      # Chromosome: "chr17"
print(gene.strand)     # Strand: "-"

# Gene coordinates
print(gene.tss_coordinate)  # Primary TSS coordinate
print(gene.tes)             # Primary TES coordinate

# Full gene body
chrom, start, end, strand = gene.genomic_range
```

## Transcription Start Sites (TSS)

```python
# Get list of TSS objects (one per transcript)
for tss in gene.tss:
    print(tss.chrom, tss.start, tss.strand)

# Get primary TSS coordinate
tss_coord = gene.tss_coordinate

# Access TSS DataFrame for detailed info
print(gene.tss_list)  # DataFrame with Chromosome, Start, End, Strand, gene_name, gene_id
```

## Query Genes by Region

```python
# Find genes in a genomic region
result = gencode.query_region("chr17", 41196312, 41277500)

# Returns DataFrame with matching genes
print(result[['gene_name', 'Chromosome', 'Start', 'End', 'Strand']])
```

## Gencode Lookup Properties

```python
# Quick lookups without creating Gene objects
strand = gencode.gene_to_strand["TP53"]  # "-"
chrom = gencode.gene_to_chrom["TP53"]    # "chr17"
tss = gencode.gene_to_tss["TP53"]        # 7687538
tes = gencode.gene_to_tes["TP53"]        # 7668421
gene_type = gencode.gene_to_type["TP53"] # "protein_coding"
gene_id = gencode.gene_to_id["TP53"]     # "ENSG00000141510"
```

## Key Classes

| Class | Purpose |
|-------|---------|
| `Gencode` | GENCODE annotation database |
| `Gene` | Gene with coordinates and TSS/TES |
| `TSS` | Transcription start site object |
| `GeneSets` | Collection of Gene objects |

## Data Location

- Annotations: `~/.gcell_data/annotations/`
- Override: `GCELL_ANNOTATION_DIR` environment variable
