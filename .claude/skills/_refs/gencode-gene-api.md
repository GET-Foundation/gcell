# Gene and Gencode API Reference

This is a shared reference for the `Gencode`, `Gene`, and `TSS` classes used across gcell skills.

## Gencode Class

```python
from gcell.rna.gencode import Gencode

# Initialize with genome assembly
gencode = Gencode(assembly="hg38")  # Human GRCh38
gencode = Gencode(assembly="hg19")  # Human GRCh37
gencode = Gencode(assembly="mm10")  # Mouse GRCm38
gencode = Gencode(assembly="mm39")  # Mouse GRCm39
```

### Getting Gene Objects

```python
# Get Gene object by symbol (returns Gene)
gene = gencode.get_gene("TP53")
gene = gencode.get_gene("MYC")
gene = gencode.get_gene("BRCA1")

# Get Gene object by Ensembl ID
gene = gencode.get_gene_id("ENSG00000141510")

# Get multiple genes (returns GeneSets)
genes = gencode.get_genes(["TP53", "MYC", "BRCA1"])
```

### Gencode Lookup Properties

```python
# Dictionary mappings (gene_name -> value)
gencode.gene_to_strand       # {"TP53": "-", "MYC": "+", ...}
gencode.gene_to_strand_numeric  # {"TP53": 1, "MYC": 0, ...}
gencode.gene_to_chrom        # {"TP53": "chr17", ...}
gencode.gene_to_tss          # {"TP53": 7687538, ...}
gencode.gene_to_tes          # {"TP53": 7668421, ...}
gencode.gene_to_id           # {"TP53": "ENSG00000141510", ...}
gencode.gene_to_type         # {"TP53": "protein_coding", ...}
gencode.id_to_gene           # {"ENSG00000141510": "TP53", ...}
```

## Gene Class Attributes

```python
gene = gencode.get_gene("TP53")

# Basic attributes
gene.name       # str: Gene symbol, e.g., "TP53"
gene.id         # str: Ensembl ID, e.g., "ENSG00000141510"
gene.chrom      # str: Chromosome, e.g., "chr17"
gene.strand     # str: "+" or "-"

# TSS (Transcription Start Site)
gene.tss            # list[TSS]: List of TSS objects for all transcripts
gene.tss_coordinate # int: Primary TSS coordinate
gene.tss_list       # DataFrame: TSS data with Chromosome, Start, End, Strand, gene_name, gene_id

# TES (Transcription End Site)
gene.tes            # int: Primary TES coordinate
gene.tes_list       # DataFrame: TES data with same structure as tss_list

# Full gene body coordinates
gene.genomic_range  # tuple[str, int, int, str]: (chrom, start, end, strand)
```

## TSS Class Attributes

```python
for tss in gene.tss:
    tss.name      # str: Gene symbol
    tss.peak_id   # int: TSS ID
    tss.chrom     # str: Chromosome
    tss.start     # int: TSS coordinate
    tss.strand    # str: "+" or "-"
```

## Common Patterns

### Get coordinates around gene TSS

```python
gene = gencode.get_gene("TP53")

# Get primary TSS coordinate
tss_coord = gene.tss_coordinate  # or gene.tss[0].start for first transcript

# Build promoter region (strand-aware)
if gene.strand == "+":
    # Upstream is lower coordinates
    start = tss_coord - 2000  # upstream
    end = tss_coord + 500     # downstream
else:
    # Upstream is higher coordinates (gene is on minus strand)
    start = tss_coord - 500   # downstream
    end = tss_coord + 2000    # upstream

# Use chrom for queries
chrom = gene.chrom
```

### Get full gene body

```python
gene = gencode.get_gene("TP53")
chrom, start, end, strand = gene.genomic_range
```

### Access individual TSS coordinates

```python
gene = gencode.get_gene("MYC")

# Get all TSS coordinates
for tss in gene.tss:
    print(f"{tss.chrom}:{tss.start} ({tss.strand})")

# Or from the DataFrame
print(gene.tss_list[['Chromosome', 'Start', 'Strand']])
```
