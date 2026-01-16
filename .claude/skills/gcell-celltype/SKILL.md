---
name: gcell-celltype
description: |
  Cell type regulatory analysis using gcell. Use this skill when users ask about:
  - Loading pre-inferred cell types (GET model outputs)
  - Gene-by-motif matrices showing TF influence on genes
  - Gene Jacobian analysis for regulatory importance
  - Motif subnet visualization
  - Cell type-specific gene expression patterns
  Triggers: cell type, regulatory analysis, gene expression, Jacobian, motif subnet, GET model, TF influence
---

# Cell Type Regulatory Analysis

## Loading Cell Types

```python
from gcell.cell.celltype import GETDemoLoader

# Initialize loader
loader = GETDemoLoader()

# List available pre-inferred cell types
print(loader.available_celltypes)

# Load a specific cell type
ct = loader.load_celltype('Plasma Cell')
ct = loader.load_celltype('CD4+ T Cell')
ct = loader.load_celltype('Monocyte')
```

## Gene-by-Motif Analysis

The gene-by-motif matrix shows how transcription factor motifs influence gene expression in a cell type.

```python
# Get gene-by-motif matrix
gbm = ct.get_gene_by_motif()

# gbm is a DataFrame with genes as rows, motifs as columns
# Values represent regulatory influence scores
print(gbm.shape)
print(gbm.loc['MYC'])  # TF influences on MYC
```

## Gene Jacobian Analysis

Jacobian analysis reveals which regulatory elements most influence a gene's expression.

```python
# Get Jacobian summary for a specific gene
jacobian = ct.get_gene_jacobian_summary('MYC')
jacobian = ct.get_gene_jacobian_summary('TP53')

# Jacobian shows importance of each motif for the gene
print(jacobian.sort_values(ascending=False).head(20))
```

## Motif Subnet Visualization

Visualize the regulatory network around a specific motif/TF.

```python
# Interactive plotly visualization
ct.plotly_motif_subnet(motif_name='STAT3', top_genes=20)
ct.plotly_motif_subnet(motif_name='PU.1', top_genes=30)
ct.plotly_motif_subnet(motif_name='GATA1', top_genes=15)

# Parameters:
# - motif_name: Name of the motif/TF to center the network on
# - top_genes: Number of most influenced genes to show
```

## Key Classes

| Class | Purpose |
|-------|---------|
| `GETDemoLoader` | Load pre-inferred cell types |
| `GETCellType` | Cell type analysis container |
| `GETHydraCellType` | Multi-cell type analysis |

## Data Location

Pre-inferred cell type data is downloaded automatically to `~/.gcell_data/` on first use.
