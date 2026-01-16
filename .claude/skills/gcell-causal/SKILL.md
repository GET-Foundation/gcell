---
name: gcell-causal
description: |
  Causal network analysis using gcell. Use this skill when users ask about:
  - Inferring causal relationships from expression data
  - LiNGAM causal discovery algorithm
  - Regulatory network construction
  - Visualizing causal/regulatory networks
  Triggers: causal network, LiNGAM, causal inference, regulatory network, causal discovery, gene regulatory network
---

# Causal Network Analysis

## LiNGAM Causal Discovery

LiNGAM (Linear Non-Gaussian Acyclic Model) infers causal structure from observational data by exploiting non-Gaussianity.

```python
from gcell.utils.lingam import run_lingam
import pandas as pd

# Prepare expression data
# Rows = samples, Columns = genes
expression_df = pd.DataFrame(...)

# Run LiNGAM to infer causal structure
causal_matrix = run_lingam(expression_df)

# causal_matrix[i,j] represents causal effect of gene j on gene i
# Non-zero values indicate direct causal relationships
```

## Interpreting Results

```python
import numpy as np

# Get significant causal edges
threshold = 0.1
edges = np.where(np.abs(causal_matrix) > threshold)

for i, j in zip(edges[0], edges[1]):
    gene_from = expression_df.columns[j]
    gene_to = expression_df.columns[i]
    effect = causal_matrix[i, j]
    print(f"{gene_from} -> {gene_to}: {effect:.3f}")
```

## Network Visualization

```python
from gcell.utils.causal_lib import visualize_network

# Visualize the regulatory network
visualize_network(causal_matrix, top_edges=50)

# Parameters:
# - causal_matrix: The adjacency matrix from LiNGAM
# - top_edges: Number of strongest edges to display
```

## Building Gene Regulatory Networks

```python
# Typical workflow for GRN construction

# 1. Filter to genes of interest (e.g., TFs and targets)
tf_genes = ['STAT3', 'MYC', 'TP53', 'GATA1']
target_genes = ['BCL2', 'CCND1', 'BAX', 'MDM2']
all_genes = tf_genes + target_genes

# 2. Subset expression data
expr_subset = expression_df[all_genes]

# 3. Run causal discovery
grn_matrix = run_lingam(expr_subset)

# 4. Visualize
visualize_network(grn_matrix, top_edges=30)
```

## Combining with Cell Type Analysis

```python
from gcell.cell.celltype import GETDemoLoader

# Load cell type data
loader = GETDemoLoader()
ct = loader.load_celltype('Monocyte')

# Get gene-by-motif regulatory matrix
gbm = ct.get_gene_by_motif()

# Use as prior for causal analysis
# or compare with data-driven causal network
```

## Key Functions

| Function | Purpose |
|----------|---------|
| `run_lingam()` | Infer causal structure from data |
| `visualize_network()` | Plot causal/regulatory network |

## Tips

- LiNGAM assumes: linear relationships, non-Gaussian noise, acyclic graph
- More samples improve causal discovery accuracy
- Consider subsetting to relevant genes to reduce computational cost
- Validate inferred edges with known biology
