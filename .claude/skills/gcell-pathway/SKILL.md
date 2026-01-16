---
name: gcell-pathway
description: |
  Pathway enrichment analysis using gcell. Use this skill when users ask about:
  - Gene set enrichment analysis
  - GO (Gene Ontology) enrichment
  - KEGG pathway analysis
  - Reactome pathway enrichment
  - Custom pathway/gene set analysis
  Triggers: pathway enrichment, GO enrichment, KEGG, Reactome, gene set analysis, functional enrichment, ontology
---

# Pathway Enrichment Analysis

## Quick Enrichment with gprofiler

```python
from gcell.ontology.pathway import gprofiler_enrichment

# Basic enrichment analysis
gene_list = ['TP53', 'BRCA1', 'MYC', 'EGFR', 'KRAS']
results = gprofiler_enrichment(gene_list, organism='hsapiens')

# Specify data sources
results = gprofiler_enrichment(
    gene_list,
    organism='hsapiens',
    sources=['GO:BP', 'GO:MF', 'GO:CC', 'KEGG', 'REAC']
)

# Sources available:
# - GO:BP (Biological Process)
# - GO:MF (Molecular Function)
# - GO:CC (Cellular Component)
# - KEGG (KEGG pathways)
# - REAC (Reactome)
# - WP (WikiPathways)
# - TF (Transcription factors)
# - MIRNA (microRNA targets)
# - HPA (Human Protein Atlas)
# - CORUM (Protein complexes)
# - HP (Human Phenotype Ontology)
```

## Working with Results

```python
# Results is a pandas DataFrame
print(results.columns)
# ['source', 'term_id', 'term_name', 'p_value', 'significant',
#  'term_size', 'query_size', 'intersection_size', 'intersections']

# Filter significant results
significant = results[results['p_value'] < 0.05]

# Sort by p-value
top_terms = results.sort_values('p_value').head(20)

# Get genes in each term
for _, row in top_terms.iterrows():
    print(f"{row['term_name']}: {row['intersections']}")
```

## Mouse and Other Organisms

```python
# Mouse
results = gprofiler_enrichment(gene_list, organism='mmusculus')

# Rat
results = gprofiler_enrichment(gene_list, organism='rnorvegicus')

# Other organisms: use Ensembl species codes
```

## Custom Pathways from GMT Files

```python
from gcell.ontology.pathway import Pathways

# Load custom gene sets from GMT file
pathways = Pathways.from_gmt('custom_pathways.gmt')

# Run enrichment against custom pathways
background_genes = [...]  # All expressed genes
enriched = pathways.enrichment(gene_list, background_genes)
```

## Key Functions and Classes

| Name | Purpose |
|------|---------|
| `gprofiler_enrichment()` | Quick enrichment via g:Profiler |
| `Pathways` | Custom pathway collections |
| `Pathways.from_gmt()` | Load GMT format gene sets |
| `Pathways.enrichment()` | Run enrichment analysis |

## Tips

- Always use appropriate background genes when possible
- Multiple testing correction is applied automatically
- Use specific sources (e.g., just 'GO:BP') to reduce multiple testing burden
- Gene symbols should match the organism (human: HUGO symbols)
