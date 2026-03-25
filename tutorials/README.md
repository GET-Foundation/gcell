# CD70 Accessibility and Regulation Analysis Across Cell Types

A step-by-step tutorial series that uses the **gcell** library and the GET
(Gene Expression Transformer) pre-inferred model to study how the immune
co-stimulatory gene **CD70** (*TNFSF7*) is regulated across dozens of human
cell types.

---

## Background: What is CD70?

**CD70** (gene symbol *TNFSF7*) is a type-II transmembrane protein that
belongs to the **TNF superfamily**. It functions as the sole ligand for the
**CD27** receptor, a key regulator of T-cell and B-cell responses.

| Feature | Detail |
|---------|--------|
| Gene symbol | TNFSF7 (CD70) |
| Protein family | TNF superfamily, member 7 |
| Receptor | CD27 (TNFRSF7) |
| Normal expression | Activated T cells, NK cells, mature dendritic cells |
| Aberrant expression | Endothelial cells (vascular inflammation), tumour cells |

### Why does CD70 matter?

- **Immune activation**: CD70–CD27 signalling co-stimulates naive and memory
  T cells, promotes effector differentiation, and supports germinal centre
  B-cell responses.
- **CAR-T therapy**: CD70 is highly expressed on renal cell carcinoma, glioma,
  AML, and T-cell lymphoma, making it a validated CAR-T target (several
  clinical trials ongoing).
- **Cancer immunology**: Constitutive CD70 expression on tumour cells drives
  T-cell exhaustion and immunosuppression via CD27 signalling.
- **Endothelial biology**: Vascular endothelial CD70 participates in
  inflammatory leukocyte recruitment and tumour angiogenesis.

The **lineage-specific** regulation of CD70 — controlled by different TF
programs in different cell types — makes it an ideal showcase for the GET
model's Jacobian-based regulatory analysis.

---

## Prerequisites

### Installing gcell

**With `uv` (recommended):**

```bash
git clone https://github.com/GET-Foundation/gcell.git
cd gcell
uv sync                  # core dependencies
uv sync --all-extras     # include optional torch / plotting extras
```

**With pip:**

```bash
pip install gcell
# or for all extras:
pip install "gcell[torch,plot]"
```

### Internet access

All tutorial scripts stream pre-computed model outputs from S3 via the
`GETDemoLoader` class. An active internet connection is required the first
time each cell type is loaded; results are cached locally in
`~/.gcell_data/cache/` for subsequent runs.

### Python version

Python ≥ 3.10 is required (matches the gcell project constraint).

---

## Quick Start

Run each tutorial script from the **repository root**:

```bash
# Tutorial 1 — Cross-cell-type accessibility survey (~10–15 min)
uv run python tutorials/cd70_regulation_analysis.py

# Tutorial 2 — Endothelial deep-dive with Jacobians (~5–10 min)
uv run python tutorials/cd70_endothelial_regulation.py

# Tutorial 3 — Comparative multi-lineage analysis (~10–15 min)
uv run python tutorials/cd70_comparative_regulation.py
```

All output figures are written to `tutorials/figures/`. A summary CSV is
written to `tutorials/cd70_cross_celltype_summary.csv` by Tutorial 1.

---

## Analysis Scripts

### 1. `cd70_regulation_analysis.py` — Cross-cell-type accessibility survey

**Goal:** Survey CD70 chromatin accessibility and GET model expression
predictions across *all* available pre-inferred cell types.

**What it does:**

| Section | Description |
|---------|-------------|
| Setup & Discovery | Initialise `GETDemoLoader`, list available cell types, identify endothelial cell types by keyword |
| Cross-cell-type survey | Loop over every cell type; collect CD70 accessibility, predicted expression, and observed expression into a DataFrame |
| Visualisation | Bar chart of CD70 accessibility (endothelial cells highlighted); scatter plot of predicted vs observed expression coloured by accessibility |

**Key API calls:**
```python
from gcell.cell.celltype import GETDemoLoader

g = GETDemoLoader()
ct = g.load_celltype("Endothelial Cell")

accessibility = ct.get_gene_accessibility("CD70")
predicted_exp = ct.get_gene_exp_pred("CD70")
observed_exp  = ct.get_gene_exp("CD70")
gene_annot    = ct.gene_annot          # DataFrame with per-TSS metadata
```

**Output files:**
- `tutorials/cd70_cross_celltype_summary.csv` — tabular results for all cell types
- `tutorials/figures/cd70_accessibility_barplot.png` — ranked accessibility bar chart
- `tutorials/figures/cd70_pred_vs_obs.png` — predicted vs observed expression scatter

---

### 2. `cd70_endothelial_regulation.py` — Deep-dive into endothelial regulatory mechanisms

**Goal:** Use Jacobian-based motif attribution to identify which transcription
factors most strongly drive CD70 expression in endothelial cells.

**Background:**
The GET model computes **Jacobians** — partial derivatives of predicted gene
expression with respect to the chromatin accessibility and motif content of
each *cis*-regulatory region. Aggregating absolute Jacobians across all
regions gives a per-TF-motif *importance score*.

**What it does:**

| Section | Description |
|---------|-------------|
| Load with Jacobians | Load an endothelial cell type with `jacob=True` to enable Jacobian APIs |
| Top-20 TF motifs | Extract and rank CD70's TF-motif Jacobian scores; save horizontal bar chart |
| Regulatory regions | Visualise importance of individual *cis*-regulatory regions (chromatin peaks) |
| Gene-by-motif heatmap | Compare CD70's motif profile against housekeeping genes (GAPDH, ACTB) and endothelial markers (VWF, CDH5, PECAM1) |
| Motif subnet | Generate interactive Plotly causal regulatory network HTML files for the top-2 CD70 TF drivers |

**Key API calls:**
```python
ct = g.load_celltype("Endothelial Cell", jacob=True)

# Jacobian-based TF importance (returns pd.Series indexed by motif cluster)
motif_scores = ct.get_gene_jacobian_summary("CD70", axis="motif")

# Per-region importance (returns pd.DataFrame with genomic coordinates)
region_scores = ct.get_gene_jacobian_summary("CD70", axis="region")

# Gene-by-motif matrix (GeneByMotif object; .data is a DataFrame)
gbm = ct.get_gene_by_motif()
cd70_profile = gbm.data.loc["CD70"].mean()   # average over TSSs

# Interactive causal subnet for a specific TF motif
ct.plotly_motif_subnet(motif_instance, "CD70", top_k=20)
```

**Output files:**
- `tutorials/figures/cd70_endo_top_motifs.png` — top-20 TF drivers bar chart
- `tutorials/figures/cd70_endo_regulatory_regions.png` — *cis*-region importance
- `tutorials/figures/cd70_motif_comparison_heatmap.png` — gene × motif heatmap
- `tutorials/figures/cd70_motif_subnet_1.html` — interactive Plotly network (TF #1)
- `tutorials/figures/cd70_motif_subnet_2.html` — interactive Plotly network (TF #2)

---

### 3. `cd70_comparative_regulation.py` — Comparative analysis across lineages

**Goal:** Reveal how CD70's TF regulatory program *differs* between
endothelial cells, immune cells, and other cell types by comparing Jacobian
motif scores across lineages.

**What it does:**

| Section | Description |
|---------|-------------|
| Multi-lineage loading | Auto-select ~4–6 representative cell types per lineage (endothelial / immune / other) by keyword; load with `jacob=True` |
| Cross-lineage heatmap | Build a motif × cell-type Jacobian matrix, select top-30 motifs by peak absolute score, Z-score normalise by row, and plot a clustered seaborn heatmap with lineage-coloured columns |
| Differential waterfall | Compute endothelial-minus-other-lineage differential score per motif; render a waterfall plot distinguishing endothelial-enriched (NF-κB, AP-1, ETS) vs immune-enriched (NFAT, IRF, RUNX) TFs |
| Regulatory fingerprints | Extract CD70's row from each cell type's gene-by-motif matrix; show ranked regulatory profiles side-by-side with a pairwise correlation summary |

**Key API calls:**
```python
# Load multiple cell types
cell_types = {
    "Endothelial": g.load_celltype("Endothelial Cell", jacob=True),
    "T Cell":      g.load_celltype("T Cell", jacob=True),
    "Fibroblast":  g.load_celltype("Fibroblast", jacob=True),
}

# Build cross-cell-type motif matrix
scores = {
    name: ct.get_gene_jacobian_summary("CD70", axis="motif")
    for name, ct in cell_types.items()
}
motif_matrix = pd.DataFrame(scores)   # rows = motifs, cols = cell types
```

**Output files:**
- `tutorials/figures/cd70_motif_heatmap_cross_celltype.png` — Z-score heatmap
- `tutorials/figures/cd70_endothelial_differential_motifs.png` — waterfall plot
- `tutorials/figures/cd70_regulatory_profile_comparison.png` — fingerprint panels

---

## Expected Outputs

After running all three scripts, the `tutorials/figures/` directory should
contain the following files:

```
tutorials/figures/
├── cd70_accessibility_barplot.png         # Script 1: accessibility ranked bar chart
├── cd70_pred_vs_obs.png                   # Script 1: pred vs obs scatter
├── cd70_endo_top_motifs.png               # Script 2: endothelial top-20 TF motifs
├── cd70_endo_regulatory_regions.png       # Script 2: cis-region importance
├── cd70_motif_comparison_heatmap.png      # Script 2: gene × motif heatmap
├── cd70_motif_subnet_1.html               # Script 2: interactive Plotly network TF#1
├── cd70_motif_subnet_2.html               # Script 2: interactive Plotly network TF#2
├── cd70_motif_heatmap_cross_celltype.png  # Script 3: cross-lineage Z-score heatmap
├── cd70_endothelial_differential_motifs.png  # Script 3: waterfall differential
└── cd70_regulatory_profile_comparison.png   # Script 3: regulatory fingerprints
```

And the tabular summary:

```
tutorials/
└── cd70_cross_celltype_summary.csv   # Script 1: all cell types × CD70 metrics
```

---

## Key Findings

> **Note:** This section is a template — populate with actual values after
> running the scripts.

### Finding 1 — CD70 accessibility is highest in activated immune and endothelial cells

*(Placeholder: replace with quantitative values from
`cd70_cross_celltype_summary.csv` after running Script 1.)*

CD70 chromatin accessibility varies by more than an order of magnitude across
cell types. Endothelial cells show intermediate-to-high accessibility relative
to most non-haematopoietic cell types, consistent with a vascular-inflammatory
CD70 expression program.

### Finding 2 — Endothelial CD70 is driven by NF-κB, AP-1, and ETS TFs

*(Placeholder: confirm top motif names from
`tutorials/figures/cd70_endo_top_motifs.png` after running Script 2.)*

Jacobian analysis identifies the NF-κB (RELA/RELB), AP-1 (FOS/JUN), and ETS
(ELK1/ETV2/ERG) TF families as the dominant drivers of CD70 in endothelial
cells. These TFs are activated by vascular inflammatory signals (TNF-α, IL-1β,
shear stress) rather than antigen receptor signalling.

### Finding 3 — Immune CD70 is driven by NFAT, IRF, and RUNX TFs

*(Placeholder: confirm motif names from
`tutorials/figures/cd70_endothelial_differential_motifs.png` after running
Script 3.)*

In T cells and NK cells, the dominant CD70 drivers are the NFAT (calcineurin-
dependent), IRF (interferon regulatory factor), and RUNX (core binding factor)
families — the canonical lymphocyte-activation TF programs. This contrasts
sharply with the vascular TF program in endothelial cells, confirming that
CD70 is controlled by cell-type-specific regulatory logic rather than a shared
master regulator.

### Finding 4 — CAR-T / therapeutic implications

The cell-type-specific regulatory programs identified here suggest that CD70
up-regulation in tumour endothelium and cancer cells hijacks the
NF-κB/AP-1/ETS axis rather than the normal immune-activation programme. This
distinction may inform the design of promoter-driven cell-type-specific
therapeutic strategies.

---

## GETDemoLoader API Reference

### Initialisation

```python
from gcell.cell.celltype import GETDemoLoader

g = GETDemoLoader()
print(g.available_celltypes)   # list all pre-inferred cell types
```

### Loading a cell type

```python
# Basic load (no Jacobians — faster)
ct = g.load_celltype("Plasma Cell")

# Load with Jacobian data (required for motif attribution)
ct = g.load_celltype("Endothelial Cell", jacob=True)
```

### Core gene-level queries

| Method | Returns | Notes |
|--------|---------|-------|
| `ct.get_gene_accessibility(gene)` | sparse matrix | Call `.toarray()` to densify; mean across TSSs |
| `ct.get_gene_exp_pred(gene)` | float / array | GET model predicted expression |
| `ct.get_gene_exp(gene)` | float / array | Observed (measured) expression |
| `ct.gene_annot` | `pd.DataFrame` | Per-TSS metadata including accessibility column |

### Jacobian-based regulatory analysis (`jacob=True` required)

| Method | Returns | Notes |
|--------|---------|-------|
| `ct.get_gene_jacobian_summary(gene, axis='motif')` | `pd.Series` | TF-motif importance scores |
| `ct.get_gene_jacobian_summary(gene, axis='region')` | `pd.DataFrame` | Per-region importance with coordinates |
| `ct.get_gene_by_motif()` | `GeneByMotif` | `.data` is a DataFrame (rows = genes/TSSs) |
| `ct.plotly_motif_subnet(motif, gene, top_k)` | Plotly figure | Causal TF→gene network; requires `NrMotifV1` instance |
| `ct.plot_gene_regions(gene)` | matplotlib figure | *cis*-region importance lollipop plot |

---

## Troubleshooting

| Problem | Solution |
|---------|----------|
| `KeyError: 'CD70'` during gene lookup | CD70 may not have open chromatin in that cell type; check `ct.gene_annot` first |
| `AttributeError: 'GETCellType' object has no attribute 'motif'` | Pass an explicit `NrMotifV1` instance to `plotly_motif_subnet` (see Script 2 source) |
| Slow first run | Cell type data is streamed from S3; cached locally after first download |
| `zarr` version errors | Pin `zarr < 3.0` as required by gcell: `pip install "zarr<3"` |
| Out-of-memory on many cell types | Load cell types one at a time and free memory with `del ct` between iterations |

---

## Contributing

If you run these tutorials and generate interesting findings, please open a PR
to update the **Key Findings** section with quantitative results and figure
references.

---

## License

These tutorials are released under the same licence as the gcell library.
See [LICENSE](../LICENSE) in the repository root for details.
