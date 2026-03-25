"""CD70 Endothelial Cell Jacobian and Motif Regulatory Analysis Tutorial.

This tutorial demonstrates how to use the GET model's Jacobian-based analysis
to understand the transcription factor (TF) regulatory logic controlling CD70
expression in endothelial cells.

CD70 (TNFSF7) is a TNF superfamily ligand for CD27. Its aberrant expression
in endothelial cells has been linked to tumour angiogenesis, immunomodulation,
and autoimmune vascular disease. Understanding which TFs control CD70 in
endothelial cells can reveal potential therapeutic targets.

The GET model computes Jacobians — partial derivatives of predicted gene
expression with respect to chromatin accessibility and motif content of each
regulatory region. Summing absolute mean Jacobians across regions gives a
per-TF-motif importance score that quantifies which TFs most strongly drive
CD70 transcription in a particular cell type.

Sections
--------
1. Load endothelial cell type with Jacobian data.
2. CD70 Jacobian motif analysis: top 20 TF motifs driving CD70.
3. CD70 regulatory regions: which genomic elements matter most.
4. Gene-by-motif matrix: compare CD70's TF profile with housekeeping and
   endothelial marker genes.
5. Motif subnet visualisation: interactive causal regulatory networks for
   the top 2 TF drivers of CD70.

Usage
-----
Run from the repository root::

    uv run python tutorials/cd70_endothelial_regulation.py

Output files
------------
tutorials/figures/cd70_endo_top_motifs.png
tutorials/figures/cd70_endo_regulatory_regions.png
tutorials/figures/cd70_motif_comparison_heatmap.png
tutorials/figures/cd70_motif_subnet_1.html
tutorials/figures/cd70_motif_subnet_2.html
"""

import logging
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")  # non-interactive backend — must come before pyplot import

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
log = logging.getLogger(__name__)

# Output directory for all generated figures
FIGURES_DIR = Path(__file__).parent / "figures"
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

# ===========================================================================
# Section 1: Load endothelial cell type with Jacobian data
# ===========================================================================
# GETDemoLoader provides access to pre-inferred GET model outputs for a set
# of cell types stored on S3. Loading with jacob=True enables Jacobian-based
# regulatory attribution methods (get_gene_jacobian_summary, plotly_motif_subnet,
# etc.). Jacobians are computed once by the GET model and cached in zarr format.

log.info("Section 1: Loading endothelial cell type with Jacobians")

from gcell.cell.celltype import GETDemoLoader
from gcell.dna.nr_motif_v1 import NrMotifV1

loader = GETDemoLoader()

print("Available cell types:")
for ct_name in loader.available_celltypes:
    print(f"  {ct_name}")

# Identify endothelial cell types (case-insensitive search)
endo_names = [
    c for c in loader.available_celltypes if "Endothelial" in c or "endothelial" in c
]
if not endo_names:
    raise ValueError(
        "No endothelial cell type found. Available types: "
        + str(loader.available_celltypes)
    )

endo_name = endo_names[0]
log.info("Loading: %s (jacob=True)", endo_name)

# jacob=True loads the GET Jacobian tensor from zarr and makes it available
# via self.jacobs. Without this flag the methods used below will raise errors.
ct_endo = loader.load_celltype(endo_name, jacob=True)
log.info("Loaded successfully: %s", ct_endo.celltype_name)

# ===========================================================================
# Section 2: CD70 Jacobian analysis — which TF motifs regulate CD70?
# ===========================================================================
# get_gene_jacobian_summary(gene, axis='motif') returns a pd.Series indexed
# by NrMotifV1 cluster names (+ "Accessibility"). Each value is the sum of
# |absmean| Jacobian scores across all TSSs of that gene, reflecting how
# strongly each TF's binding motif predicts CD70 expression.
#
# Biological context: high-scoring activating motifs (positive Jacobian) are
# TFs that, when their binding sites are accessible, correlate with increased
# CD70 expression. These are the regulatory drivers we want to discover.

log.info("Section 2: CD70 Jacobian motif analysis")

motif_scores = ct_endo.get_gene_jacobian_summary("CD70", axis="motif")
log.info("Motif attribution scores — shape: %s", motif_scores.shape)

# Drop the Accessibility pseudo-feature to focus on TF motifs
motif_scores_tf = motif_scores.drop("Accessibility", errors="ignore")
top20 = motif_scores_tf.sort_values(ascending=False).head(20)

print(f"\nTop 20 TF motifs driving CD70 expression in {endo_name}:")
print(top20.to_string())

# --- Figure: Horizontal bar chart of top 20 TF motifs ---
fig, ax = plt.subplots(figsize=(10, 8))
colors = ["#d62728" if v >= 0 else "#1f77b4" for v in top20.values]
ax.barh(range(len(top20)), top20.values, color=colors, edgecolor="white", height=0.7)
ax.set_yticks(range(len(top20)))
ax.set_yticklabels(top20.index, fontsize=9)
ax.invert_yaxis()  # highest score at the top
ax.axvline(0, color="black", linewidth=0.8, linestyle="--")
ax.set_xlabel("Jacobian attribution score (Σ |absmean| across TSSs)", fontsize=11)
ax.set_title(
    f"Top 20 TF Motifs Regulating CD70 in {endo_name}\n"
    "(red = activating influence, blue = repressive)",
    fontsize=12,
)
# Add value labels on each bar
for idx, (val, label) in enumerate(zip(top20.values, top20.index)):
    ax.text(
        val + 0.001 * (1 if val >= 0 else -1),
        idx,
        f"{val:.3f}",
        va="center",
        ha="left" if val >= 0 else "right",
        fontsize=7.5,
        color="gray",
    )
sns.despine(ax=ax)
plt.tight_layout()
out_motifs = FIGURES_DIR / "cd70_endo_top_motifs.png"
fig.savefig(out_motifs, dpi=150, bbox_inches="tight")
plt.close(fig)
log.info("Saved: %s", out_motifs)

# ===========================================================================
# Section 3: CD70 regulatory regions
# ===========================================================================
# get_gene_jacobian_summary(gene, axis='region') returns a pd.DataFrame with
# columns ['index', 'Chromosome', 'Start', 'End', 'Score']. Each row is a
# genomic region; Score is the sum of absolute-mean Jacobian norms across all
# TSSs of CD70, quantifying each region's contribution to predicted expression.
#
# Biological context: high-score regions are the putative cis-regulatory
# elements (enhancers, promoters) that the model predicts are functionally
# important for CD70 regulation in this cell type.

log.info("Section 3: CD70 cis-regulatory regions")

region_scores = ct_endo.get_gene_jacobian_summary("CD70", axis="region")
log.info(
    "Region scores: %d regions, chromosomes: %s",
    len(region_scores),
    region_scores.Chromosome.unique().tolist(),
)

# Attempt to use the built-in plot_gene_regions method.
# It filters regions to the gene's chromosome (using self.focus), then draws
# each region as a rectangle whose height is proportional to normalised score.
try:
    result = ct_endo.plot_gene_regions("CD70", plotly=False)
    # The method returns (fig, ax)
    fig_r, ax_r = result
    ax_r.set_title(
        f"CD70 Cis-Regulatory Elements — {endo_name}\n"
        "(rectangle height = normalised Jacobian importance score)",
        fontsize=11,
    )
    out_regions = FIGURES_DIR / "cd70_endo_regulatory_regions.png"
    fig_r.savefig(out_regions, dpi=150, bbox_inches="tight")
    plt.close(fig_r)
    log.info("Saved: %s", out_regions)
except Exception as exc:
    # Fallback: manually plot from the region_scores DataFrame
    log.warning("plot_gene_regions raised %s — using manual fallback", exc)
    df_r = region_scores.copy()
    df_r = df_r[df_r.Score > 0].sort_values("Score", ascending=False)
    if df_r.empty:
        df_r = region_scores.copy()
    df_r["NormHeight"] = df_r.Score.abs() / df_r.Score.abs().max()

    fig_r, ax_r = plt.subplots(figsize=(14, 3))
    cmap = plt.cm.get_cmap("Reds")
    for _, row in df_r.iterrows():
        fc = cmap(row.NormHeight * 0.85 + 0.15)
        ax_r.add_patch(
            plt.Rectangle(
                (row.Start, 0),
                row.End - row.Start,
                row.NormHeight,
                facecolor=fc,
                edgecolor="none",
                alpha=0.85,
            )
        )
    xmin = region_scores.Start.min()
    xmax = region_scores.End.max()
    ax_r.set_xlim(xmin, xmax)
    ax_r.plot([xmin, xmax], [0, 0], color="black", linewidth=0.8)
    ax_r.set_ylim(-0.05, 1.15)
    ax_r.set_yticks([])
    chrom = region_scores.Chromosome.iloc[0]
    ax_r.set_xlabel(f"Genomic position ({chrom})", fontsize=10)
    ax_r.set_title(
        f"CD70 Cis-Regulatory Elements — {endo_name}\n"
        "(rectangle height = normalised Jacobian importance score)",
        fontsize=11,
    )
    # Annotate top 3 regions by score
    top3 = df_r.nlargest(3, "NormHeight")
    for _, row in top3.iterrows():
        mid = (row.Start + row.End) / 2
        ax_r.annotate(
            f"Score={row.Score:.2f}",
            xy=(mid, row.NormHeight),
            xytext=(mid, row.NormHeight + 0.1),
            fontsize=7,
            ha="center",
            arrowprops=dict(arrowstyle="-", color="gray", lw=0.5),
            color="gray",
        )
    sns.despine(ax=ax_r, left=True)
    plt.tight_layout()
    out_regions = FIGURES_DIR / "cd70_endo_regulatory_regions.png"
    fig_r.savefig(out_regions, dpi=150, bbox_inches="tight")
    plt.close(fig_r)
    log.info("Saved: %s", out_regions)

# ===========================================================================
# Section 4: Gene-by-motif matrix — CD70 vs housekeeping & endothelial markers
# ===========================================================================
# get_gene_by_motif() returns a GeneByMotif object. Its .data attribute is a
# pd.DataFrame indexed by gene_name. Each row contains mean Jacobian × input
# (attribution) scores per motif feature for one TSS. Genes with multiple
# TSSs appear multiple times; we average across TSSs.
#
# Biological context: comparing CD70's motif attribution profile with
# housekeeping genes (GAPDH, ACTB — should be broadly accessible, driven by
# universal TFs) and endothelial markers (VWF, CDH5, PECAM1 — cell-type
# identity genes) reveals whether CD70 in endothelial cells is controlled by
# the same TFs as cell-identity genes or by distinct regulatory logic.

log.info("Section 4: Gene-by-motif comparison heatmap")

gbm = ct_endo.get_gene_by_motif()
log.info("gene_by_motif data shape: %s", gbm.data.shape)

GENES_OF_INTEREST = ["CD70", "GAPDH", "ACTB", "VWF", "CDH5", "PECAM1"]
ENDO_MARKERS = {"VWF", "CDH5", "PECAM1"}

gene_profiles: dict = {}
for gene in GENES_OF_INTEREST:
    if gene in gbm.data.index:
        profile = gbm.data.loc[gene]
        if isinstance(profile, pd.DataFrame):
            # Multiple TSSs for this gene — take the mean across TSSs
            profile = profile.mean(axis=0)
        gene_profiles[gene] = profile
        log.info("  Found %s in gene_by_motif", gene)
    else:
        log.warning("  Gene %s not found in gene_by_motif — skipping", gene)

if len(gene_profiles) < 2:
    log.warning("Fewer than 2 genes found; skipping Section 4 heatmap")
else:
    # Reference gene for selecting top 15 motifs
    ref_gene = "CD70" if "CD70" in gene_profiles else list(gene_profiles.keys())[0]
    top15_motifs = (
        gene_profiles[ref_gene]
        .drop("Accessibility", errors="ignore")
        .abs()
        .sort_values(ascending=False)
        .head(15)
        .index
    )

    # Build matrix: rows = genes, columns = top 15 motifs
    heatmap_df = pd.DataFrame(gene_profiles).T[top15_motifs]

    # Z-score normalise each gene's row so that visual contrast reflects
    # relative motif specificity rather than absolute score magnitude.
    row_mean = heatmap_df.mean(axis=1)
    row_std = heatmap_df.std(axis=1).replace(0, 1)
    heatmap_norm = heatmap_df.subtract(row_mean, axis=0).divide(row_std, axis=0)

    # --- Figure: heatmap ---
    n_genes = len(heatmap_norm)
    fig_h, ax_h = plt.subplots(figsize=(15, max(4, n_genes * 0.85 + 2.5)))
    vmax = float(np.nanmax(np.abs(heatmap_norm.values)))
    im = ax_h.imshow(
        heatmap_norm.values,
        aspect="auto",
        cmap="RdBu_r",
        vmin=-vmax,
        vmax=vmax,
        interpolation="nearest",
    )
    ax_h.set_xticks(range(len(top15_motifs)))
    ax_h.set_xticklabels(top15_motifs, rotation=45, ha="right", fontsize=8.5)

    # Mark endothelial marker genes with a star
    ylabels = [
        f"★  {g}" if g in ENDO_MARKERS else ("→  " + g if g == ref_gene else f"    {g}")
        for g in heatmap_norm.index
    ]
    ax_h.set_yticks(range(n_genes))
    ax_h.set_yticklabels(ylabels, fontsize=10)

    cbar = plt.colorbar(im, ax=ax_h, fraction=0.03, pad=0.02)
    cbar.set_label("Z-scored motif attribution (row-normalised)", fontsize=9)
    ax_h.set_title(
        f"TF Motif Influence on Gene Expression — {endo_name}\n"
        f"★ = endothelial markers  |  → = reference gene  |  "
        f"top-15 motifs by {ref_gene} attribution",
        fontsize=11,
    )
    plt.tight_layout()
    out_heatmap = FIGURES_DIR / "cd70_motif_comparison_heatmap.png"
    fig_h.savefig(out_heatmap, dpi=150, bbox_inches="tight")
    plt.close(fig_h)
    log.info("Saved: %s", out_heatmap)

    print(f"\nTop-15 motifs (by {ref_gene} attribution):")
    for rank, (motif_name, val) in enumerate(
        gene_profiles[ref_gene]
        .drop("Accessibility", errors="ignore")
        .abs()
        .sort_values(ascending=False)
        .head(15)
        .items(),
        start=1,
    ):
        print(f"  {rank:2d}. {motif_name:<45s}  |score|={val:.4f}")

# ===========================================================================
# Section 5: Motif subnet visualisation for the top 2 CD70 drivers
# ===========================================================================
# plotly_motif_subnet(motif_db, m) generates an interactive directed network
# graph showing the local causal neighbourhood of motif `m` in the LiNGAM
# causal graph over all TF motifs (computed from the gene-by-motif matrix).
#
# This reveals:
#   - Which upstream TFs (parents) are predicted to causally activate/repress m.
#   - Which downstream TFs (children) are regulated by m.
#   - The broader co-regulatory context of the top CD70 drivers.
#
# Biological context: a TF appearing as a strong CD70 driver AND as a hub with
# many downstream targets in the causal graph is a good candidate master
# regulator worth experimental validation.
#
# NOTE: The NrMotifV1 motif database must be loaded explicitly here because
# GETCellType (returned by GETDemoLoader) does not auto-load it. NrMotifV1
# is used to annotate network nodes with the TF genes associated with each
# motif cluster and their mean predicted expression levels.

log.info("Section 5: Motif subnet visualisation")

motif_db = NrMotifV1.load_from_pickle()
log.info("NrMotifV1 loaded: %d clusters", len(motif_db.cluster_names))

# Use the top 2 TF motifs from the CD70 Jacobian summary (already computed
# in Section 2, Accessibility excluded). These are the two TFs most predicted
# to drive CD70 expression in this endothelial cell type.
top2_motifs = motif_scores_tf.sort_values(ascending=False).head(2).index.tolist()
log.info("Top 2 motifs selected for subnet: %s", top2_motifs)

for i, m in enumerate(top2_motifs, start=1):
    out_html = FIGURES_DIR / f"cd70_motif_subnet_{i}.html"
    log.info("  Generating subnet %d/%d: %s", i, len(top2_motifs), m)
    try:
        # plotly_motif_subnet internally:
        #   1. Obtains the LiNGAM causal graph via self.gene_by_motif.get_causal()
        #      (loaded from zarr cache if available — fast; computed fresh if not).
        #   2. Extracts the local 'neighbors' subgraph of motif m.
        #   3. Annotates each node with TF gene names and mean predicted expression.
        # The returned plotly Figure is fully interactive HTML — hover for details.
        fig_net = ct_endo.plotly_motif_subnet(motif_db, m, type="neighbors")
        fig_net.update_layout(
            title=dict(
                text=(
                    f"Causal Motif Regulatory Subnet: {m}"
                    f"<br><sup>CD70 rank-{i} driver in {endo_name}"
                    f" (neighbours of {m})</sup>"
                ),
                font=dict(size=13),
            ),
            margin=dict(l=20, r=20, t=80, b=20),
        )
        fig_net.write_html(str(out_html), include_plotlyjs="cdn")
        log.info("  Saved: %s", out_html)
    except Exception as exc:
        log.warning("  Could not generate subnet for %s: %s", m, exc)

# ===========================================================================
# Summary
# ===========================================================================
print("\n" + "=" * 60)
print("CD70 Endothelial Regulatory Analysis — Complete")
print("=" * 60)
print(f"Cell type  : {endo_name}")
print(f"Output dir : {FIGURES_DIR}")
print("\nFiles produced:")
produced = sorted(FIGURES_DIR.glob("cd70_endo*")) + sorted(
    FIGURES_DIR.glob("cd70_motif*")
)
for f in produced:
    try:
        rel = f.relative_to(Path(__file__).parent.parent)
    except ValueError:
        rel = f
    print(f"  {rel}")
print("=" * 60)
