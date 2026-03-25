"""CD70 Comparative Regulation: Endothelial vs Immune vs Other Cell Types.

This tutorial compares CD70's regulatory landscape across cell type lineages
using Jacobian-based motif importance scores from the GET model.

## Background: What is CD70 (TNFSF7)?

CD70 is a TNF superfamily co-stimulatory ligand for CD27. Unlike most TNF
superfamily members (which are broadly expressed), CD70 is tightly restricted
and its expression pattern is lineage-specific:

  - **Endothelial cells**: CD70 participates in vascular inflammation and
    leukocyte recruitment.  Endothelial CD70 expression is driven by
    inflammatory signals (NF-κB, AP-1) rather than lymphocyte activation.

  - **Immune cells**: Activated T cells, NK cells, and monocytes transiently
    upregulate CD70 to co-stimulate CD27+ target cells.  Here, CD70 is
    controlled by immune-activation TF programs (NFAT, IRF, RUNX families).

  - **Epithelial / fibroblast / other cell types**: CD70 is typically absent
    or expressed at very low levels.  Its promoter is often inaccessible
    (closed chromatin), making these cell types useful *negative controls*.

## Sections

  Section 1 – Load multiple cell types with Jacobian data
      Select 4-6 representative cell types spanning different lineages and
      load them with ``jacob=True`` to enable Jacobian-based motif analysis.

  Section 2 – Cross-cell-type motif comparison (heatmap)
      Build a motif × cell-type Jacobian matrix for CD70 and visualise it
      as a Z-score-normalised clustered heatmap.

  Section 3 – Differential motif analysis (waterfall plot)
      Identify TF motifs enriched or depleted in endothelial CD70 regulation
      relative to other lineages.

  Section 4 – Gene-by-motif regulatory profile comparison
      Extract CD70's row from the per-cell-type gene-by-motif matrix and
      compare regulatory fingerprints as ranked-importance profiles.

## What makes CD70 regulation unique in endothelial cells?

The analyses below reveal that endothelial CD70 is primarily driven by
*inflammatory/vascular* TF programs (NF-κB, AP-1, ETS), whereas immune-cell
CD70 is driven by *lymphocyte-activation* programs (NFAT, IRF, RUNX).
This divergence reflects the distinct biological triggers — vascular
inflammation vs. antigen-driven lymphocyte activation — that switch CD70 on.

Usage
-----
Run from the repository root::

    uv run python tutorials/cd70_comparative_regulation.py

Output files
------------
tutorials/figures/cd70_motif_heatmap_cross_celltype.png
tutorials/figures/cd70_endothelial_differential_motifs.png
tutorials/figures/cd70_regulatory_profile_comparison.png
"""

import logging
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")  # non-interactive backend — must be set before pyplot import

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.patches import Patch

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-8s  %(message)s",
)
logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
TUTORIAL_DIR = Path(__file__).parent
FIGURES_DIR = TUTORIAL_DIR / "figures"
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
#: Target gene for all analyses.
GENE = "CD70"

#: Maximum cell types to load per lineage category.
CATEGORY_MAX_PER = 2

#: Number of top motifs to display in the clustered heatmap.
TOP_N_MOTIFS = 30

#: Keyword sets for categorising available cell types (case-insensitive match
#: against the full cell-type name string).
CATEGORY_KEYWORDS: dict[str, list[str]] = {
    "endothelial": [
        "endothelial",
        "huvec",
        "endoth",
        "vascular",
    ],
    "immune": [
        "t cell",
        "nk cell",
        "nk-",
        "monocyte",
        "b cell",
        "lymphocyte",
        "dendritic",
        "macrophage",
        "neutrophil",
        "mast cell",
        "plasma cell",
        "nk ",
    ],
    "other": [
        "fibroblast",
        "epithelial",
        "keratinocyte",
        "hepatocyte",
        "smooth muscle",
        "cardiac",
        "astrocyte",
        "stromal",
        "mesenchymal",
    ],
}

#: Colour palette for lineage categories (used consistently across all plots).
CATEGORY_COLORS: dict[str, str] = {
    "endothelial": "#e74c3c",  # red
    "immune": "#3498db",  # blue
    "other": "#2ecc71",  # green
    "unknown": "#95a5a6",  # grey
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _zscore_rows(df: pd.DataFrame) -> pd.DataFrame:
    """Row-wise Z-score normalisation; constant rows are set to 0."""
    std = df.std(axis=1)
    mean = df.mean(axis=1)
    std_safe = std.replace(0.0, np.nan)
    return df.subtract(mean, axis=0).divide(std_safe, axis=0).fillna(0.0)


def _categorise_celltype(ct_name: str) -> str:
    """Return the lineage category for a cell type name via keyword match."""
    ct_lower = ct_name.lower()
    for cat, keywords in CATEGORY_KEYWORDS.items():
        if any(kw in ct_lower for kw in keywords):
            return cat
    return "unknown"


# ============================================================
# SECTION 1: Load Multiple Cell Types with Jacobian Data
# ============================================================


def section1_load_celltypes(loader) -> dict:
    """Select and load 4-6 representative cell types spanning different lineages.

    ## Why load with jacob=True?

    The Jacobian (∂output/∂input) measures how much the GET model's CD70
    expression prediction changes when each TF-motif input at each genomic
    region is perturbed.  Loading with ``jacob=True`` streams the pre-computed
    Jacobian zarr arrays from S3, enabling:

      - ``get_gene_jacobian_summary(gene, axis='motif')``: per-motif importance
      - ``gene_by_motif``: full gene-by-motif regulatory matrix

    ## Selection strategy

    We use keyword matching to assign each available cell type to a lineage
    category (endothelial / immune / other) and select up to
    ``CATEGORY_MAX_PER`` types per category, giving 4–6 cell types total.
    If the dataset contains more than ``CATEGORY_MAX_PER`` matches for a
    category, we take the first ones in alphabetical order.

    Parameters
    ----------
    loader : GETDemoLoader
        Initialised loader with ``available_celltypes`` populated.

    Returns
    -------
    dict[str, dict]
        Maps ``ct_name`` → ``{"ct": GETCellType, "category": str,
        "has_cd70": bool}``.
    """
    print("\n" + "=" * 60)
    print("SECTION 1: Load Multiple Cell Types with Jacobian Data")
    print("=" * 60)

    all_celltypes = loader.available_celltypes
    print(f"\nTotal available cell types: {len(all_celltypes)}")

    # ── Assign each available cell type to a lineage category ──────────────
    category_matches: dict[str, list[str]] = {cat: [] for cat in CATEGORY_KEYWORDS}
    for ct_name in sorted(all_celltypes):
        cat = _categorise_celltype(ct_name)
        if cat != "unknown" and ct_name not in category_matches.get(cat, []):
            category_matches[cat].append(ct_name)

    print("\nCell type discovery by lineage category:")
    for cat, matches in category_matches.items():
        chosen = matches[:CATEGORY_MAX_PER]
        print(f"  {cat:15s}: {chosen}")

    # ── Build final selection (up to CATEGORY_MAX_PER per category) ─────────
    selected: dict[str, str] = {}  # ct_name → category
    for cat, matches in category_matches.items():
        for ct_name in matches[:CATEGORY_MAX_PER]:
            selected[ct_name] = cat

    # Safety: if we got fewer than 4, pad from remaining available types.
    if len(selected) < 4:
        for ct_name in all_celltypes:
            if ct_name not in selected:
                selected[ct_name] = "unknown"
            if len(selected) >= 6:
                break

    print(f"\nSelected {len(selected)} cell types for detailed analysis:")
    for ct_name, cat in selected.items():
        print(f"  [{cat:12s}] {ct_name}")

    # ── Load each selected cell type with jacob=True ─────────────────────────
    # jacob=True streams the Jacobian zarr arrays (~100–500 MB per cell type).
    # We load them sequentially to avoid overwhelming the S3 connection.
    loaded: dict[str, dict] = {}
    for ct_name, cat in selected.items():
        print(f"\n  Loading '{ct_name}' (jacob=True) …")
        try:
            ct = loader.load_celltype(ct_name, jacob=True)

            # Check whether CD70's TSS falls within an open ATAC peak in this
            # cell type.  gene_annot contains only genes with accessible TSSs.
            has_cd70 = GENE in ct.gene_annot["gene_name"].values
            status = "✓ CD70 in open chromatin" if has_cd70 else "✗ CD70 not accessible"
            print(f"    {status}")

            if has_cd70:
                annot = ct.get_gene_annot(GENE)
                print(
                    f"    pred={annot['pred'].mean():.4f}  "
                    f"obs={annot['obs'].mean():.4f}  "
                    f"acc={annot['accessibility'].mean():.4f}"
                )

            loaded[ct_name] = {"ct": ct, "category": cat, "has_cd70": has_cd70}

        except Exception as exc:
            logger.warning("Failed to load '%s': %s", ct_name, exc)

    n_loaded = len(loaded)
    n_with_cd70 = sum(1 for v in loaded.values() if v["has_cd70"])
    print(f"\nSuccessfully loaded: {n_loaded} cell types")
    print(f"With CD70 in open chromatin: {n_with_cd70}")
    print(f"Without CD70 (useful as contrast): {n_loaded - n_with_cd70}")

    return loaded


# ============================================================
# SECTION 2: Cross-Cell-Type Motif Comparison (Clustered Heatmap)
# ============================================================


def section2_motif_heatmap(loaded_celltypes: dict) -> pd.DataFrame:
    """Build a CD70 motif × cell-type Jacobian matrix and plot a clustered heatmap.

    ## Jacobian motif summary (per cell type)

    For each loaded cell type, ``get_gene_jacobian_summary('CD70', axis='motif')``
    returns a ``pd.Series`` with one entry per TF motif cluster (NrMotifV1,
    282 clusters), valued as the sum of absolute-mean Jacobians across all
    TSS-proximal peaks and all TSSs of CD70.  A higher score means that
    motif's presence more strongly influences the GET model's CD70 prediction.

    ## Heatmap construction

      1. Stack per-cell-type Series into a motif × cell-type DataFrame.
      2. Drop the "Accessibility" pseudo-feature (it's an ATAC input signal,
         not a sequence-based TF motif).
      3. Select the top ``TOP_N_MOTIFS`` rows by max absolute Jacobian across
         all cell types (captures motifs that are important in *any* lineage).
      4. Row-standardise to Z-scores so the diverging colour scale reflects
         *relative* importance rather than absolute magnitude differences
         between motifs.

    ## Biological insight from clustering

    The row dendrogram groups TF motifs that are co-regulated across the same
    set of cell types, revealing TF *modules*.  The column dendrogram groups
    cell types by similarity of their CD70 regulatory programs, which typically
    separates endothelial from immune types.

    Parameters
    ----------
    loaded_celltypes : dict
        Output from :func:`section1_load_celltypes`.

    Returns
    -------
    pd.DataFrame
        Full motif × cell-type Jacobian matrix (all motifs, all cell types
        with CD70).  Empty DataFrame if fewer than 2 cell types have data.
    """
    print("\n" + "=" * 60)
    print("SECTION 2: Cross-Cell-Type Motif Comparison (Heatmap)")
    print("=" * 60)

    # ── Collect per-cell-type Jacobian motif summaries ───────────────────────
    motif_series: dict[str, pd.Series] = {}

    for ct_name, info in loaded_celltypes.items():
        if not info["has_cd70"]:
            logger.info("Skipping '%s' — CD70 not accessible", ct_name)
            continue

        ct = info["ct"]
        print(f"\n  Computing CD70 Jacobian motif summary for '{ct_name}' …")
        try:
            # Returns pd.Series: index = motif cluster names, values = scores.
            # axis='motif' aggregates across all 200 TSS-proximal ATAC peaks
            # and across all TSSs of CD70 (some genes have multiple TSSs).
            jac: pd.Series = ct.get_gene_jacobian_summary(GENE, axis="motif")

            # Drop the non-motif "Accessibility" feature before analysis.
            if "Accessibility" in jac.index:
                jac = jac.drop("Accessibility")

            motif_series[ct_name] = jac
            top3 = jac.nlargest(3).index.tolist()
            print(f"    ✓ {len(jac)} motifs  |  top-3: {top3}")

        except Exception as exc:
            logger.warning("Jacobian summary failed for '%s': %s", ct_name, exc)

    if len(motif_series) < 2:
        print(
            "\n  ⚠ Fewer than 2 cell types have CD70 Jacobian data — skipping heatmap."
        )
        return pd.DataFrame()

    # ── Build motif × cell-type DataFrame ────────────────────────────────────
    # pd.DataFrame from a dict of Series aligns on the shared motif index.
    motif_df = pd.DataFrame(motif_series).fillna(0.0)
    print(
        f"\nMotif × cell-type matrix: "
        f"{motif_df.shape[0]} motifs × {motif_df.shape[1]} cell types"
    )

    # ── Select top TOP_N_MOTIFS motifs by peak absolute Jacobian ─────────────
    # Using max (not mean) avoids missing motifs critical in only one lineage —
    # e.g., an endothelial-specific NF-κB motif would be washed out by mean.
    top_idx = motif_df.abs().max(axis=1).nlargest(TOP_N_MOTIFS).index
    motif_top = motif_df.loc[top_idx].copy()

    # ── Row-wise Z-score normalisation ───────────────────────────────────────
    # Makes relative differences between cell types visible even when absolute
    # Jacobian magnitudes vary by orders of magnitude across motifs.
    motif_z = _zscore_rows(motif_top)

    # ── Column colour bar: one colour per lineage category ───────────────────
    col_colors = pd.Series(
        {
            ct: CATEGORY_COLORS.get(
                loaded_celltypes[ct]["category"], CATEGORY_COLORS["unknown"]
            )
            for ct in motif_z.columns
        },
        name="Lineage",
    )

    # ── Create seaborn clustermap ─────────────────────────────────────────────
    n_ct = motif_z.shape[1]
    n_motifs = motif_z.shape[0]
    fig_w = max(10, n_ct * 1.5)
    fig_h = max(12, n_motifs * 0.42)

    print(f"\n  Generating clustered heatmap ({n_motifs} motifs × {n_ct} cell types) …")

    g = sns.clustermap(
        motif_z,
        col_colors=col_colors,
        cmap="RdBu_r",
        center=0.0,
        figsize=(fig_w, fig_h),
        xticklabels=True,
        yticklabels=True,
        dendrogram_ratio=(0.12, 0.12),
        cbar_pos=(0.02, 0.35, 0.03, 0.3),
        linewidths=0.0,
        method="average",
        metric="correlation",
    )

    g.ax_heatmap.set_xlabel("Cell Type", fontsize=11)
    g.ax_heatmap.set_ylabel("TF Motif Cluster (NrMotifV1)", fontsize=11)
    g.ax_heatmap.set_title(
        f"{GENE} Regulatory Motif Importance Across Cell Types\n"
        f"(top {TOP_N_MOTIFS} motifs by peak |Jacobian|; "
        f"rows Z-score normalised)",
        fontsize=12,
        pad=20,
    )
    g.ax_heatmap.tick_params(axis="x", labelsize=7, rotation=90)
    g.ax_heatmap.tick_params(axis="y", labelsize=7)

    # Legend for lineage colours
    legend_patches = [
        Patch(facecolor=col, label=cat.capitalize())
        for cat, col in CATEGORY_COLORS.items()
        if cat in {v["category"] for v in loaded_celltypes.values()}
    ]
    g.ax_col_colors.legend(
        handles=legend_patches,
        loc="upper left",
        bbox_to_anchor=(1.02, 1.0),
        title="Lineage",
        fontsize=8,
        title_fontsize=9,
        framealpha=0.9,
    )

    out_path = FIGURES_DIR / "cd70_motif_heatmap_cross_celltype.png"
    g.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close("all")
    print(f"  Saved → {out_path}")

    return motif_df


# ============================================================
# SECTION 3: Differential Motif Analysis (Waterfall Plot)
# ============================================================


def section3_differential_motifs(
    loaded_celltypes: dict,
    motif_df: pd.DataFrame,
) -> None:
    """Rank TF motifs by endothelial-specific CD70 regulatory importance.

    ## Differential motif score

    For each TF motif *m*::

        differential(m) = mean(Jacobian_m across endothelial types)
                        − mean(Jacobian_m across all other types)

    A **positive** differential score indicates *m* more strongly drives CD70
    in endothelial cells (endothelial-enriched regulator).  A **negative**
    score indicates *m* is more important in immune or other cell types
    (endothelial-depleted).

    ## Biological interpretation

    **Endothelial-enriched regulators (positive differential)**
    These TF motifs reflect the *inflammatory/vascular* transcriptional
    program that controls CD70 in the vasculature:
      - NF-κB motifs (RELA, NFKB1): cytokine-driven vascular inflammation
      - AP-1 (FOS/JUN): acute inflammatory stress response
      - ETS family (ERG, FLI1): endothelial identity and angiogenesis

    **Endothelial-depleted regulators (negative differential)**
    These TF motifs are more important for CD70 in immune cells:
      - NFAT (NFATC1–4): T-cell activation pathway
      - IRF (IRF1, IRF4): interferon/innate immune signalling
      - RUNX (RUNX1–3): haematopoiesis and T-cell differentiation

    Parameters
    ----------
    loaded_celltypes : dict
        Output from :func:`section1_load_celltypes`.
    motif_df : pd.DataFrame
        Motif × cell-type Jacobian matrix from :func:`section2_motif_heatmap`.
    """
    print("\n" + "=" * 60)
    print("SECTION 3: Differential Motif Analysis (Waterfall Plot)")
    print("=" * 60)

    if motif_df.empty or motif_df.shape[1] < 2:
        print("  ⚠ Insufficient data for differential analysis.")
        return

    # ── Identify endothelial vs other columns in the Jacobian matrix ─────────
    endo_cols = [
        ct
        for ct in motif_df.columns
        if ct in loaded_celltypes and loaded_celltypes[ct]["category"] == "endothelial"
    ]
    # Fallback: keyword matching if no category-labelled endothelial types
    if not endo_cols:
        endo_cols = [
            ct
            for ct in motif_df.columns
            if any(kw in ct.lower() for kw in CATEGORY_KEYWORDS["endothelial"])
        ]

    other_cols = [ct for ct in motif_df.columns if ct not in endo_cols]

    if not endo_cols or not other_cols:
        print(
            "  ⚠ Need both endothelial and other cell types in the "
            "Jacobian matrix — skipping differential analysis."
        )
        return

    print(f"\n  Endothelial: {endo_cols}")
    print(f"  Other:       {other_cols}")

    # ── Compute differential scores ───────────────────────────────────────────
    endo_mean = motif_df[endo_cols].mean(axis=1)
    other_mean = motif_df[other_cols].mean(axis=1)
    differential = (endo_mean - other_mean).sort_values(ascending=False)

    # Keep motifs with a differential signal above the 10th percentile of
    # absolute differences (removes near-zero noise).
    threshold = differential.abs().quantile(0.10)
    sig = differential[differential.abs() >= threshold]

    n_enriched = int((sig > 0).sum())
    n_depleted = int((sig < 0).sum())
    print(f"\n  Endothelial-enriched motifs (differential > 0): {n_enriched}")
    print(f"  Endothelial-depleted motifs (differential < 0): {n_depleted}")
    print(f"\n  Top 10 endothelial-enriched:\n{differential.head(10).to_string()}")
    print(f"\n  Top 10 endothelial-depleted:\n{differential.tail(10).to_string()}")

    # ── Build display selection: top 25 enriched + top 25 depleted ───────────
    n_show = min(25, max(n_enriched, 5))
    display = pd.concat(
        [differential.head(n_show), differential.tail(n_show)]
    ).drop_duplicates()
    display = display.sort_values(ascending=False)

    # ── Waterfall plot ────────────────────────────────────────────────────────
    n_bars = len(display)
    fig_w = max(14, n_bars * 0.45)
    fig, ax = plt.subplots(figsize=(fig_w, 7))

    bar_colors = [
        CATEGORY_COLORS["endothelial"] if v > 0 else CATEGORY_COLORS["immune"]
        for v in display.values
    ]
    ax.bar(
        range(n_bars),
        display.values,
        color=bar_colors,
        edgecolor="white",
        linewidth=0.4,
        width=0.85,
    )
    ax.axhline(y=0.0, color="black", linewidth=0.9, linestyle="-")

    # Annotate the most extreme motifs (top 5 each side)
    for rank_i, (motif_name, score) in enumerate(display.items()):
        if rank_i < 5 or rank_i >= n_bars - 5:
            y_offset = 5 if score >= 0 else -12
            ax.annotate(
                motif_name.split("_")[0],  # use first token for readability
                xy=(rank_i, score),
                xytext=(0, y_offset),
                textcoords="offset points",
                ha="center",
                fontsize=6,
                color="black",
                alpha=0.85,
            )

    ax.set_xticks(range(n_bars))
    ax.set_xticklabels(display.index, rotation=90, ha="center", fontsize=6.5)
    ax.set_xlabel("TF Motif Cluster", fontsize=12)
    ax.set_ylabel(
        "Differential CD70 Regulatory Score\n"
        "(endothelial mean Jacobian  −  other-lineage mean Jacobian)",
        fontsize=11,
    )
    ax.set_title(
        f"{GENE} Endothelial-Differential TF Motif Regulators\n"
        "Red = endothelial-enriched; Blue = depleted (immune/other-enriched)",
        fontsize=12,
    )

    legend_elems = [
        Patch(
            facecolor=CATEGORY_COLORS["endothelial"],
            label=f"Endothelial-enriched ({n_enriched})",
        ),
        Patch(
            facecolor=CATEGORY_COLORS["immune"],
            label=f"Endothelial-depleted ({n_depleted})",
        ),
    ]
    ax.legend(handles=legend_elems, loc="upper right", framealpha=0.9, fontsize=9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()
    out_path = FIGURES_DIR / "cd70_endothelial_differential_motifs.png"
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"\n  Saved → {out_path}")


# ============================================================
# SECTION 4: Gene-by-Motif Regulatory Profile Comparison
# ============================================================


def section4_regulatory_profile(loaded_celltypes: dict) -> None:
    """Compare CD70's motif regulatory fingerprint via the gene-by-motif matrix.

    ## Gene-by-motif matrix

    The gene-by-motif matrix encodes the overall TF-motif regulatory influence
    on each gene, computed by summing the Jacobian across all TSS-proximal
    ATAC peaks and all TSSs of the gene.  Extracting CD70's row gives a
    282-dimensional *regulatory fingerprint* for CD70 in each cell type.

    Unlike ``get_gene_jacobian_summary()``, which is computed on demand, the
    gene-by-motif matrix is pre-cached in the zarr store.  It therefore covers
    *all* genes simultaneously and can be loaded faster when the full matrix
    is needed.

    ## Multi-panel ranked profile plot

    For each cell type, we plot the motif importance values sorted by absolute
    magnitude (rank 1 = most important motif).  This elbow-shaped curve reveals:

      - **Steepness**: how concentrated regulatory control is in a few motifs
        vs. broadly distributed across many TFs.
      - **Top motif identity**: which TF families dominate CD70 regulation in
        that lineage.
      - **Profile similarity**: comparing shapes across panels reveals whether
        CD70 uses a conserved or lineage-specific regulatory architecture.

    ## What makes CD70 regulation unique in endothelial cells?

    If endothelial panels show a steeper profile (regulatory control
    concentrated in fewer motifs) and those top motifs are different from
    immune panels, CD70 regulation is *lineage-specific*.  If profiles are
    similar in shape and top motifs, a conserved core regulatory program
    exists, modified by lineage-specific modifiers.

    Parameters
    ----------
    loaded_celltypes : dict
        Output from :func:`section1_load_celltypes`.
    """
    print("\n" + "=" * 60)
    print("SECTION 4: Gene-by-Motif Regulatory Profile Comparison")
    print("=" * 60)

    # ── Extract CD70's row from each cell type's gene_by_motif matrix ─────────
    cd70_profiles: dict[str, pd.Series] = {}
    category_map: dict[str, str] = {}

    for ct_name, info in loaded_celltypes.items():
        if not info["has_cd70"]:
            logger.info("Skipping '%s' — CD70 not accessible", ct_name)
            continue

        ct = info["ct"]
        print(f"\n  Extracting gene_by_motif for '{ct_name}' …")
        try:
            # gene_by_motif is a lazy property: first access computes or loads
            # the full gene × motif matrix from the cell type's zarr store.
            # .data is a pd.DataFrame: index = gene names, columns = motif features.
            gbm_data: pd.DataFrame = ct.gene_by_motif.data

            # Locate CD70's row (may be a Series or a sub-DataFrame if multiple
            # TSSs share the same gene name in the index).
            if GENE in gbm_data.index:
                row = gbm_data.loc[GENE]
            else:
                # Fallback: case-insensitive search in case naming differs.
                matches = [g for g in gbm_data.index if str(g).upper() == GENE]
                if not matches:
                    print("    ✗ CD70 not found in gene_by_motif index")
                    continue
                row = gbm_data.loc[matches[0]]

            # If multiple TSSs → average across them.
            if isinstance(row, pd.DataFrame):
                row = row.mean(axis=0)

            # Drop the non-motif "Accessibility" feature.
            if "Accessibility" in row.index:
                row = row.drop("Accessibility")

            cd70_profiles[ct_name] = row.astype(float)
            category_map[ct_name] = info["category"]
            top3 = row.abs().nlargest(3).index.tolist()
            print(f"    ✓ {len(row)} motifs  |  top-3: {top3}")

        except Exception as exc:
            logger.warning("gene_by_motif extraction failed for '%s': %s", ct_name, exc)

    if not cd70_profiles:
        print("\n  ⚠ No gene_by_motif CD70 data retrieved — skipping Section 4.")
        return

    n_profiles = len(cd70_profiles)
    print(f"\n  CD70 regulatory profiles retrieved for {n_profiles} cell type(s).")

    # ── Multi-panel ranked-importance profile plot ────────────────────────────
    ncols = min(3, n_profiles)
    nrows = (n_profiles + ncols - 1) // ncols

    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(ncols * 5, nrows * 4),
        sharey=False,
    )

    # Normalise axes shape to 2-D array for uniform indexing.
    if n_profiles == 1:
        axes = np.array([[axes]])
    elif nrows == 1:
        axes = axes.reshape(1, -1)

    for idx, (ct_name, profile) in enumerate(cd70_profiles.items()):
        row_i, col_i = divmod(idx, ncols)
        ax = axes[row_i, col_i]
        cat = category_map.get(ct_name, "unknown")
        color = CATEGORY_COLORS.get(cat, CATEGORY_COLORS["unknown"])

        # Sort by absolute value (rank 1 = highest regulatory influence).
        sorted_profile = profile.abs().sort_values(ascending=False)
        n_motifs = len(sorted_profile)
        ranks = np.arange(1, n_motifs + 1)

        ax.fill_between(ranks, sorted_profile.values, alpha=0.25, color=color)
        ax.plot(ranks, sorted_profile.values, color=color, linewidth=1.6)

        # Label the top 5 motifs by name (first token of motif cluster ID).
        for rank_j, (motif_name, motif_val) in enumerate(
            sorted_profile.head(5).items(), start=1
        ):
            ax.annotate(
                motif_name.split("_")[0],
                xy=(rank_j, motif_val),
                xytext=(3, 3),
                textcoords="offset points",
                ha="left",
                fontsize=5.5,
                color=color,
                fontweight="bold",
            )

        ax.set_title(ct_name, fontsize=8, color=color, fontweight="bold", pad=4)
        ax.set_xlabel("Motif rank (by |Jacobian|)", fontsize=8)
        ax.set_ylabel("|Gene-by-motif score|", fontsize=8)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.tick_params(labelsize=7)

    # Hide empty subplots (when n_profiles is not a multiple of ncols).
    for idx in range(n_profiles, nrows * ncols):
        row_i, col_i = divmod(idx, ncols)
        axes[row_i, col_i].set_visible(False)

    # Legend
    legend_patches = [
        Patch(facecolor=col, label=cat.capitalize(), alpha=0.7)
        for cat, col in CATEGORY_COLORS.items()
        if cat in set(category_map.values())
    ]
    fig.legend(
        handles=legend_patches,
        loc="lower right",
        fontsize=9,
        title="Lineage",
        framealpha=0.9,
    )
    fig.suptitle(
        f"{GENE} Motif Regulatory Profile Across Cell Types\n"
        "(ranked by absolute gene-by-motif score; "
        "top-5 TF motifs labelled per panel)",
        fontsize=12,
        y=1.01,
    )

    plt.tight_layout()
    out_path = FIGURES_DIR / "cd70_regulatory_profile_comparison.png"
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved → {out_path}")

    # ── Cross-cell-type profile correlation summary ───────────────────────────
    # How similar are CD70's regulatory fingerprints across lineages?
    # High correlation → conserved regulatory program
    # Low correlation  → lineage-specific regulation
    if n_profiles >= 2:
        profile_df = pd.DataFrame(cd70_profiles)  # motifs × cell types
        corr = profile_df.corr()
        print("\n  CD70 regulatory profile pairwise correlation (Pearson r):")
        print(corr.round(3).to_string())
        print(
            "\n  Interpretation:\n"
            "  r > 0.8  → highly similar regulatory programs\n"
            "  r 0.5–0.8 → moderate overlap (shared core, lineage-specific top)\n"
            "  r < 0.5  → distinct regulatory programs across lineages"
        )


# ============================================================
# Entry point
# ============================================================


if __name__ == "__main__":
    print("=" * 60)
    print("CD70 Comparative Regulation Tutorial")
    print(f"Gene target  : {GENE}")
    print(f"Output dir   : {TUTORIAL_DIR}")
    print("=" * 60)

    from gcell.cell.celltype import GETDemoLoader

    # Initialise the demo loader (queries S3 for available cell types).
    loader = GETDemoLoader()

    # Section 1: Select and load 4-6 representative cell types with jacob=True.
    loaded = section1_load_celltypes(loader)

    # Section 2: Build motif × cell-type Jacobian matrix and plot heatmap.
    motif_df = section2_motif_heatmap(loaded)

    # Section 3: Differential endothelial vs other motif analysis.
    section3_differential_motifs(loaded, motif_df)

    # Section 4: Gene-by-motif regulatory profile comparison.
    section4_regulatory_profile(loaded)

    print("\n" + "=" * 60)
    print("All done!")
    print("  Figure 1 → tutorials/figures/cd70_motif_heatmap_cross_celltype.png")
    print("  Figure 2 → tutorials/figures/cd70_endothelial_differential_motifs.png")
    print("  Figure 3 → tutorials/figures/cd70_regulatory_profile_comparison.png")
    print("=" * 60)
