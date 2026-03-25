"""CD70 Cross-Cell-Type Regulatory Analysis Tutorial.

This tutorial demonstrates how to use GETDemoLoader to analyze CD70 gene
chromatin accessibility and expression patterns across all available cell
types using the GET (Gene Expression Transformer) pre-inferred model.

CD70 (also known as TNFSF7) is a TNF superfamily member that acts as a
co-stimulatory ligand for CD27.  Its expression is tightly regulated and
varies considerably across hematopoietic and non-hematopoietic cell types,
making it an interesting target for cross-cell-type regulatory analysis.

The tutorial is divided into three sections:

  Section 1 – Setup & Discovery
      Initialize GETDemoLoader, print all available cell types, identify
      endothelial cell types, and verify that CD70 is accessible in a
      representative endothelial cell type.

  Section 2 – Cross-cell-type CD70 Accessibility & Expression
      Loop over all available cell types, retrieve CD70 accessibility,
      predicted expression, and observed expression from the GET model,
      collect the results into a DataFrame, and save it to CSV.

  Section 3 – Visualization
      (a) Bar chart of CD70 accessibility across cell types, sorted by
          value, endothelial cells highlighted in red.
      (b) Scatter plot of CD70 predicted vs observed expression across
          cell types, coloured by accessibility level.

Usage
-----
Run from the repository root::

    uv run python tutorials/cd70_regulation_analysis.py

Output files
------------
tutorials/cd70_cross_celltype_summary.csv
tutorials/figures/cd70_accessibility_barplot.png
tutorials/figures/cd70_pred_vs_obs.png
"""

import logging
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")  # non-interactive backend — must be set before pyplot import
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
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

#: Lower-case keywords used to identify endothelial cell types.
ENDOTHELIAL_KEYWORDS = ["endothelial", "huvec", "ec ", "endoth"]


# ============================================================
# SECTION 1: Setup & Discovery
# ============================================================


def section1_setup():
    """Initialize GETDemoLoader and identify endothelial cell types.

    Steps
    -----
    1. Construct a :class:`~gcell.cell.celltype.GETDemoLoader` which fetches
       the list of available pre-inferred cell types from S3.
    2. Print the complete list of available cell types.
    3. Filter the list to retain only endothelial-related entries using a
       set of case-insensitive keywords.
    4. Load a representative endothelial cell type (first match, or the
       first available cell type when no endothelial types are found) and
       verify that CD70 appears in its open-chromatin gene annotation.

    Returns
    -------
    loader : GETDemoLoader
        Initialised loader with ``available_celltypes`` populated.
    endothelial_celltypes : list[str]
        Subset of ``available_celltypes`` matching endothelial keywords.
    """
    print("\n" + "=" * 60)
    print("SECTION 1: Setup & Discovery")
    print("=" * 60)

    from gcell.cell.celltype import GETDemoLoader

    # -- 1. Initialise loader -----------------------------------------------
    # GETDemoLoader queries S3 for the interpretation directory and builds
    # available_celltypes from the folder names found there.
    loader = GETDemoLoader()

    total = len(loader.available_celltypes)
    print(f"\nTotal available cell types: {total}")
    print("\nAll available cell types:")
    for ct in loader.available_celltypes:
        print(f"  - {ct}")

    # -- 2. Identify endothelial cell types ----------------------------------
    endothelial_celltypes = [
        ct
        for ct in loader.available_celltypes
        if any(kw in ct.lower() for kw in ENDOTHELIAL_KEYWORDS)
    ]

    print(f"\nEndothelial-related cell types ({len(endothelial_celltypes)}):")
    if endothelial_celltypes:
        for ct in endothelial_celltypes:
            print(f"  - {ct}")
    else:
        print("  (none found with the current keyword list)")

    # -- 3. Load representative cell type & verify CD70 ---------------------
    if endothelial_celltypes:
        representative = endothelial_celltypes[0]
    else:
        representative = loader.available_celltypes[0]
        print(
            f"\nNo endothelial cell type found — using '{representative}' as representative."
        )

    print(f"\nLoading '{representative}' to verify CD70 accessibility …")
    ct = loader.load_celltype(representative)

    gene_names = ct.gene_annot["gene_name"].values

    if GENE in gene_names:
        annot = ct.get_gene_annot(GENE)
        print(f"  ✓ {GENE} found in '{representative}':")
        print(f"    Chromosomes : {annot['Chromosome'].values}")
        print(f"    Accessibility : {annot['accessibility'].values}")
        print(f"    Pred expression : {annot['pred'].values}")
        print(f"    Obs  expression : {annot['obs'].values}")
    else:
        print(
            f"  ✗ {GENE} is NOT in '{representative}' open-chromatin gene annotation."
        )
        print("    (The promoter region may not be accessible in this cell type.)")

    return loader, endothelial_celltypes


# ============================================================
# SECTION 2: Cross-cell-type CD70 Accessibility & Expression
# ============================================================


def section2_cross_celltype(loader, endothelial_celltypes):
    """Collect CD70 accessibility and expression data across all cell types.

    For each cell type in ``loader.available_celltypes``:

    1. Load the cell type via :meth:`~GETDemoLoader.load_celltype`.
    2. Check whether CD70 is present in the cell type's open-chromatin gene
       annotation (``gene_annot``).  If not, record ``NaN`` values and
       ``has_cd70=False``.
    3. If CD70 is present, call:
       - :meth:`~Celltype.get_gene_accessibility` — ATAC signal at TSS.
       - :meth:`~Celltype.get_gene_pred` — GET model predicted expression.
       - :meth:`~Celltype.get_gene_obs` — observed RNA expression.
       - :meth:`~Celltype.get_gene_annot` — full per-TSS annotation.
    4. Average all values across TSSs (most genes have one, some have two).

    Parameters
    ----------
    loader : GETDemoLoader
        Initialised loader.
    endothelial_celltypes : list[str]
        Cell types identified as endothelial in Section 1.

    Returns
    -------
    pd.DataFrame
        One row per cell type with columns:
        ``celltype``, ``accessibility``, ``pred_exp``, ``obs_exp``,
        ``has_cd70``, ``is_endothelial``.
        Also written to ``tutorials/cd70_cross_celltype_summary.csv``.
    """
    print("\n" + "=" * 60)
    print("SECTION 2: Cross-cell-type CD70 Accessibility & Expression")
    print("=" * 60)

    all_celltypes = loader.available_celltypes
    print(f"\nAnalysing {GENE} across {len(all_celltypes)} cell types …\n")

    records = []

    for i, ct_name in enumerate(all_celltypes):
        prefix = f"  [{i + 1:3d}/{len(all_celltypes)}] {ct_name:<45s}"

        try:
            ct = loader.load_celltype(ct_name)

            # ----------------------------------------------------------------
            # Check whether CD70 is in this cell type's open chromatin.
            # gene_annot contains only genes whose TSS falls within an open
            # ATAC-seq peak, so absence here means the promoter is inaccessible.
            # ----------------------------------------------------------------
            if GENE not in ct.gene_annot["gene_name"].values:
                print(f"{prefix} → {GENE} not in open chromatin")
                records.append(
                    {
                        "celltype": ct_name,
                        "accessibility": np.nan,
                        "pred_exp": np.nan,
                        "obs_exp": np.nan,
                        "has_cd70": False,
                        "is_endothelial": ct_name in endothelial_celltypes,
                    }
                )
                continue

            # ----------------------------------------------------------------
            # Retrieve per-TSS annotation DataFrame.
            # gene_annot already contains pre-computed 'accessibility', 'pred',
            # and 'obs' columns, providing a reliable fallback for all metrics.
            # ----------------------------------------------------------------
            annot = ct.get_gene_annot(GENE)

            # Accessibility ––––––––––––––––––––––––––––––––––––––––––––––––
            # Prefer get_gene_accessibility() which returns the raw
            # tss_accessibility sparse matrix slice (available when
            # cfg.celltype.input=True, which is the GETDemoLoader default).
            # Fall back to the pre-computed column in gene_annot otherwise.
            acc_raw = ct.get_gene_accessibility(GENE)
            if acc_raw is not None and hasattr(acc_raw, "__len__") and len(acc_raw) > 0:
                try:
                    # tss_accessibility is a sparse matrix; dense-ify it.
                    acc_vals = np.asarray(acc_raw.toarray()).flatten()
                except AttributeError:
                    acc_vals = np.asarray(acc_raw).flatten()
                accessibility = float(np.mean(acc_vals))
            else:
                # Fallback: use the 'accessibility' column already in gene_annot
                accessibility = float(annot["accessibility"].mean())

            # Predicted expression –––––––––––––––––––––––––––––––––––––––––
            # get_gene_pred returns preds[gene_idx]: a 1-D float array with
            # one element per TSS.  Mean across TSSs gives a single summary.
            pred_arr = ct.get_gene_pred(GENE)
            pred_exp = float(np.mean(pred_arr)) if len(pred_arr) > 0 else np.nan

            # Observed expression ––––––––––––––––––––––––––––––––––––––––––
            # get_gene_obs returns obs[gene_idx].  Note: not all cell types
            # have measured RNA-seq; in those cases obs is set to 0.0.
            obs_arr = ct.get_gene_obs(GENE)
            obs_exp = float(np.mean(obs_arr)) if len(obs_arr) > 0 else np.nan

            print(
                f"{prefix} → acc={accessibility:.4f}  "
                f"pred={pred_exp:.4f}  obs={obs_exp:.4f}"
            )
            records.append(
                {
                    "celltype": ct_name,
                    "accessibility": accessibility,
                    "pred_exp": pred_exp,
                    "obs_exp": obs_exp,
                    "has_cd70": True,
                    "is_endothelial": ct_name in endothelial_celltypes,
                }
            )

        except Exception as exc:
            logger.warning(f"Failed to process '{ct_name}': {exc}")
            records.append(
                {
                    "celltype": ct_name,
                    "accessibility": np.nan,
                    "pred_exp": np.nan,
                    "obs_exp": np.nan,
                    "has_cd70": False,
                    "is_endothelial": ct_name in endothelial_celltypes,
                }
            )

    # ----------------------------------------------------------------
    # Assemble DataFrame and print summary statistics
    # ----------------------------------------------------------------
    df = pd.DataFrame(records)

    n_with = int(df["has_cd70"].sum())
    print("\nResults summary:")
    print(f"  Total cell types analysed : {len(df)}")
    print(f"  With CD70 in open chromatin : {n_with}")
    print(f"  Without CD70               : {len(df) - n_with}")

    if n_with > 0:
        stats = df[df["has_cd70"]][["accessibility", "pred_exp", "obs_exp"]].describe()
        print(f"\nDescriptive statistics (cell types with CD70):\n{stats}")

    # ----------------------------------------------------------------
    # Save to CSV
    # ----------------------------------------------------------------
    csv_path = TUTORIAL_DIR / "cd70_cross_celltype_summary.csv"
    df.to_csv(csv_path, index=False)
    print(f"\nSaved summary CSV → {csv_path}")

    return df


# ============================================================
# SECTION 3: Visualization
# ============================================================


def section3_visualization(df):
    """Generate two publication-style figures from the cross-cell-type data.

    Figure 1 — CD70 Accessibility Bar Chart
        Horizontal bar chart of mean TSS ATAC signal per cell type, sorted
        from highest to lowest accessibility.  Endothelial cell types are
        highlighted in red (#e74c3c); other cell types are shown in blue
        (#3498db).  A dashed horizontal line marks the cross-cell-type mean.
        Saved to ``tutorials/figures/cd70_accessibility_barplot.png``.

    Figure 2 — CD70 Predicted vs Observed Expression Scatter
        Scatter plot of GET-predicted vs RNA-observed CD70 expression across
        cell types.  Each point is coloured by its TSS accessibility score
        (viridis colour map).  Endothelial cell types are labelled by name
        in red.  A dashed diagonal (y = x) indicates perfect agreement
        between model prediction and observation.
        Saved to ``tutorials/figures/cd70_pred_vs_obs.png``.

    Parameters
    ----------
    df : pd.DataFrame
        Output DataFrame from :func:`section2_cross_celltype`.
    """
    print("\n" + "=" * 60)
    print("SECTION 3: Visualization")
    print("=" * 60)

    # Keep only cell types where CD70 is in open chromatin
    df_plot = df[df["has_cd70"]].dropna(subset=["accessibility"]).copy()

    if df_plot.empty:
        print("No cell types have CD70 in open chromatin — skipping plots.")
        return

    print(f"\nGenerating figures for {len(df_plot)} cell types with CD70 …")

    # ----------------------------------------------------------------
    # Figure 1: Bar chart — CD70 accessibility across cell types
    # ----------------------------------------------------------------
    df_sorted = df_plot.sort_values("accessibility", ascending=False).reset_index(
        drop=True
    )
    n = len(df_sorted)

    fig_width = max(12, n * 0.5)  # scale width with number of bars
    fig, ax = plt.subplots(figsize=(fig_width, 7))

    bar_colors = [
        "#e74c3c" if is_endo else "#3498db" for is_endo in df_sorted["is_endothelial"]
    ]
    ax.bar(
        range(n),
        df_sorted["accessibility"],
        color=bar_colors,
        edgecolor="white",
        linewidth=0.5,
    )

    # Mean reference line
    mean_acc = float(df_sorted["accessibility"].mean())
    ax.axhline(
        y=mean_acc,
        color="#7f8c8d",
        linestyle="--",
        linewidth=1.2,
        label=f"Mean = {mean_acc:.3f}",
    )

    ax.set_xticks(range(n))
    ax.set_xticklabels(df_sorted["celltype"], rotation=90, ha="right", fontsize=7)
    ax.set_xlabel("Cell Type", fontsize=12)
    ax.set_ylabel("CD70 TSS Accessibility\n(mean ATAC signal)", fontsize=12)
    ax.set_title(
        f"{GENE} Chromatin Accessibility Across Cell Types\n"
        "(sorted by accessibility; red = endothelial, blue = other)",
        fontsize=13,
    )

    legend_elements = [
        Patch(facecolor="#e74c3c", label="Endothelial"),
        Patch(facecolor="#3498db", label="Other cell types"),
    ]
    ax.legend(handles=legend_elements, loc="upper right", framealpha=0.9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()
    bar_path = FIGURES_DIR / "cd70_accessibility_barplot.png"
    fig.savefig(bar_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved Figure 1 → {bar_path}")

    # ----------------------------------------------------------------
    # Figure 2: Scatter — predicted vs observed expression
    # ----------------------------------------------------------------
    df_scatter = df_plot.dropna(subset=["pred_exp", "obs_exp", "accessibility"])

    if df_scatter.empty:
        print("  No complete pred/obs/accessibility data — skipping scatter plot.")
        return

    fig, ax = plt.subplots(figsize=(8, 7))

    sc = ax.scatter(
        df_scatter["obs_exp"],
        df_scatter["pred_exp"],
        c=df_scatter["accessibility"],
        cmap="viridis",
        s=90,
        alpha=0.85,
        edgecolors="none",
        zorder=3,
    )
    cbar = plt.colorbar(sc, ax=ax)
    cbar.set_label("CD70 TSS Accessibility (mean ATAC)", fontsize=10)

    # Label endothelial cell types
    for _, row in df_scatter[df_scatter["is_endothelial"]].iterrows():
        ax.annotate(
            row["celltype"],
            xy=(row["obs_exp"], row["pred_exp"]),
            fontsize=6.5,
            xytext=(5, 5),
            textcoords="offset points",
            color="#c0392b",
            fontweight="bold",
        )

    # Diagonal reference line (y = x → perfect prediction)
    all_vals = pd.concat([df_scatter["obs_exp"], df_scatter["pred_exp"]])
    margin = (all_vals.max() - all_vals.min()) * 0.05
    lim = (all_vals.min() - margin, all_vals.max() + margin)
    ax.plot(lim, lim, "k--", linewidth=1, alpha=0.4, label="y = x", zorder=2)
    ax.set_xlim(lim)
    ax.set_ylim(lim)

    ax.set_xlabel("Observed CD70 Expression", fontsize=12)
    ax.set_ylabel("Predicted CD70 Expression", fontsize=12)
    ax.set_title(
        f"{GENE} Predicted vs Observed Expression Across Cell Types\n"
        "(colour = TSS accessibility; red labels = endothelial)",
        fontsize=13,
    )
    ax.legend(loc="upper left", framealpha=0.9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()
    scatter_path = FIGURES_DIR / "cd70_pred_vs_obs.png"
    fig.savefig(scatter_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved Figure 2 → {scatter_path}")


# ============================================================
# Entry point
# ============================================================

if __name__ == "__main__":
    print("=" * 60)
    print("CD70 Cross-Cell-Type Regulatory Analysis")
    print(f"Gene target  : {GENE}")
    print(f"Output dir   : {TUTORIAL_DIR}")
    print("=" * 60)

    # Section 1: Initialise loader, identify endothelial cell types,
    #             verify CD70 accessibility in a representative cell type.
    loader, endothelial_celltypes = section1_setup()

    # Section 2: Collect CD70 data across all available cell types and
    #             save the results to CSV.
    df = section2_cross_celltype(loader, endothelial_celltypes)

    # Section 3: Generate accessibility bar chart and pred-vs-obs scatter plot.
    section3_visualization(df)

    print("\n" + "=" * 60)
    print("All done!")
    print("  CSV    → tutorials/cd70_cross_celltype_summary.csv")
    print("  Plot 1 → tutorials/figures/cd70_accessibility_barplot.png")
    print("  Plot 2 → tutorials/figures/cd70_pred_vs_obs.png")
    print("=" * 60)
