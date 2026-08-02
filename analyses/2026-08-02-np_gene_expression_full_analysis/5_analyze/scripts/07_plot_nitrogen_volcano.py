"""Volcano plot: all 53 usable genes from the original 54-gene N/P list,
differential expression under nitrogen starvation only, faceted by
platform (RNASEQ/PROTEOMICS/MICROARRAY -- padj granularity differs
sharply by platform, see caveat below), colored by each gene's
*original* researcher-assigned N/P annotation (not by what this
analysis found -- the color is the source-spreadsheet label).

x = log2FC (treatment vs control), y = -log10(padj). padj=0 rows
(floating-point underflow in the source study, 16 rows, all RNASEQ)
are floored to one order of magnitude below the smallest observed
nonzero padj (so they plot as the most-significant points, not tied at
an arbitrary mid-range value -- flooring to a fixed constant like 1e-12
would have placed them *below* several genuinely-observed p-values).
The top 3 points per panel by |log2FC| among significant (padj<0.05)
genes are labeled with their source gene name -- chosen over "top by
-log10(padj)" because many points tie at the padj floor, making that
selection arbitrary; |log2FC| ranking is unambiguous.

Caveat surfaced directly in the figure: MICROARRAY padj values in this
KG build are coarsely discretized (predominantly 0.01 or 1.0 -- see
5_analyze/notebook.md) rather than continuous, so the microarray panel
reads as two horizontal bands rather than a smooth cloud; RNASEQ and
PROTEOMICS have fine-grained padj and look like conventional volcano
plots.

Input:  data/06_nitrogen_volcano_data.csv
Output: figures/04_nitrogen_volcano.png

Usage:
  uv run python analyses/2026-08-02-np_gene_expression_full_analysis/5_analyze/scripts/07_plot_nitrogen_volcano.py
"""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

STEP_DIR = Path(__file__).resolve().parents[1]
DATA_DIR = STEP_DIR / "data"
FIG_DIR = STEP_DIR / "figures"

N_COLOR = "#2a78d6"   # blue -- categorical slot 1
P_COLOR = "#eb6834"   # orange -- categorical slot 8
INK = "#0b0b0b"
MUTED = "#898781"
GRID = "#e1e0d9"
PLATFORM_ORDER = ["RNASEQ", "PROTEOMICS", "MICROARRAY"]


def main() -> None:
    df = pd.read_csv(DATA_DIR / "06_nitrogen_volcano_data.csv")
    n_zero = int((df["padj"] == 0).sum())
    padj_floor = df.loc[df["padj"] > 0, "padj"].min() / 10
    df["padj_floored"] = df["padj"].clip(lower=padj_floor)
    df["neg_log10_padj"] = -np.log10(df["padj_floored"])

    fig, axes = plt.subplots(1, 3, figsize=(14, 5.5), dpi=300, sharey=True)

    for ax, platform in zip(axes, PLATFORM_ORDER):
        sub = df[df["omics_type"] == platform]
        for label, color in [("N", N_COLOR), ("P", P_COLOR)]:
            g = sub[sub["n_or_p"] == label]
            ax.scatter(
                g["log2fc"], g["neg_log10_padj"], s=16, alpha=0.65,
                color=color, edgecolors="none", label=f"{label}-annotated (n={len(g)})",
            )

        ax.axhline(-np.log10(0.05), color=MUTED, linewidth=1, linestyle="--")
        ax.axvline(0, color=MUTED, linewidth=0.8, linestyle=":")

        sig = sub[sub["padj"] < 0.05]
        top = sig.reindex(sig["log2fc"].abs().sort_values(ascending=False).index).head(3)
        for i, (_, r) in enumerate(top.iterrows()):
            ax.annotate(
                r["gene_name_raw"], xy=(r["log2fc"], r["neg_log10_padj"]),
                xytext=(4, 3 + 11 * i), textcoords="offset points", fontsize=7.5, color=INK,
            )

        ax.set_title(f"{platform}", fontsize=10, color=INK)
        ax.set_xlabel("log2 fold change", color=INK)
        for spine in ("top", "right"):
            ax.spines[spine].set_visible(False)
        ax.grid(color=GRID, linewidth=0.8, zorder=0)
        ax.set_axisbelow(True)
        ax.tick_params(colors=MUTED)

    axes[0].set_ylabel("-log10(padj)", color=INK)
    axes[0].legend(loc="lower right", frameon=False, fontsize=8)

    fig.suptitle(
        "Nitrogen-starvation DE, all 53 genes, colored by original N/P annotation (not by result)\n"
        f"Dashed line: padj=0.05. {n_zero} padj=0 points (RNASEQ) floored to {padj_floor:.1e}, just below "
        "the smallest observed nonzero padj. Microarray padj is coarsely discretized in this KG build "
        "(mostly 0.01/1.0), unlike RNA-seq/proteomics.",
        fontsize=9, color=INK, x=0.01, y=0.98, ha="left",
    )
    fig.tight_layout(rect=(0, 0, 1, 0.86))

    out_path = FIG_DIR / "04_nitrogen_volcano.png"
    fig.savefig(out_path, dpi=300)
    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
