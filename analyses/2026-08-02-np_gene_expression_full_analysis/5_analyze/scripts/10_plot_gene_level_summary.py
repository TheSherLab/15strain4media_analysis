"""One point per gene per platform (not per measurement) -- answers the
researcher's question directly: the earlier volcano/strip plots showed
one point per (gene x experiment x timepoint), so a heavily-measured
gene like ntcA (19 nitrogen rows) visually dominated a lightly-measured
gene (2-3 rows) despite both being "one gene." This figure collapses
that: x = each gene's median log2FC across all its measurements in a
platform, y = fraction of its measurements that were significant,
point size = how many measurements back it (so the weighting is shown,
not hidden), colored by original N/P annotation. Faceted by platform
(magnitudes still not pooled across platforms).

Input:  data/09_gene_level_summary.csv
Output: figures/06_gene_level_summary.png

Usage:
  uv run python analyses/2026-08-02-np_gene_expression_full_analysis/5_analyze/scripts/10_plot_gene_level_summary.py
"""
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

STEP_DIR = Path(__file__).resolve().parents[1]
DATA_DIR = STEP_DIR / "data"
FIG_DIR = STEP_DIR / "figures"

N_COLOR = "#2a78d6"
P_COLOR = "#eb6834"
INK = "#0b0b0b"
MUTED = "#898781"
GRID = "#e1e0d9"
PLATFORM_ORDER = ["RNASEQ", "PROTEOMICS", "MICROARRAY"]


def main() -> None:
    df = pd.read_csv(DATA_DIR / "09_gene_level_summary.csv")

    fig, axes = plt.subplots(1, 3, figsize=(14, 5.5), dpi=300, sharey=True)

    for ax, platform in zip(axes, PLATFORM_ORDER):
        sub = df[df["omics_type"] == platform]
        for label, color in [("N", N_COLOR), ("P", P_COLOR)]:
            g = sub[sub["n_or_p"] == label]
            ax.scatter(
                g["median_log2fc"], g["fraction_significant"],
                s=g["n_measurements"] * 12, alpha=0.55, color=color,
                edgecolors="none", label=f"{label}-annotated (n={len(g)} genes)",
            )

        ax.axvline(0, color=MUTED, linewidth=0.8, linestyle=":")
        ax.set_title(f"{platform}  ({sub['n_measurements'].sum()} total measurements)", fontsize=9.5, color=INK)
        ax.set_xlabel("median log2FC (across all timepoints)", color=INK)
        for spine in ("top", "right"):
            ax.spines[spine].set_visible(False)
        ax.grid(color=GRID, linewidth=0.8, zorder=0)
        ax.set_axisbelow(True)
        ax.tick_params(colors=MUTED)
        ax.legend(loc="upper left", frameon=False, fontsize=8)

    axes[0].set_ylabel("fraction of timepoints significant", color=INK)
    axes[0].set_ylim(-0.05, 1.05)

    fig.suptitle(
        "One point per gene (not per measurement) -- nitrogen starvation, all 53 genes\n"
        "Point size = number of (experiment x timepoint) measurements behind that gene's summary "
        "(e.g. ntcA=19 in RNASEQ+PROTEOMICS+MICROARRAY combined) -- shows the weighting instead of hiding it.",
        fontsize=9, color=INK, x=0.01, y=0.99, ha="left",
    )
    fig.tight_layout(rect=(0, 0, 1, 0.86))

    out_path = FIG_DIR / "06_gene_level_summary.png"
    fig.savefig(out_path, dpi=300)
    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
