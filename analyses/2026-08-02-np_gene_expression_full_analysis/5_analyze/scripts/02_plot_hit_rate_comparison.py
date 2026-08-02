"""Grouped bar figure comparing significant-hit rates (up/down/not
significant) across all 8 group x experiment-set combinations from
01_compute_hit_rates.py, grouped by which experiment set was queried
(nitrogen experiments on top, phosphorus experiments below) so the
target-vs-positive-control-vs-background comparison within each nutrient
side is visually adjacent.

Colors: blue = significant_up, red = significant_down (validated diverging
pair, see dataviz skill), light gray = not significant (deliberately
recedes; each segment carries a direct value label as the required
relief for its sub-3:1 contrast).

Input:  data/01_hit_rate_summary.csv
Output: figures/01_hit_rate_comparison.png

Usage:
  uv run python analyses/2026-08-02-np_gene_expression_full_analysis/5_analyze/scripts/02_plot_hit_rate_comparison.py
"""
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

STEP_DIR = Path(__file__).resolve().parents[1]
DATA_DIR = STEP_DIR / "data"
FIG_DIR = STEP_DIR / "figures"

BLUE = "#2a78d6"
RED = "#e34948"
GRAY = "#c3c2b7"
INK = "#0b0b0b"
MUTED = "#898781"

# Display order (nitrogen-experiment groups first, then phosphorus-experiment groups)
GROUP_ORDER = [
    "H1_N",
    "H2_Pgenes_Nexp",
    "positive_control_N",
    "background_N",
    "H1_P",
    "H2_Ngenes_Pexp",
    "positive_control_P",
]
LABELS = {
    "H1_N": "N genes x N starvation (matched, H1)",
    "H2_Pgenes_Nexp": "P genes x N starvation (cross, H2)",
    "positive_control_N": "N positive controls x N starvation",
    "background_N": "Background genes x N starvation (H3)",
    "H1_P": "P genes x P starvation (matched, H1)*",
    "H2_Ngenes_Pexp": "N genes x P starvation (cross, H2)*",
    "positive_control_P": "P positive controls x P starvation*",
}


def main() -> None:
    df = pd.read_csv(DATA_DIR / "01_hit_rate_summary.csv").set_index("group").loc[GROUP_ORDER]

    fig, ax = plt.subplots(figsize=(10.5, 5.8), dpi=300)
    y = range(len(GROUP_ORDER))[::-1]

    for yi, group in zip(y, GROUP_ORDER):
        row = df.loc[group]
        up, down, ns = row["rate_up"], row["rate_down"], row["n_not_significant"] / row["n_tests"]
        ax.barh(yi, up, color=BLUE, height=0.6)
        ax.barh(yi, down, left=up, color=RED, height=0.6)
        ax.barh(yi, ns, left=up + down, color=GRAY, height=0.6)
        ax.text(up + down + ns + 0.02, yi, f"n={int(row['n_tests'])}", va="center", fontsize=8, color=MUTED)
        if up > 0.03:
            ax.text(up / 2, yi, f"{up:.0%}", va="center", ha="center", fontsize=8, color="white")
        if down > 0.03:
            ax.text(up + down / 2, yi, f"{down:.0%}", va="center", ha="center", fontsize=8, color="white")

    ax.set_yticks(list(y))
    ax.set_yticklabels([LABELS[g] for g in GROUP_ORDER], fontsize=9, color=INK)
    ax.set_xlim(0, 1.12)
    ax.set_xlabel("Fraction of (gene x experiment x timepoint) tests", color=INK)
    ax.set_title(
        "Significant-hit rate: blue = up, red = down, gray = not significant\n"
        "* phosphorus tables are pre-filtered to already-significant genes (no phosphorus H3)",
        fontsize=9.5,
        color=INK,
        loc="left",
    )
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.tick_params(axis="x", colors=MUTED)
    fig.tight_layout()

    out_path = FIG_DIR / "01_hit_rate_comparison.png"
    fig.savefig(out_path, dpi=300)
    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
