"""Log2FC magnitude, not just hit/no-hit: box + strip plot of fold-change
values among *significant* (gene x experiment x timepoint) hits, per
group -- shows effect size, which 01_compute_hit_rates.py's rate-based
figure never did (it only used the categorical up/down/not-significant
call).

Only significant rows are plotted (not_significant log2fc values are
mostly near-zero noise and would wash out the comparison). Points are
colored by direction (blue = up, red = down, same validated diverging
pair as figure 01) since sign is the primary thing being compared
across groups.

Input:  data/03_log2fc_raw.csv
Output: figures/02_log2fc_distribution.png

Usage:
  uv run python analyses/2026-08-02-np_gene_expression_full_analysis/5_analyze/scripts/04_plot_log2fc_distribution.py
"""
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

STEP_DIR = Path(__file__).resolve().parents[1]
DATA_DIR = STEP_DIR / "data"
FIG_DIR = STEP_DIR / "figures"

BLUE = "#2a78d6"
RED = "#e34948"
INK = "#0b0b0b"
MUTED = "#898781"
GRID = "#e1e0d9"

GROUP_ORDER = ["H1_N", "positive_control_N", "background_N", "H1_P", "positive_control_P"]
LABELS = {
    "H1_N": "N target genes\n(N starvation)",
    "positive_control_N": "N positive\ncontrols",
    "background_N": "Background\n(N starvation)",
    "H1_P": "P target genes\n(P starvation)*",
    "positive_control_P": "P positive\ncontrols*",
}


def main() -> None:
    df = pd.read_csv(DATA_DIR / "03_log2fc_raw.csv")
    sig = df[df["expression_status"].isin(["significant_up", "significant_down"])].copy()
    sig["group"] = pd.Categorical(sig["group"], categories=GROUP_ORDER, ordered=True)
    sig["direction"] = sig["expression_status"].map({"significant_up": "Up", "significant_down": "Down"})

    fig, ax = plt.subplots(figsize=(11, 5.5), dpi=300)

    sns.boxplot(
        data=sig, x="group", y="log2fc", order=GROUP_ORDER, ax=ax,
        color="none", showfliers=False, width=0.5,
        boxprops={"edgecolor": MUTED, "facecolor": "none"},
        medianprops={"color": MUTED, "linewidth": 1.5},
        whiskerprops={"color": MUTED}, capprops={"color": MUTED},
    )
    sns.stripplot(
        data=sig, x="group", y="log2fc", order=GROUP_ORDER, hue="direction",
        palette={"Up": BLUE, "Down": RED}, ax=ax, alpha=0.6, size=4, jitter=0.25,
        hue_order=["Up", "Down"],
    )

    # 3 known outlier genes (log2fc 30-162, MED4 phosphate-starvation
    # microarray, PMM0707/PMM0708/PMM1416; PMM0708=phoA appears in both
    # H1_P and positive_control_P) -- flagged as likely artifacts by the
    # prior triage analysis, not filtered out there either. Clip the axis
    # so the bulk of the distribution stays readable; note in the title
    # rather than an in-plot annotation (avoids colliding with data/legend).
    n_clipped = int((sig["log2fc"] > 20).sum())
    ax.set_ylim(-10, 20)

    ax.axhline(0, color=MUTED, linewidth=1, linestyle="--")
    ax.set_xticks(range(len(GROUP_ORDER)))
    ax.set_xticklabels([LABELS[g] for g in GROUP_ORDER], fontsize=9, color=INK)
    ax.set_xlabel("")
    ax.set_ylabel("log2 fold change (significant hits only)", color=INK)
    ax.set_title(
        "Fold-change magnitude among significant hits, by group\n"
        "* phosphorus: no valid background to compare against (see notebook)  |  "
        f"{n_clipped} points off-scale (log2FC 30-162, MED4 phosphate microarray, known outliers)",
        fontsize=9, color=INK, loc="left",
    )
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.tick_params(axis="y", colors=MUTED)
    ax.grid(axis="y", color=GRID, linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)
    ax.legend(title="", loc="upper left", frameon=False)
    fig.tight_layout()

    out_path = FIG_DIR / "02_log2fc_distribution.png"
    fig.savefig(out_path, dpi=300)
    print(f"Wrote {out_path}")

    print("\nMedian |log2FC| per group (significant hits):")
    print(sig.groupby("group", observed=True)["log2fc"].apply(lambda s: s.abs().median()).reindex(GROUP_ORDER))


if __name__ == "__main__":
    main()
