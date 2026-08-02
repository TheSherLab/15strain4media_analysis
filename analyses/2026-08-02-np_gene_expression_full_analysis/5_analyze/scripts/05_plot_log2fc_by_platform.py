"""Log2FC distribution split by omics platform (RNA-seq, proteomics,
microarray) instead of pooled -- this methodology's statistical-rigor
rule says platform log2FC magnitudes are not directly comparable
(different dynamic ranges, detection limits, noise profiles), so a
figure that pools all three would be misleading. One panel per
platform, same 5 groups and colors as 02_log2fc_distribution.png.

Input:  data/03_log2fc_raw.csv
Output: figures/03_log2fc_by_platform.png

Usage:
  uv run python analyses/2026-08-02-np_gene_expression_full_analysis/5_analyze/scripts/05_plot_log2fc_by_platform.py
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
    "H1_N": "N target",
    "positive_control_N": "N pos.\nctrl",
    "background_N": "Back-\nground",
    "H1_P": "P target*",
    "positive_control_P": "P pos.\nctrl*",
}
PLATFORM_ORDER = ["RNASEQ", "PROTEOMICS", "MICROARRAY"]


def main() -> None:
    df = pd.read_csv(DATA_DIR / "03_log2fc_raw.csv")
    sig = df[df["expression_status"].isin(["significant_up", "significant_down"])].copy()
    sig["group"] = pd.Categorical(sig["group"], categories=GROUP_ORDER, ordered=True)
    sig["direction"] = sig["expression_status"].map({"significant_up": "Up", "significant_down": "Down"})

    n_clipped_total = int((sig["log2fc"] > 20).sum())
    sig_plot = sig[sig["log2fc"] <= 20]

    fig, axes = plt.subplots(1, 3, figsize=(13, 5), dpi=300, sharey=True)

    for ax, platform in zip(axes, PLATFORM_ORDER):
        sub = sig_plot[sig_plot["omics_type"] == platform]
        n_per_group = sub.groupby("group", observed=True).size().reindex(GROUP_ORDER, fill_value=0)

        if not sub.empty:
            sns.boxplot(
                data=sub, x="group", y="log2fc", order=GROUP_ORDER, ax=ax,
                showfliers=False, width=0.5,
                boxprops={"edgecolor": MUTED, "facecolor": "none"},
                medianprops={"color": MUTED, "linewidth": 1.5},
                whiskerprops={"color": MUTED}, capprops={"color": MUTED},
            )
            sns.stripplot(
                data=sub, x="group", y="log2fc", order=GROUP_ORDER, hue="direction",
                palette={"Up": BLUE, "Down": RED}, ax=ax, alpha=0.6, size=4, jitter=0.25,
                hue_order=["Up", "Down"], legend=(platform == PLATFORM_ORDER[0]),
            )

        ax.axhline(0, color=MUTED, linewidth=1, linestyle="--")
        ax.set_xticks(range(len(GROUP_ORDER)))
        ax.set_xticklabels([LABELS[g] for g in GROUP_ORDER], fontsize=8, color=INK)
        ax.set_xlabel("")
        ax.set_title(f"{platform}  (n={len(sub)})", fontsize=10, color=INK)
        for spine in ("top", "right"):
            ax.spines[spine].set_visible(False)
        ax.grid(axis="y", color=GRID, linewidth=0.8, zorder=0)
        ax.set_axisbelow(True)
        if ax is axes[0]:
            ax.set_ylabel("log2 fold change (significant hits only)", color=INK)
            ax.legend(title="", loc="upper left", frameon=False, fontsize=8)
        else:
            ax.set_ylabel("")

    axes[0].set_ylim(-10, 20)
    fig.suptitle(
        "Fold-change magnitude by platform (not pooled -- magnitudes are not comparable across platforms)\n"
        f"* phosphorus: no valid background  |  {n_clipped_total} points off-scale "
        "(log2FC 30-162, MED4 phosphate microarray, known outliers)",
        fontsize=9.5, color=INK, x=0.01, y=0.98, ha="left",
    )
    fig.tight_layout(rect=(0, 0, 1, 0.86))

    out_path = FIG_DIR / "03_log2fc_by_platform.png"
    fig.savefig(out_path, dpi=300)
    print(f"Wrote {out_path}")

    print("\nn (significant hits) per group x platform:")
    print(sig.groupby(["omics_type", "group"], observed=True).size().unstack(fill_value=0).reindex(columns=GROUP_ORDER))


if __name__ == "__main__":
    main()
