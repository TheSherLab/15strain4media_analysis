"""
Step 5 -- percentage heatmap, collapsed by experiment nutrient (not by
individual experiment): 2 rows -- "N starvation experiments" and
"P starvation experiments" -- one column per gene (all genes, both N-
and P-acquisition annotated). Each cell = that gene's % of tests
significant under that row's nutrient, colored red (up) / blue (down) /
gray (tested, not significant) / hatched (no data), intensity scaled by
the percentage. Directly visualizes matched (H1, diagonal blocks) vs.
cross (H2, off-diagonal blocks) response for every gene at once.

Inputs: data/01_target_gene_experiment_matrix.csv
Outputs: data/06_pct_by_nutrient.csv
         figures/04_pct_heatmap_by_nutrient.png

Usage: uv run analyses/2026-08-10-np_starvation_expression_walkthrough/5_analyze/scripts/06_pct_heatmap_by_nutrient.py
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap

BASE = Path(__file__).resolve().parent.parent
DATA_DIR = BASE / "data"
FIG_DIR = BASE / "figures"

TESTED_STATUSES = ["significant_up", "significant_down", "not_significant"]

INK = "#0b0b0b"
SECONDARY_INK = "#52514e"
MUTED = "#898781"
SURFACE = "#fcfcfb"
NA_COLOR = "#efeee9"
COLOR_GRAY = "#c3c2b7"

RED_RAMP = LinearSegmentedColormap.from_list("red_ramp", ["#fbe4e3", "#e34948", "#7a1d1c"])
BLUE_RAMP = LinearSegmentedColormap.from_list("blue_ramp", ["#dce9fa", "#2a78d6", "#123a66"])


def main() -> None:
    matrix = pd.read_csv(DATA_DIR / "01_target_gene_experiment_matrix.csv")
    tested = matrix[matrix["status"].isin(TESTED_STATUSES)]

    summary = tested.groupby(["gene_name", "n_or_p", "nutrient"]).agg(
        n_tested=("status", "size"),
        n_up=("status", lambda s: (s == "significant_up").sum()),
        n_down=("status", lambda s: (s == "significant_down").sum()),
    ).reset_index()
    summary["pct_up"] = 100 * summary["n_up"] / summary["n_tested"]
    summary["pct_down"] = 100 * summary["n_down"] / summary["n_tested"]
    summary.to_csv(DATA_DIR / "06_pct_by_nutrient.csv", index=False)
    print(f"Wrote {len(summary)} rows to data/06_pct_by_nutrient.csv")

    all_genes = matrix.drop_duplicates("gene_name").set_index("gene_name")["n_or_p"]
    n_genes = sorted(all_genes[all_genes == "N"].index)
    p_genes = sorted(all_genes[all_genes == "P"].index)
    gene_order = n_genes + p_genes
    n_gene_count = len(n_genes)

    # value: signed magnitude (positive = %up dominant, negative = %down dominant, 0 = not-sig dominant, NaN = no data)
    lookup = summary.set_index(["gene_name", "nutrient"])
    rows = ["nitrogen", "phosphorus"]
    row_labels = ["N starvation\nexperiments", "P starvation\nexperiments"]
    value = np.full((2, len(gene_order)), np.nan)
    for i, nutrient in enumerate(rows):
        for j, g in enumerate(gene_order):
            if (g, nutrient) not in lookup.index:
                continue
            r = lookup.loc[(g, nutrient)]
            if r["n_up"] >= r["n_down"] and r["n_up"] > 0:
                value[i, j] = r["pct_up"]
            elif r["n_down"] > r["n_up"]:
                value[i, j] = -r["pct_down"]
            else:
                value[i, j] = 0  # tested, not significant

    cell = 0.30
    fig_width = cell * len(gene_order) + 2.6
    fig_height = 5.0
    top_frac = 1 - 2.7 / fig_height
    bottom_frac = 0.75 / fig_height
    left_frac = 2.15 / fig_width
    right_frac = 1 - 0.3 / fig_width

    fig, ax = plt.subplots(figsize=(fig_width, fig_height), dpi=300)
    fig.patch.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)

    for i in range(2):
        for j in range(len(gene_order)):
            v = value[i, j]
            x0, y0 = j - 0.5, i - 0.5
            if np.isnan(v):
                ax.add_patch(plt.Rectangle((x0, y0), 1, 1, facecolor=NA_COLOR, hatch="////",
                                             edgecolor="#d5d3cb", linewidth=0))
            elif v == 0:
                ax.add_patch(plt.Rectangle((x0, y0), 1, 1, facecolor=COLOR_GRAY, linewidth=0))
            elif v > 0:
                color = RED_RAMP(0.15 + 0.85 * v / 100)
                ax.add_patch(plt.Rectangle((x0, y0), 1, 1, facecolor=color, linewidth=0))
                if v >= 55:
                    ax.text(j, i, f"{v:.0f}", ha="center", va="center", fontsize=7, color="white", fontweight="bold")
            else:
                color = BLUE_RAMP(0.15 + 0.85 * (-v) / 100)
                ax.add_patch(plt.Rectangle((x0, y0), 1, 1, facecolor=color, linewidth=0))
                if -v >= 55:
                    ax.text(j, i, f"{-v:.0f}", ha="center", va="center", fontsize=7, color="white", fontweight="bold")

    ax.set_xlim(-0.5, len(gene_order) - 0.5)
    ax.set_ylim(1.5, -0.5)

    ax.set_xticks(range(len(gene_order)))
    ax.set_xticklabels(gene_order, rotation=60, ha="left", fontsize=9, fontfamily="monospace")
    col_colors = ["#2a78d6" if all_genes[g] == "N" else "#e34948" for g in gene_order]
    for tick, color in zip(ax.get_xticklabels(), col_colors):
        tick.set_color(color)
    ax.xaxis.tick_top()

    ax.set_yticks([0, 1])
    ax.set_yticklabels(row_labels, fontsize=11, fontweight="bold", color=INK)

    ax.axvline(n_gene_count - 0.5, color=INK, linewidth=1.4)
    ax.axhline(0.5, color=INK, linewidth=1.4)

    ax.set_xticks(np.arange(-0.5, len(gene_order), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, 2, 1), minor=True)
    ax.grid(which="minor", color=SURFACE, linewidth=1.8)
    ax.tick_params(which="minor", length=0)
    ax.tick_params(which="major", length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)

    # Figure-fraction x for each block's center, computed from the axes' actual
    # left/right edges -- more robust than guessing a data-coordinate y offset
    # (which is sensitive to how much vertical space the rotated column labels take).
    n_gene_frac = left_frac + (n_gene_count / 2) / len(gene_order) * (right_frac - left_frac)
    p_gene_frac = left_frac + (n_gene_count + (len(gene_order) - n_gene_count) / 2) / len(gene_order) * (right_frac - left_frac)
    group_label_y = 0.995 - 0.62 / fig_height
    fig.text(n_gene_frac, group_label_y, "N acquisition genes", ha="center", va="top",
              fontsize=11.5, fontweight="bold", color="#2a78d6")
    fig.text(p_gene_frac, group_label_y, "P acquisition genes", ha="center", va="top",
              fontsize=11.5, fontweight="bold", color="#e34948")

    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor=RED_RAMP(0.9), label="Upregulated (darker = higher %)"),
        Patch(facecolor=BLUE_RAMP(0.9), label="Downregulated (darker = higher %)"),
        Patch(facecolor=COLOR_GRAY, label="Tested, not significant"),
        Patch(facecolor=NA_COLOR, edgecolor=MUTED, hatch="////", label="No data"),
    ]
    fig.legend(handles=legend_elements, loc="lower center",
               bbox_to_anchor=(0.5 * (left_frac + right_frac), 0.06), ncol=2, frameon=False,
               fontsize=8.5, handlelength=1.3, handleheight=1.3)

    fig.suptitle("Gene response, by which nutrient was starved",
                  x=left_frac, y=0.98, ha="left", fontsize=13, fontweight="bold", color=INK)
    fig.text(left_frac, 0.98 - 0.26 / fig_height,
              "cell = % of that gene's tests significant under N (row 1) or P (row 2) starvation; color intensity = the percentage",
              ha="left", va="top", fontsize=8.5, color=SECONDARY_INK)

    fig.subplots_adjust(left=left_frac, right=right_frac, top=top_frac, bottom=bottom_frac + 0.09)
    out_path = FIG_DIR / "04_pct_heatmap_by_nutrient.png"
    fig.savefig(out_path, dpi=300, facecolor=SURFACE)
    print(f"Wrote {out_path}")
    print(f"{len(gene_order)} genes ({n_gene_count} N, {len(gene_order)-n_gene_count} P)")


if __name__ == "__main__":
    main()
