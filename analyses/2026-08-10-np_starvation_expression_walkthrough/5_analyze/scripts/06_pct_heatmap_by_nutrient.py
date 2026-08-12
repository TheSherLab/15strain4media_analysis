"""
Step 5 -- percentage heatmap, collapsed by experiment nutrient, split into
two stacked panels by the gene's own N/P annotation: top panel = N-acquisition
genes, bottom panel = P-acquisition genes. Each panel has its own pair of
rows -- "N starvation experiments" / "P starvation experiments" -- so the
matched response (H1, top row of each panel) and cross response (H2, bottom
row of each panel) are each fully visible within one gene group instead of
sharing a row with the other group.

Each cell = that gene's % of tests significant under that row's nutrient,
colored red (up) / blue (down) / gray (tested, not significant) / hatched
(no data), intensity scaled by the percentage.

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
from matplotlib.patches import Patch

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

ROWS = ["nitrogen", "phosphorus"]
ROW_LABELS = ["N starvation\nexperiments", "P starvation\nexperiments"]


def panel_values(summary: pd.DataFrame, genes: list[str]) -> np.ndarray:
    """2 x len(genes) matrix: signed magnitude (positive = %up dominant,
    negative = %down dominant, 0 = not-sig dominant, NaN = no data)."""
    lookup = summary.set_index(["gene_name", "nutrient"])
    value = np.full((2, len(genes)), np.nan)
    for i, nutrient in enumerate(ROWS):
        for j, g in enumerate(genes):
            if (g, nutrient) not in lookup.index:
                continue
            r = lookup.loc[(g, nutrient)]
            if r["n_up"] >= r["n_down"] and r["n_up"] > 0:
                value[i, j] = r["pct_up"]
            elif r["n_down"] > r["n_up"]:
                value[i, j] = -r["pct_down"]
            else:
                value[i, j] = 0  # tested, not significant
    return value


def draw_panel(ax, value: np.ndarray, genes: list[str], show_col_labels: bool) -> None:
    for i in range(2):
        for j in range(len(genes)):
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

    ax.set_xlim(-0.5, len(genes) - 0.5)
    ax.set_ylim(1.5, -0.5)

    ax.set_xticks(range(len(genes)))
    if show_col_labels:
        ax.set_xticklabels(genes, rotation=60, ha="left", fontsize=9, fontfamily="monospace", color=INK)
        ax.xaxis.tick_top()
    else:
        ax.set_xticklabels([])

    ax.set_yticks([0, 1])
    ax.set_yticklabels(ROW_LABELS, fontsize=10.5, fontweight="bold", color=INK)
    ax.text(0, -0.62, "Gene", transform=ax.get_yaxis_transform(), ha="right", va="bottom",
            fontsize=9, fontweight="bold", style="italic", color=SECONDARY_INK, clip_on=False)

    ax.axhline(0.5, color=INK, linewidth=1.4)

    ax.set_xticks(np.arange(-0.5, len(genes), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, 2, 1), minor=True)
    ax.grid(which="minor", color=SURFACE, linewidth=1.8)
    ax.tick_params(which="minor", length=0)
    ax.tick_params(which="major", length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)


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

    value_n = panel_values(summary, n_genes)
    value_p = panel_values(summary, p_genes)

    cell = 0.30
    fig_width = cell * max(len(n_genes), len(p_genes)) + 2.6
    left_frac = 2.15 / fig_width
    right_frac = 1 - 0.3 / fig_width

    # Vertical geometry, in inches, laid out top-to-bottom. Each panel needs
    # headroom above its own axis for the rotated (60-degree) gene-name column
    # labels, which matplotlib draws *outside* the axes' bounding box.
    TITLE_H = 0.55
    GROUPLABEL_H = 0.28
    COLLABEL_H = 0.85
    ROWS_H = 0.85
    PANEL_GAP = 0.20
    LEGEND_H = 0.65

    fig_height = TITLE_H + 2 * (GROUPLABEL_H + COLLABEL_H + ROWS_H) + PANEL_GAP + LEGEND_H

    cursor = 0.0  # inches consumed from the top

    def take(h: float) -> tuple[float, float]:
        nonlocal cursor
        top = 1 - cursor / fig_height
        cursor += h
        bottom = 1 - cursor / fig_height
        return top, bottom

    _title_top, _title_bottom = take(TITLE_H)
    n_group_top, _ = take(GROUPLABEL_H)
    _, n_collabel_bottom = take(COLLABEL_H)
    n_ax_top, n_ax_bottom = take(ROWS_H)
    take(PANEL_GAP)
    p_group_top, _ = take(GROUPLABEL_H)
    _, p_collabel_bottom = take(COLLABEL_H)
    p_ax_top, p_ax_bottom = take(ROWS_H)

    fig = plt.figure(figsize=(fig_width, fig_height), dpi=300)
    fig.patch.set_facecolor(SURFACE)

    ax_n = fig.add_axes((left_frac, n_ax_bottom, right_frac - left_frac, n_ax_top - n_ax_bottom))
    ax_p = fig.add_axes((left_frac, p_ax_bottom, right_frac - left_frac, p_ax_top - p_ax_bottom))

    draw_panel(ax_n, value_n, n_genes, show_col_labels=True)
    draw_panel(ax_p, value_p, p_genes, show_col_labels=True)

    fig.text(left_frac, n_group_top, "N acquisition genes",
              ha="left", va="top", fontsize=12, fontweight="bold", color=INK)
    fig.text(left_frac, p_group_top, "P acquisition genes",
              ha="left", va="top", fontsize=12, fontweight="bold", color=INK)

    legend_elements = [
        Patch(facecolor=RED_RAMP(0.9), label="Upregulated (darker = higher %)"),
        Patch(facecolor=BLUE_RAMP(0.9), label="Downregulated (darker = higher %)"),
        Patch(facecolor=COLOR_GRAY, label="Tested, not significant"),
        Patch(facecolor=NA_COLOR, edgecolor=MUTED, hatch="////", label="No data"),
    ]
    fig.legend(handles=legend_elements, loc="center",
               bbox_to_anchor=(0.5 * (left_frac + right_frac), p_ax_bottom / 2), ncol=2, frameon=False,
               fontsize=8.5, handlelength=1.3, handleheight=1.3)

    fig.suptitle("Gene response, by which nutrient was starved",
                  x=left_frac, y=_title_top, ha="left", va="top", fontsize=13, fontweight="bold", color=INK)
    fig.text(left_frac, _title_top - 0.28 / fig_height,
              "cell = % of that gene's tests significant under N or P starvation; color intensity = the percentage. "
              "N-acquisition genes (top) and P-acquisition genes (bottom) shown as separate panels.",
              ha="left", va="top", fontsize=8.5, color=SECONDARY_INK)

    out_path = FIG_DIR / "04_pct_heatmap_by_nutrient.png"
    fig.savefig(out_path, dpi=300, facecolor=SURFACE)
    print(f"Wrote {out_path}")
    print(f"{len(n_genes)} N genes, {len(p_genes)} P genes")


if __name__ == "__main__":
    main()
