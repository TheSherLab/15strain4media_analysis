"""
Step 5 -- main supplementary figure: experiment x gene heatmap.

Rows: all 10 in-scope experiments (grouped N then P), each at its single
chosen starvation timepoint (no pooling/collapsing across experiments or
timepoints -- see 5_analyze/notebook.md for the timepoint table and the
rationale for not merging independent studies).
Columns: only genes that showed up- or down-regulation in >=1 experiment
(genes that were never significant anywhere, or never had data anywhere,
are dropped from this figure -- full 61-gene data stays in
data/01_target_gene_experiment_matrix.csv).
Cell color: 4 categories -- up, down, tested-not-significant, no data.

Inputs: data/01_target_gene_experiment_matrix.csv
Outputs: figures/01_gene_experiment_heatmap.png

Usage: uv run analyses/2026-08-10-np_starvation_expression_walkthrough/5_analyze/scripts/03_heatmap_figure.py
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch

BASE = Path(__file__).resolve().parent.parent
DATA_DIR = BASE / "data"
FIG_DIR = BASE / "figures"

# Diverging pair (dataviz palette.md): blue<->red, neutral gray midpoint.
COLOR_UP = "#e34948"       # red -- upregulated
COLOR_DOWN = "#2a78d6"     # blue -- downregulated
COLOR_MIXED = "#c3c2b7"    # neutral gray -- tested, not significant
COLOR_NODATA = "#f4f3ef"   # near-surface, distinct from gray -- no data
INK = "#0b0b0b"
SECONDARY_INK = "#52514e"
MUTED = "#898781"
SURFACE = "#fcfcfb"

STATUS_TO_CODE = {
    "significant_up": 0,
    "significant_down": 1,
    "not_significant": 2,
    "no_locus_in_strain": 3,
    "no_data_at_timepoint": 3,
}
CMAP = ListedColormap([COLOR_UP, COLOR_DOWN, COLOR_MIXED, COLOR_NODATA])

EXPERIMENT_LABELS = {
    "10.1101/2025.11.24.690089_growth_state_pro99lown_nutrient_starvation_med4_proteomics_axenic": "N: Weissberg -- Proteomics (MED4)",
    "10.1101/2025.11.24.690089_growth_state_pro99lown_nutrient_starvation_med4_rnaseq_axenic": "N: Weissberg -- RNA-seq (MED4)",
    "10.1038/ismej.2017.88_nitrogen_stress_ndepleted_pro99_medium_med4_rnaseq": "N: Read -- RNA-seq (MED4)",
    "10.1038/msb4100087_nitrogen_nitrogen_deprivation_med4_med4_microarray": "N: Tolonen -- Microarray (MED4)",
    "10.1038/msb4100087_nitrogen_nitrogen_deprivation_mit9313_mit9313_microarray": "N: Tolonen -- Microarray (MIT9313)",
    "10.1111/1462-2920.13104_phosphorus_plimited_natl2a_rnaseq_uninfected": "P: Lin -- RNA-seq (NATL2A)",
    "10.1186/2046-9063-8-7_pi_limitation_mit9312_itraq": "P: Fuszard -- Proteomics (MIT9312)",
    "10.1186/2046-9063-8-7_pi_limitation_natl2a_itraq": "P: Fuszard -- Proteomics (NATL2A)",
    "10.1073/pnas.0601301103_phosphorus_phosphate_starvation_med4_microarray": "P: Martiny -- Microarray (MED4)",
    "10.1073/pnas.0601301103_phosphorus_phosphate_starvation_mit9313_microarray": "P: Martiny -- Microarray (MIT9313)*",
}
ROW_ORDER = list(EXPERIMENT_LABELS.keys())


def main() -> None:
    matrix = pd.read_csv(DATA_DIR / "01_target_gene_experiment_matrix.csv")
    matrix["code"] = matrix["status"].map(STATUS_TO_CODE)

    has_signal = matrix.groupby("gene_name")["code"].apply(lambda s: s.isin([0, 1]).any())
    kept_genes = set(has_signal[has_signal].index)
    dropped = matrix["gene_name"].nunique() - len(kept_genes)
    print(f"Keeping {len(kept_genes)} of {matrix['gene_name'].nunique()} genes "
          f"(dropped {dropped} genes with no up/down signal in any experiment)")

    n_or_p_lookup = matrix.drop_duplicates("gene_name").set_index("gene_name")["n_or_p"]
    gene_order = sorted(kept_genes, key=lambda g: (n_or_p_lookup[g] != "N", g))  # N block first, alpha within

    pivot = matrix.pivot_table(index="experiment_id", columns="gene_name", values="code", aggfunc="first")
    pivot = pivot.reindex(index=ROW_ORDER, columns=gene_order)

    col_colors = ["#2a78d6" if n_or_p_lookup[g] == "N" else "#e34948" for g in gene_order]

    cell = 0.42
    fig_width = cell * len(gene_order) + 3.2
    fig_height = cell * len(ROW_ORDER) + 3.6
    top_frac = 1 - 1.95 / fig_height
    bottom_frac = 0.85 / fig_height
    left_frac = 3.7 / fig_width
    right_frac = 1 - 0.35 / fig_width

    fig, ax = plt.subplots(figsize=(fig_width, fig_height), dpi=300)
    fig.patch.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)

    ax.imshow(pivot.values, cmap=CMAP, vmin=-0.5, vmax=3.5, aspect="auto")

    ax.set_xticks(range(len(gene_order)))
    ax.set_xticklabels(gene_order, rotation=60, ha="left", fontsize=11, fontfamily="monospace")
    for tick, color in zip(ax.get_xticklabels(), col_colors):
        tick.set_color(color)
    ax.xaxis.tick_top()

    # Group labels above each gene block (N acquisition genes / P acquisition genes).
    n_gene_count_preview = sum(1 for g in gene_order if n_or_p_lookup[g] == "N")
    ax.text(n_gene_count_preview / 2 - 0.5, -3.0, "N acquisition genes", ha="center", va="bottom",
             fontsize=12.5, fontweight="bold", color="#2a78d6", clip_on=False)
    ax.text(n_gene_count_preview + (len(gene_order) - n_gene_count_preview) / 2 - 0.5, -3.0,
             "P acquisition genes", ha="center", va="bottom", fontsize=12.5, fontweight="bold",
             color="#e34948", clip_on=False)

    ax.set_yticks(range(len(ROW_ORDER)))
    ax.set_yticklabels([EXPERIMENT_LABELS[c] for c in ROW_ORDER], fontsize=11.5, color=SECONDARY_INK)

    # Horizontal divider between nitrogen and phosphorus experiment blocks.
    n_count = sum(1 for c in ROW_ORDER if EXPERIMENT_LABELS[c].startswith("N:"))
    ax.axhline(n_count - 0.5, color=INK, linewidth=1.4)
    # Vertical divider between N and P gene blocks.
    n_gene_count = sum(1 for g in gene_order if n_or_p_lookup[g] == "N")
    ax.axvline(n_gene_count - 0.5, color=INK, linewidth=1.4)

    ax.set_xticks(np.arange(-0.5, len(gene_order), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(ROW_ORDER), 1), minor=True)
    ax.grid(which="minor", color=SURFACE, linewidth=1.8)
    ax.tick_params(which="minor", length=0)
    ax.tick_params(which="major", length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)

    legend_elements = [Patch(facecolor=COLOR_UP, label="Upregulated"),
                        Patch(facecolor=COLOR_DOWN, label="Downregulated"),
                        Patch(facecolor=COLOR_MIXED, label="Tested, not significant"),
                        Patch(facecolor=COLOR_NODATA, edgecolor=MUTED, label="No data")]
    fig.legend(handles=legend_elements, loc="lower center", bbox_to_anchor=(0.5 * (left_frac + right_frac), 0.08),
               ncol=4, frameon=False, fontsize=11, handlelength=1.4, handleheight=1.4)

    fig.suptitle("Nitrogen- and phosphorus-acquisition gene response to nutrient starvation",
                  x=0.5 * (left_frac + right_frac), y=0.995, ha="center", va="top",
                  fontsize=14, fontweight="bold", color=INK)
    fig.text(0.5 * (left_frac + right_frac), 0.995 - 0.3 / fig_height,
              "genes with a response in >=1 experiment; each experiment at its single chosen starvation timepoint",
              ha="center", va="top", fontsize=10.5, color=SECONDARY_INK)
    fig.text(0.5 * (left_frac + right_frac), 0.20 / fig_height,
              "* MIT9313 phosphorus (Martiny): significance criterion mismatch flagged -- see notebook.md",
              ha="center", va="bottom", fontsize=9, color=MUTED, style="italic")

    fig.subplots_adjust(left=left_frac, right=right_frac, top=top_frac, bottom=bottom_frac)
    out_path = FIG_DIR / "01_gene_experiment_heatmap.png"
    fig.savefig(out_path, dpi=300, facecolor=SURFACE)
    print(f"Wrote {out_path}")
    print(f"{len(ROW_ORDER)} experiments x {len(gene_order)} genes")


if __name__ == "__main__":
    main()
