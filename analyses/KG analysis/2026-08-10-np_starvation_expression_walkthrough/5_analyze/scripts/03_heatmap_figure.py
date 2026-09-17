"""
Step 5 -- main supplementary figure: experiment x gene heatmap.

Rows: all 10 in-scope experiments (grouped N then P), each at its single
chosen starvation timepoint (no pooling/collapsing across experiments or
timepoints -- see 5_analyze/notebook.md for the timepoint table and the
rationale for not merging independent studies). Row labels are the
Oxford-style citation "First author et al., YEAR - analysis type - strain";
a left-side bracket marks the nitrogen vs phosphorus block.
Columns: only genes that showed up- or down-regulation in >=1 experiment
(genes that were never significant anywhere, or never had data anywhere,
are dropped from this figure -- full 61-gene data stays in
data/01_target_gene_experiment_matrix.csv).
Cell color: 5 categories -- up, down, tested-not-significant, no data
(gene present in the strain but not reported in this experiment's table),
and gene-absent-from-strain (no locus tag in this strain's genome, from
the step-2 resolution -- drawn as a single diagonal line, not a colored
cell).

Style (2026-09-09, researcher-requested): Arial throughout, no bold,
larger fonts, darker "tested-not-significant" grey.

Inputs: data/01_target_gene_experiment_matrix.csv
Outputs: figures/01_gene_experiment_heatmap.png

Usage: uv run analyses/2026-08-10-np_starvation_expression_walkthrough/5_analyze/scripts/03_heatmap_figure.py
"""

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch

mpl.rcParams["font.family"] = "sans-serif"
mpl.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "DejaVu Sans"]
mpl.rcParams["font.weight"] = "normal"
mpl.rcParams["axes.titleweight"] = "normal"
mpl.rcParams["figure.titleweight"] = "normal"

BASE = Path(__file__).resolve().parent.parent
DATA_DIR = BASE / "data"
FIG_DIR = BASE / "figures"

# Diverging pair (dataviz palette.md): blue<->red, neutral gray midpoint.
COLOR_UP = "#e34948"       # red -- upregulated
COLOR_DOWN = "#2a78d6"     # blue -- downregulated
COLOR_MIXED = "#a8a69b"    # darker neutral gray -- tested, not significant
COLOR_NODATA = "#e7e4d8"   # distinctly darker beige than the near-white absent cell
COLOR_ABSENT = "#fdfcf9"   # near-white -- gene not in strain genome (3 diagonal strokes on top)
ABSENT_LINE = "#9a988e"    # diagonal strokes for absent cells
GRIDLINE = "#d7d5c9"       # visible cell borders
INK = "#0b0b0b"
MUTED = "#898781"
SURFACE = "#fcfcfb"

STATUS_TO_CODE = {
    "significant_up": 0,
    "significant_down": 1,
    "not_significant": 2,
    "no_data_at_timepoint": 3,
    "no_locus_in_strain": 4,
}
# Code 4 (gene absent from strain) is drawn as one diagonal line, not a fill;
# its colormap entry is just the surface color so the line reads cleanly.
CMAP = ListedColormap([COLOR_UP, COLOR_DOWN, COLOR_MIXED, COLOR_NODATA, COLOR_ABSENT])

EXPERIMENT_LABELS = {
    "10.1101/2025.11.24.690089_growth_state_pro99lown_nutrient_starvation_med4_proteomics_axenic": "Weissberg et al., 2025 - Proteomics - MED4",
    "10.1101/2025.11.24.690089_growth_state_pro99lown_nutrient_starvation_med4_rnaseq_axenic": "Weissberg et al., 2025 - RNA-seq - MED4",
    "10.1038/ismej.2017.88_nitrogen_stress_ndepleted_pro99_medium_med4_rnaseq": "Read et al., 2017 - RNA-seq - MED4",
    "10.1038/msb4100087_nitrogen_nitrogen_deprivation_med4_med4_microarray": "Tolonen et al., 2006 - Microarray - MED4",
    "10.1038/msb4100087_nitrogen_nitrogen_deprivation_mit9313_mit9313_microarray": "Tolonen et al., 2006 - Microarray - MIT9313",
    "10.1111/1462-2920.13104_phosphorus_plimited_natl2a_rnaseq_uninfected": "Lin et al., 2015 - RNA-seq - NATL2A",
    "10.1186/2046-9063-8-7_pi_limitation_mit9312_itraq": "Fuszard et al., 2012 - Proteomics - MIT9312",
    "10.1186/2046-9063-8-7_pi_limitation_natl2a_itraq": "Fuszard et al., 2012 - Proteomics - NATL2A",
    "10.1073/pnas.0601301103_phosphorus_phosphate_starvation_med4_microarray": "Martiny et al., 2006 - Microarray - MED4",
    "10.1073/pnas.0601301103_phosphorus_phosphate_starvation_mit9313_microarray": "Martiny et al., 2006 - Microarray - MIT9313*",
}
ROW_ORDER = list(EXPERIMENT_LABELS.keys())
ROW_NUTRIENT = ["N", "N", "N", "N", "N", "P", "P", "P", "P", "P"]


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

    cell = 0.46
    fig_width = cell * len(gene_order) + 9.2
    fig_height = cell * len(ROW_ORDER) + 3.6
    left_frac = 7.9 / fig_width
    right_frac = 1 - 0.35 / fig_width
    top_frac = 1 - 2.05 / fig_height
    bottom_frac = 1.55 / fig_height

    fig, ax = plt.subplots(figsize=(fig_width, fig_height), dpi=300)
    fig.patch.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)

    ax.imshow(pivot.values, cmap=CMAP, vmin=-0.5, vmax=4.5, aspect="auto")

    # Gene-absent-from-strain cells (code 4): three diagonal strokes (a sparse
    # hatch that reads clearly against the near-white fill), so they stay
    # distinct from the "No data" cells (code 3, pale beige) where the gene is
    # present but untested.
    values = pivot.values
    for (yi, xi), code in np.ndenumerate(values):
        if code == 4:
            x0, x1 = xi - 0.5, xi + 0.5
            yb, yt = yi + 0.5, yi - 0.5
            for (px0, py0), (px1, py1) in (((x0, yb), (x1, yt)),
                                           ((x0, yi), (xi, yt)),
                                           ((xi, yb), (x1, yi))):
                ax.plot([px0, px1], [py0, py1], color=ABSENT_LINE, linewidth=0.8, clip_on=True)

    ax.set_xticks(range(len(gene_order)))
    ax.set_xticklabels(gene_order, rotation=60, ha="left", fontsize=12, color=INK)
    ax.xaxis.tick_top()

    # Group labels above each gene block (N acquisition genes / P acquisition genes).
    n_gene_count_preview = sum(1 for g in gene_order if n_or_p_lookup[g] == "N")
    ax.text(n_gene_count_preview / 2 - 0.5, -2.9, "N acquisition genes", ha="center", va="bottom",
            fontsize=15, color="#2a78d6", clip_on=False)
    ax.text(n_gene_count_preview + (len(gene_order) - n_gene_count_preview) / 2 - 0.5, -2.9,
            "P acquisition genes", ha="center", va="bottom", fontsize=15,
            color="#e34948", clip_on=False)

    ax.set_yticks(range(len(ROW_ORDER)))
    ax.set_yticklabels([EXPERIMENT_LABELS[c] for c in ROW_ORDER], fontsize=15, color=INK)

    # Horizontal divider between nitrogen and phosphorus experiment blocks.
    n_count = ROW_NUTRIENT.index("P")
    ax.axhline(n_count - 0.5, color=INK, linewidth=1.4)
    # Vertical divider between N and P gene blocks.
    n_gene_count = sum(1 for g in gene_order if n_or_p_lookup[g] == "N")
    ax.axvline(n_gene_count - 0.5, color=INK, linewidth=1.4)

    ax.set_xticks(np.arange(-0.5, len(gene_order), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(ROW_ORDER), 1), minor=True)
    ax.grid(which="minor", color=GRIDLINE, linewidth=1.0)
    ax.tick_params(which="minor", length=0)
    ax.tick_params(which="major", length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)

    # left-side bracket: nitrogen block vs phosphorus block
    ytrans = ax.get_yaxis_transform()
    x_bar, x_txt = -0.25, -0.272
    for name, a, b in [("Nitrogen starvation", 0, n_count - 1),
                       ("Phosphorus starvation", n_count, len(ROW_ORDER) - 1)]:
        ax.plot([x_bar, x_bar], [a - 0.4, b + 0.4], color=INK, linewidth=1.4,
                transform=ytrans, clip_on=False)
        for y_end in (a - 0.4, b + 0.4):
            ax.plot([x_bar, x_bar + 0.008], [y_end, y_end], color=INK, linewidth=1.4,
                    transform=ytrans, clip_on=False)
        ax.text(x_txt, (a + b) / 2, name, rotation=90, ha="center", va="center",
                fontsize=15, color=INK, transform=ytrans, clip_on=False)

    legend_elements = [Patch(facecolor=COLOR_UP, label="Upregulated"),
                       Patch(facecolor=COLOR_DOWN, label="Downregulated"),
                       Patch(facecolor=COLOR_MIXED, label="Tested, not significant"),
                       Patch(facecolor=COLOR_NODATA, edgecolor=MUTED, label="No data (gene present, not in table)"),
                       Patch(facecolor=COLOR_ABSENT, edgecolor=ABSENT_LINE, hatch="///",
                             label="Gene absent from strain genome")]
    fig.legend(handles=legend_elements, loc="lower center",
               bbox_to_anchor=(0.5 * (left_frac + right_frac), 0.055),
               ncol=5, frameon=False, fontsize=13.5, handlelength=1.4, handleheight=1.4)

    fig.suptitle("Nitrogen- and phosphorus-acquisition gene response to nutrient starvation",
                 x=0.5 * (left_frac + right_frac), y=0.988, ha="center", va="top",
                 fontsize=19, color=INK)
    fig.text(0.5 * (left_frac + right_frac), 0.988 - 0.36 / fig_height,
             "genes with a response in >=1 experiment; each experiment at its single chosen starvation timepoint",
             ha="center", va="top", fontsize=12.5, color=INK)
    fig.text(0.5 * (left_frac + right_frac), 0.012,
             "* MIT9313 phosphorus (Martiny): significance-criterion / timepoint mismatch (24h used, not 48h). "
             "Martiny 'tested, not significant' cells include genes absent from that publication's filtered table. "
             "See 5_analyze/notebook.md.",
             ha="center", va="bottom", fontsize=10.5, color=MUTED)

    fig.subplots_adjust(left=left_frac, right=right_frac, top=top_frac, bottom=bottom_frac)
    out_path = FIG_DIR / "01_gene_experiment_heatmap.png"
    fig.savefig(out_path, dpi=300, facecolor=SURFACE)
    print(f"Wrote {out_path}")
    print(f"{len(ROW_ORDER)} experiments x {len(gene_order)} genes")


if __name__ == "__main__":
    main()
