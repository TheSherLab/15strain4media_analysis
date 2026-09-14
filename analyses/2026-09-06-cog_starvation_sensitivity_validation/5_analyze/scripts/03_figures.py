"""
Step 5 (redo 2026-09-09) -- two figures.

1. Experiment x COG heatmap
   - rows: the 10 in-scope experiments (5 nitrogen block, then 5 phosphorus
     block), each at its single chosen starvation timepoint. Row labels are
     the Oxford-style citation "First author et al., YEAR - analysis type -
     strain"; a left-side bracket marks the nitrogen vs phosphorus block.
   - columns: the 41 COGs in the researcher's own functional ordering
     (`00_source_normixed_cogs.csv` sheet order == the ordered screenshot),
     bracketed by her "General annotation" functional category.
   - a COG resolves to several genes per strain; one cell per COG x
     experiment. Cell text k/n/m = significant / measured here / genome
     copies; colour = dominant direction, intensity ~ k/n. Cells whose
     genes split direction are drawn in orange and show "up/down" for k.

2. Bar chart: % of tests significant -- background line vs nitrogen pooled
   / N / mixed and phosphorus pooled / N / mixed.

Figure style (2026-09-09, researcher-requested): Arial throughout, no bold,
larger fonts, darker "tested-not-significant" grey, single diagonal line
for "absent from genome" cells.

Inputs: data/01_target_gene_experiment_matrix.csv, data/02_hit_rate_results.csv,
        ../../2_kg_selection/data/00_source_normixed_cogs.csv,
        ../../3_analysis_framing/data/02_negative_background_pool.csv
Outputs: figures/01_gene_experiment_heatmap.png,
         figures/02_pct_significant_vs_background.png

Usage: uv run analyses/2026-09-06-cog_starvation_sensitivity_validation/5_analyze/scripts/03_figures.py
"""

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch, Rectangle

mpl.rcParams["font.family"] = "sans-serif"
mpl.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "DejaVu Sans"]
mpl.rcParams["font.weight"] = "normal"
mpl.rcParams["axes.titleweight"] = "normal"
mpl.rcParams["figure.titleweight"] = "normal"

BASE = Path(__file__).resolve().parent.parent
DATA_DIR = BASE / "data"
FIG_DIR = BASE / "figures"
STEP2_DATA = BASE.parent / "2_kg_selection" / "data"
FRAMING_DATA = BASE.parent / "3_analysis_framing" / "data"

# walkthrough palette (5_analyze/scripts/03_heatmap_figure.py there)
COLOR_UP = "#e34948"
COLOR_DOWN = "#2a78d6"
COLOR_BOTH = "#e8913a"       # orange -- a COG's genes split up/down in one cell
COLOR_NOTSIG = "#a8a69b"     # darker grey than the walkthrough (#c3c2b7)
COLOR_NODATA = "#e7e4d8"     # distinctly darker than the near-white absent cell
COLOR_ABSENT = "#fdfcf9"     # near-white; 3 diagonal strokes drawn on top
HATCH_INK = "#9a988e"        # diagonal strokes for absent cells
GRIDLINE = "#d7d5c9"         # visible cell borders
INK = "#0b0b0b"
MUTED = "#898781"
SURFACE = "#fcfcfb"

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
ROW_ORDER = list(EXPERIMENT_LABELS)
ROW_NUTRIENT = ["N", "N", "N", "N", "N", "P", "P", "P", "P", "P"]


def col_label(row) -> str:
    name = row["Protien"]
    name = "" if pd.isna(name) or not str(name).strip() else str(name).strip()
    name = name.split(" (")[0].split(",")[0].strip()  # "prop (including AmpG)" -> "prop"
    return f"{name or row['COGs']}  {row['COGs']}"


# tidy the researcher's "General annotation" category names for display
CATEGORY_DISPLAY = {
    "Quality contol- DNA level": "Quality control\nDNA level",
    "Quality contol- Protein level": "Quality control\nProtein level",
    "Amino acid and mixotrophy related biosynthesis": "Amino acid & mixotrophy\nbiosynthesis",
    "LPS": "LPS",
    "Membrane linker": "Membrane\nlinker",
    "Exopolysaccharide": "Exopoly-\nsaccharide",
    "Biofilm/Attachment": "Biofilm /\nAttachment",
    "Energy production": "Energy\nproduction",
    "Translation, transcription and signal transduction": "Translation, transcription\n& signal transduction",
    "RTX toxins": "RTX\ntoxins",
    "mixeds": "mixed",
}

TESTED = {"significant_up", "significant_down", "not_significant"}


def make_heatmap(matrix, cogs):
    order = cogs.sort_values("order").reset_index(drop=True)
    order["General annotation"] = order["General annotation"].ffill()
    cog_order = order["COGs"].tolist()
    labels = order.apply(col_label, axis=1).tolist()

    # contiguous spans of one "General annotation" category, in column order
    groups = []
    start = 0
    cats = order["General annotation"].tolist()
    for i in range(1, len(cats) + 1):
        if i == len(cats) or cats[i] != cats[start]:
            groups.append((start, i - 1, cats[start]))
            start = i

    # One cell per (COG, experiment). A COG has several genes per strain, so
    # the cell summarises them: fraction of the COG's TESTED genes that are
    # significant, and which direction dominates. Colour = direction with more
    # hits (orange if the genes disagree); intensity = that fraction; grey if
    # no gene is significant; pale if the COG has genes present but none in
    # this experiment's table; single diagonal line if absent from the strain.
    n_cog, n_exp = len(cog_order), len(ROW_ORDER)
    rgba = np.zeros((n_exp, n_cog, 4))
    absent_mask = np.zeros((n_exp, n_cog), dtype=bool)
    frac_txt = {}
    ci = {c: i for i, c in enumerate(cog_order)}
    ri = {e: i for i, e in enumerate(ROW_ORDER)}

    from matplotlib.colors import to_rgb
    up_rgb, down_rgb = np.array(to_rgb(COLOR_UP)), np.array(to_rgb(COLOR_DOWN))
    both_rgb = np.array(to_rgb(COLOR_BOTH))
    notsig_rgb, nodata_rgb, absent_rgb = (np.array(to_rgb(c)) for c in (COLOR_NOTSIG, COLOR_NODATA, COLOR_ABSENT))

    for (cog, exp), g in matrix.groupby(["cog_number", "experiment_id"]):
        if cog not in ci or exp not in ri:
            continue
        yi, xi = ri[exp], ci[cog]
        st = g["status"].value_counts().to_dict()
        n_tested = sum(st.get(s, 0) for s in TESTED)
        n_genome = int(g["locus_tag"].notna().sum())  # COG copies in this strain
        if n_tested == 0:
            if st.get("no_data_at_timepoint", 0) > 0:
                rgba[yi, xi, :3], rgba[yi, xi, 3] = nodata_rgb, 1.0
            else:  # only no_locus_in_strain
                rgba[yi, xi, :3], rgba[yi, xi, 3] = absent_rgb, 1.0
                absent_mask[yi, xi] = True
            continue
        n_up, n_down = st.get("significant_up", 0), st.get("significant_down", 0)
        n_sig = n_up + n_down
        if n_sig == 0:
            rgba[yi, xi, :3], rgba[yi, xi, 3] = notsig_rgb, 1.0
            continue
        frac = n_sig / n_tested
        # k / n / m  = significant / measured in this experiment / copies in
        # the strain genome. Always all three.
        if n_up > 0 and n_down > 0:
            # genes disagree -- flat orange, k split as up/down arrows
            rgba[yi, xi, :3], rgba[yi, xi, 3] = both_rgb, 1.0
            frac_txt[(yi, xi)] = (f"↑{n_up}↓{n_down}/{n_tested}/{n_genome}", 10.0)
        else:
            base = up_rgb if n_up >= n_down else down_rgb
            col = notsig_rgb + (base - notsig_rgb) * min(1.0, 0.25 + 0.75 * frac)
            rgba[yi, xi, :3], rgba[yi, xi, 3] = col, 1.0
            frac_txt[(yi, xi)] = (f"{n_sig}/{n_tested}/{n_genome}", 10.0)

    cell = 0.62
    fw = cell * n_cog + 8.5
    fh = cell * n_exp + 4.8
    fig, ax = plt.subplots(figsize=(fw, fh), dpi=300)
    fig.patch.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)
    ax.imshow(rgba, aspect="auto")

    # absent-from-genome cells: three diagonal strokes (a sparse hatch that
    # reads clearly against the near-white fill and the pale "no data" cell)
    for (yi, xi), on in np.ndenumerate(absent_mask):
        if on:
            x0, x1 = xi - 0.5, xi + 0.5
            yb, yt = yi + 0.5, yi - 0.5
            for (px0, py0), (px1, py1) in (((x0, yb), (x1, yt)),
                                           ((x0, yi), (xi, yt)),
                                           ((xi, yb), (x1, yi))):
                ax.plot([px0, px1], [py0, py1], color=HATCH_INK, linewidth=0.8, clip_on=True)
    for (yi, xi), (t, fs) in frac_txt.items():
        ax.text(xi, yi, t, ha="center", va="center", fontsize=fs, color="#1c1c1c",
                linespacing=0.9)

    # per-COG labels at the BOTTOM (protein + COG), all black
    ax.set_xticks(range(n_cog))
    ax.set_xticklabels(labels, rotation=50, ha="right", rotation_mode="anchor",
                       fontsize=13, color=INK)
    ax.xaxis.tick_bottom()
    ax.xaxis.set_label_position("bottom")
    ax.set_xlim(-0.5, n_cog - 0.5)
    ax.set_yticks(range(n_exp))
    ax.set_yticklabels([EXPERIMENT_LABELS[e] for e in ROW_ORDER], fontsize=15, color=INK)

    n_count = ROW_NUTRIENT.index("P")
    ax.axhline(n_count - 0.5, color=INK, linewidth=1.6)
    ax.set_xticks(np.arange(-0.5, n_cog, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, n_exp, 1), minor=True)
    ax.grid(which="minor", color=GRIDLINE, linewidth=1.0)
    ax.tick_params(which="both", length=0)
    for s in ax.spines.values():
        s.set_visible(False)

    # left-side bracket: nitrogen block vs phosphorus block
    ytrans = ax.get_yaxis_transform()  # x in axes fraction, y in data units
    x_bar, x_txt = -0.305, -0.325
    blocks = [("Nitrogen starvation", 0, n_count - 1),
              ("Phosphorus starvation", n_count, n_exp - 1)]
    for name, a, b in blocks:
        ax.plot([x_bar, x_bar], [a - 0.4, b + 0.4], color=INK, linewidth=1.4,
                transform=ytrans, clip_on=False)
        for y_end in (a - 0.4, b + 0.4):
            ax.plot([x_bar, x_bar + 0.010], [y_end, y_end], color=INK, linewidth=1.4,
                    transform=ytrans, clip_on=False)
        ax.text(x_txt, (a + b) / 2, name, rotation=90, ha="center", va="center",
                fontsize=15, color=INK, transform=ytrans, clip_on=False)

    # "General annotation" group brackets + category names, ABOVE the heatmap,
    # with a per-COG Direction strip (N / mixed, black) just under the brackets
    ax.set_ylim(n_exp - 0.5, -0.5)  # keep data orientation
    y_strip = -0.72
    for xi, d in enumerate(order["Direction"]):
        ax.text(xi, y_strip, "N" if d == "N" else "mixed", ha="center", va="center",
                fontsize=9.5, color=INK, clip_on=False)
    y_bracket = -1.30
    small_run = 0
    for a, b, cat in groups:
        ax.plot([a - 0.42, b + 0.42], [y_bracket, y_bracket], color=INK,
                linewidth=1.3, clip_on=False, solid_capstyle="butt")
        for x_end in (a - 0.42, b + 0.42):
            ax.plot([x_end, x_end], [y_bracket, y_bracket + 0.16], color=INK,
                    linewidth=1.3, clip_on=False)  # ticks point down toward the strip
        if a > 0:  # vertical divider between groups, through the heatmap
            ax.axvline(a - 0.5, color=INK, linewidth=1.0, alpha=0.5)
        # stagger only where consecutive narrow groups would otherwise collide
        is_small = (b - a) < 3
        small_run = small_run + 1 if is_small else 0
        raised = is_small and small_run % 2 == 0
        y_text = y_bracket - 0.36 - (0.60 if raised else 0.0)
        if raised:  # connector so the lifted label still reads as this bracket's
            ax.plot([(a + b) / 2, (a + b) / 2], [y_bracket, y_text - 0.05],
                    color=MUTED, linewidth=0.7, clip_on=False)
        ax.text((a + b) / 2, y_text, CATEGORY_DISPLAY.get(cat, cat),
                ha="center", va="bottom", fontsize=10, color=INK,
                linespacing=0.95, clip_on=False)

    legend = [
        Patch(facecolor=COLOR_UP, label="Upregulated (most of COG's genes)"),
        Patch(facecolor=COLOR_DOWN, label="Downregulated"),
        Patch(facecolor=COLOR_BOTH, label="Both up- and downregulated genes"),
        Patch(facecolor=COLOR_NOTSIG, label="Tested, no gene significant"),
        Patch(facecolor=COLOR_NODATA, edgecolor=MUTED, label="No data (COG present, not in table)"),
        Patch(facecolor=COLOR_ABSENT, edgecolor=HATCH_INK, hatch="///", label="COG absent from strain genome"),
    ]
    fig.legend(handles=legend, loc="lower center", ncol=6, frameon=False, fontsize=13,
               handlelength=1.4, handleheight=1.4, bbox_to_anchor=(0.5, 0.012))
    fig.suptitle("41 candidate starvation-sensitivity COGs vs. nutrient-starvation experiments",
                 x=0.5, y=0.985, ha="center", va="top", fontsize=19, color=INK)
    fig.text(0.5, 0.945,
             "Cell text is  k / n / m :  significant genes  /  genes measured in this experiment  /  copies of the COG in that strain's genome.  "
             "Colour = dominant direction, intensity scales with k/n;\n"
             "orange = the COG's genes split direction (k shown as up/down arrows).  Columns bracketed by functional category; "
             "N / mixed strip = the COG's Direction tag.  * Martiny MIT9313: 24h used (no 48h row).",
             ha="center", va="top", fontsize=11, color=MUTED, linespacing=1.5)
    fig.subplots_adjust(left=0.30, right=0.985, top=0.80, bottom=0.20)
    fig.savefig(FIG_DIR / "01_gene_experiment_heatmap.png", dpi=300, facecolor=SURFACE)
    print(f"Wrote figures/01_gene_experiment_heatmap.png  ({n_exp} x {n_cog})")


def make_bar(hit_rates, bg_pct):
    rows = [
        ("Background\n(random genes)", bg_pct, MUTED),
        ("N pooled", None, COLOR_UP), ("N: Dir N", None, COLOR_UP), ("N: Dir mixed", None, COLOR_UP),
        ("P pooled", None, COLOR_DOWN), ("P: Dir N", None, COLOR_DOWN), ("P: Dir mixed", None, COLOR_DOWN),
    ]
    lut = hit_rates.set_index("group")["pct_significant"].to_dict()
    keymap = {
        "N pooled": "nitrogen - all 41 pooled", "N: Dir N": "nitrogen - Direction N (13)",
        "N: Dir mixed": "nitrogen - Direction mixed (28)",
        "P pooled": "phosphorus - all 41 pooled", "P: Dir N": "phosphorus - Direction N (13)",
        "P: Dir mixed": "phosphorus - Direction mixed (28)",
    }
    labels, values, colors = [], [], []
    for lab, v, c in rows:
        labels.append(lab)
        values.append(bg_pct if v is not None else lut[keymap[lab]])
        colors.append(c)

    fig, ax = plt.subplots(figsize=(8.5, 4.5), dpi=300)
    fig.patch.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)
    bars = ax.bar(labels, values, color=colors)
    ax.axhline(bg_pct, color=MUTED, linewidth=1)
    for b, v in zip(bars, values):
        ax.text(b.get_x() + b.get_width() / 2, v + 0.4, f"{v:.1f}%", ha="center", fontsize=10)
    ax.set_ylabel("% of tests significant", fontsize=12)
    ax.tick_params(labelsize=11)
    ax.set_title("41-COG candidate list: starvation-response hit rate vs. background", fontsize=13)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "02_pct_significant_vs_background.png", dpi=300, facecolor=SURFACE)
    print("Wrote figures/02_pct_significant_vs_background.png")


def main():
    matrix = pd.read_csv(DATA_DIR / "01_target_gene_experiment_matrix.csv")
    hit_rates = pd.read_csv(DATA_DIR / "02_hit_rate_results.csv")
    cogs = pd.read_csv(STEP2_DATA / "00_source_normixed_cogs.csv")
    bg = pd.read_csv(FRAMING_DATA / "02_negative_background_pool.csv")
    bg_pct = 100 * bg["expression_status"].isin(["significant_up", "significant_down"]).mean()
    make_heatmap(matrix, cogs)
    make_bar(hit_rates, bg_pct)


if __name__ == "__main__":
    main()
