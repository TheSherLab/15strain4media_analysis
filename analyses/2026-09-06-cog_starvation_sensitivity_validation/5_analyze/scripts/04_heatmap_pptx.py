"""
Step 5 -- Figure 1 as an editable PowerPoint.

Same data and layout as 03_figures.py's heatmap (experiments x 41 COGs,
grouped by the researcher's "General annotation" category), but every cell
is its own rectangle + text box so the researcher can recolour, move, or
relabel anything in PowerPoint.

Cell fill: dominant direction (up red / down blue) blended toward grey by
the significant fraction; orange = the COG's genes split direction in that
cell; grey = tested but no gene significant; pale = no data (COG present,
not in this experiment's table); single diagonal line = COG absent from
the strain genome.
Cell text: k / n / m  (significant / measured here / genome copies). In a
split-direction (orange) cell, k is written "Nup Mdn" on its own line.

Row labels are the Oxford-style citation "First author et al., YEAR -
analysis type - strain"; a left-side bracket marks the nitrogen vs
phosphorus block.

Style (2026-09-09, researcher-requested): Arial throughout, no bold,
larger fonts, darker "tested-not-significant" grey.

Inputs: data/01_target_gene_experiment_matrix.csv,
        ../../2_kg_selection/data/00_source_normixed_cogs.csv
Outputs: figures/01_gene_experiment_heatmap.pptx

Usage (python-pptx pulled just for this, not a project dep):
  uv run --with python-pptx analyses/2026-09-06-cog_starvation_sensitivity_validation/5_analyze/scripts/04_heatmap_pptx.py
"""

from pathlib import Path

import pandas as pd
from pptx import Presentation
from pptx.dml.color import RGBColor
from pptx.enum.shapes import MSO_CONNECTOR, MSO_SHAPE
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR
from pptx.util import Inches, Pt

BASE = Path(__file__).resolve().parent.parent
DATA_DIR = BASE / "data"
FIG_DIR = BASE / "figures"
STEP2_DATA = BASE.parent / "2_kg_selection" / "data"

UP = (0xE3, 0x49, 0x48)
DOWN = (0x2A, 0x78, 0xD6)
BOTH = (0xE8, 0x91, 0x3A)          # orange -- genes split up/down in one cell
NS = (0xA8, 0xA6, 0x9B)            # darker than the walkthrough (#c3c2b7)
NODATA = (0xE7, 0xE4, 0xD8)        # distinctly darker beige than the absent cell
ABSENT_BG = RGBColor(0xFD, 0xFC, 0xF9)
ABSENT_LINE = RGBColor(0x9A, 0x98, 0x8E)
GRIDLINE = RGBColor(0xD7, 0xD5, 0xC9)   # visible cell borders
INK = RGBColor(0x0B, 0x0B, 0x0B)
WHITE = RGBColor(0xFF, 0xFF, 0xFF)
MUTED = RGBColor(0x89, 0x87, 0x81)
SURFACE = RGBColor(0xFC, 0xFC, 0xFB)
TESTED = {"significant_up", "significant_down", "not_significant"}

CATEGORY_DISPLAY = {
    "Quality contol- DNA level": "Quality control - DNA level",
    "Quality contol- Protein level": "Quality control - Protein level",
    "Amino acid and mixotrophy related biosynthesis": "Amino acid & mixotrophy biosynthesis",
    "LPS": "LPS",
    "Membrane linker": "Membrane linker",
    "Exopolysaccharide": "Exopolysaccharide",
    "Biofilm/Attachment": "Biofilm / Attachment",
    "Energy production": "Energy production",
    "Translation, transcription and signal transduction": "Translation, transcription & signalling",
    "RTX toxins": "RTX toxins",
    "mixeds": "mixed",
}

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

SLIDE_W, SLIDE_H = 24.0, 9.4
BRACKET_X = 0.65
LEFT_LABELS = 5.2
GRID_LEFT = LEFT_LABELS + 0.1
GRID_TOP = 3.05
CELL_H = 0.46


def _blend(base, frac):
    f = min(1.0, 0.25 + 0.75 * frac)
    return RGBColor(*(round(NS[i] + (base[i] - NS[i]) * f) for i in range(3)))


def _line(slide, x1, y1, x2, y2, color, width_pt):
    ln = slide.shapes.add_connector(MSO_CONNECTOR.STRAIGHT, Inches(x1), Inches(y1),
                                    Inches(x2), Inches(y2))
    ln.line.color.rgb = color
    ln.line.width = Pt(width_pt)
    ln.shadow.inherit = False
    return ln


def _text(slide, left, top, w, h, text, size, *, color=INK,
          align=PP_ALIGN.LEFT, anchor=MSO_ANCHOR.MIDDLE, rot=0, font="Arial"):
    tb = slide.shapes.add_textbox(Inches(left), Inches(top), Inches(w), Inches(h))
    tf = tb.text_frame
    tf.word_wrap = False
    tf.margin_left = tf.margin_right = tf.margin_top = tf.margin_bottom = 0
    tf.vertical_anchor = anchor
    p = tf.paragraphs[0]
    p.alignment = align
    r = p.add_run()
    r.text = text
    r.font.size = Pt(size)
    r.font.bold = False
    r.font.name = font
    r.font.color.rgb = color
    if rot:
        tb.rotation = rot
    return tb


def cell_spec(g):
    """(category, fill_or_None, text) for one (COG, experiment) group."""
    st = g["status"].value_counts().to_dict()
    n_tested = sum(st.get(s, 0) for s in TESTED)
    n_genome = int(g["locus_tag"].notna().sum())
    if n_tested == 0:
        if st.get("no_data_at_timepoint", 0) > 0:
            return "nodata", RGBColor(*NODATA), ""
        return "absent", None, ""
    n_up, n_down = st.get("significant_up", 0), st.get("significant_down", 0)
    n_sig = n_up + n_down
    if n_sig == 0:
        return "ns", RGBColor(*NS), f"{n_sig}/{n_tested}/{n_genome}"
    if n_up > 0 and n_down > 0:
        return "both", RGBColor(*BOTH), f"↑{n_up}↓{n_down}/{n_tested}/{n_genome}"
    base = UP if n_up >= n_down else DOWN
    return ("up" if base is UP else "down"), _blend(base, n_sig / n_tested), f"{n_sig}/{n_tested}/{n_genome}"


def main() -> None:
    matrix = pd.read_csv(DATA_DIR / "01_target_gene_experiment_matrix.csv")
    cogs = pd.read_csv(STEP2_DATA / "00_source_normixed_cogs.csv").sort_values("order").reset_index(drop=True)
    cogs["General annotation"] = cogs["General annotation"].ffill()
    cog_order = cogs["COGs"].tolist()

    def label(row):
        name = "" if pd.isna(row["Protien"]) or not str(row["Protien"]).strip() else str(row["Protien"]).strip()
        name = name.split(" (")[0].split(",")[0].strip()
        return f"{name or row['COGs']}  {row['COGs']}"
    col_labels = cogs.apply(label, axis=1).tolist()
    col_dir = cogs["Direction"].tolist()

    groups, start = [], 0
    cats = cogs["General annotation"].tolist()
    for i in range(1, len(cats) + 1):
        if i == len(cats) or cats[i] != cats[start]:
            groups.append((start, i - 1, cats[start]))
            start = i

    specs = {}
    for (cog, exp), g in matrix.groupby(["cog_number", "experiment_id"]):
        specs[(exp, cog)] = cell_spec(g)

    cw = (SLIDE_W - GRID_LEFT - 0.3) / len(cog_order)

    prs = Presentation()
    prs.slide_width, prs.slide_height = Inches(SLIDE_W), Inches(SLIDE_H)
    slide = prs.slides.add_slide(prs.slide_layouts[6])
    bg = slide.shapes.add_shape(MSO_SHAPE.RECTANGLE, 0, 0, prs.slide_width, prs.slide_height)
    bg.fill.solid(); bg.fill.fore_color.rgb = SURFACE
    bg.line.fill.background(); bg.shadow.inherit = False

    _text(slide, 0.4, 0.10, SLIDE_W - 0.8, 0.44,
          "41 candidate starvation-sensitivity COGs vs. nutrient-starvation experiments",
          22, align=PP_ALIGN.CENTER)
    _text(slide, 0.4, 0.56, SLIDE_W - 0.8, 0.50,
          "cell text k / n / m = significant / measured in this experiment / copies in that strain's genome  ·  "
          "fill = dominant direction, blended toward grey by k/n  ·  orange = the COG's genes split direction (k shown as up/down arrows)  ·  "
          "columns grouped by functional category, N / mixed strip = that COG's Direction tag  ·  "
          "* Martiny MIT9313 at 24h (no 48h row); table-absent COGs = not significant",
          12, color=MUTED, align=PP_ALIGN.CENTER)

    # group brackets + category labels, with a per-COG Direction strip beneath
    for a, b, cat in groups:
        x0 = GRID_LEFT + a * cw
        x1 = GRID_LEFT + (b + 1) * cw
        bar = slide.shapes.add_shape(MSO_SHAPE.RECTANGLE, Inches(x0 + 0.03), Inches(GRID_TOP - 0.55),
                                     Inches(x1 - x0 - 0.06), Inches(0.022))
        bar.fill.solid(); bar.fill.fore_color.rgb = INK
        bar.line.fill.background(); bar.shadow.inherit = False
        _text(slide, x0 - 0.6, GRID_TOP - 1.25, (x1 - x0) + 1.2, 0.55,
              CATEGORY_DISPLAY.get(cat, cat), 12, color=INK, align=PP_ALIGN.CENTER,
              anchor=MSO_ANCHOR.BOTTOM)

    # Direction strip: N / mixed per COG, black, just under the brackets
    for j, d in enumerate(col_dir):
        _text(slide, GRID_LEFT + j * cw, GRID_TOP - 0.36, cw, 0.28,
              "N" if d == "N" else "mixed", 10, color=INK, align=PP_ALIGN.CENTER)

    # column labels (rotated), all black
    for j, lab in enumerate(col_labels):
        cx = GRID_LEFT + j * cw + cw / 2
        _text(slide, cx - 1.0, GRID_TOP + len(ROW_ORDER) * CELL_H + 0.15, 2.0, 0.24, lab, 12,
              color=INK, align=PP_ALIGN.LEFT, rot=300)

    # rows: label + 41 cells
    for i, exp in enumerate(ROW_ORDER):
        top = GRID_TOP + i * CELL_H
        _text(slide, BRACKET_X + 0.35, top, LEFT_LABELS - BRACKET_X - 0.5, CELL_H,
              EXPERIMENT_LABELS[exp], 14, color=INK, align=PP_ALIGN.RIGHT)
        for j, cog in enumerate(cog_order):
            cat, fill, txt = specs.get((exp, cog), ("absent", None, ""))
            left = GRID_LEFT + j * cw
            shp = slide.shapes.add_shape(MSO_SHAPE.RECTANGLE, Inches(left), Inches(top), Inches(cw), Inches(CELL_H))
            shp.name = f"cell|{EXPERIMENT_LABELS[exp]}|{cog}|{cat}"
            shp.shadow.inherit = False
            shp.line.color.rgb = GRIDLINE
            shp.line.width = Pt(1.0)
            if fill is None:
                shp.fill.solid(); shp.fill.fore_color.rgb = ABSENT_BG
                _line(slide, left, top + CELL_H, left + cw, top, ABSENT_LINE, 0.8)
                _line(slide, left, top + CELL_H / 2, left + cw / 2, top, ABSENT_LINE, 0.8)
                _line(slide, left + cw / 2, top + CELL_H, left + cw, top + CELL_H / 2, ABSENT_LINE, 0.8)
            else:
                shp.fill.solid(); shp.fill.fore_color.rgb = fill
            if txt:
                _text(slide, left, top, cw, CELL_H, txt, 10,
                      color=INK, align=PP_ALIGN.CENTER, anchor=MSO_ANCHOR.MIDDLE)

    grid_bottom = GRID_TOP + len(ROW_ORDER) * CELL_H
    grid_right = GRID_LEFT + len(cog_order) * cw

    # left bracket: nitrogen block vs phosphorus block
    n_rows = ROW_NUTRIENT.index("P")
    for name, a, b in [("Nitrogen starvation", 0, n_rows - 1),
                       ("Phosphorus starvation", n_rows, len(ROW_ORDER) - 1)]:
        y0 = GRID_TOP + a * CELL_H + 0.03
        y1 = GRID_TOP + (b + 1) * CELL_H - 0.03
        _line(slide, BRACKET_X, y0, BRACKET_X, y1, INK, 1.6)
        _line(slide, BRACKET_X, y0, BRACKET_X + 0.10, y0, INK, 1.6)
        _line(slide, BRACKET_X, y1, BRACKET_X + 0.10, y1, INK, 1.6)
        tb = _text(slide, BRACKET_X - 1.05, (y0 + y1) / 2 - 0.15, 2.1, 0.3, name, 14,
                   color=INK, align=PP_ALIGN.CENTER)
        tb.rotation = 270

    n_rows_count = n_rows
    hline = slide.shapes.add_shape(MSO_SHAPE.RECTANGLE, Inches(GRID_LEFT),
                                   Inches(GRID_TOP + n_rows_count * CELL_H - 0.012),
                                   Inches(grid_right - GRID_LEFT), Inches(0.028))
    hline.fill.solid(); hline.fill.fore_color.rgb = INK
    hline.line.fill.background(); hline.shadow.inherit = False
    for a, b, _cat in groups[1:]:
        vx = GRID_LEFT + a * cw
        vln = slide.shapes.add_shape(MSO_SHAPE.RECTANGLE, Inches(vx - 0.012), Inches(GRID_TOP),
                                     Inches(0.02), Inches(grid_bottom - GRID_TOP))
        vln.fill.solid(); vln.fill.fore_color.rgb = RGBColor(0x9A, 0x99, 0x94)
        vln.line.fill.background(); vln.shadow.inherit = False

    # legend
    ly = grid_bottom + 1.7
    legend = [(RGBColor(*UP), "Upregulated (most of COG's genes)"),
              (RGBColor(*DOWN), "Downregulated"),
              (RGBColor(*BOTH), "Both up- and downregulated genes"),
              (RGBColor(*NS), "Tested, no gene significant"),
              (RGBColor(*NODATA), "No data (COG present, not in table)"),
              (None, "COG absent from strain genome")]
    lx = GRID_LEFT
    for clr, lab in legend:
        sw = slide.shapes.add_shape(MSO_SHAPE.RECTANGLE, Inches(lx), Inches(ly), Inches(0.30), Inches(0.30))
        sw.shadow.inherit = False; sw.line.color.rgb = MUTED; sw.line.width = Pt(0.75)
        if clr is None:
            sw.fill.solid(); sw.fill.fore_color.rgb = ABSENT_BG
            _line(slide, lx, ly + 0.30, lx + 0.30, ly, ABSENT_LINE, 0.8)
            _line(slide, lx, ly + 0.15, lx + 0.15, ly, ABSENT_LINE, 0.8)
            _line(slide, lx + 0.15, ly + 0.30, lx + 0.30, ly + 0.15, ABSENT_LINE, 0.8)
        else:
            sw.fill.solid(); sw.fill.fore_color.rgb = clr
        _text(slide, lx + 0.42, ly - 0.02, 4.0, 0.34, lab, 13)
        lx += 0.42 + 0.072 * len(lab) + 0.6

    out = FIG_DIR / "01_gene_experiment_heatmap.pptx"
    prs.save(out)
    print(f"Wrote {out}  ({len(ROW_ORDER)} rows x {len(cog_order)} COGs = "
          f"{len(ROW_ORDER) * len(cog_order)} editable cells)")


if __name__ == "__main__":
    main()
