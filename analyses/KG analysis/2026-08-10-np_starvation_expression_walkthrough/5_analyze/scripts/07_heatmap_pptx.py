"""
Step 5 -- Figure 1 as an editable PowerPoint.

Same data and layout as 03_heatmap_figure.py, but every cell is its own
rectangle shape and every label its own text box, so the researcher can
recolour individual cells, move labels, or relabel in PowerPoint. Reads
the same frozen matrix -- rerun after 01_extract_target_de.py changes.

Cell fill:
  up   = red     down = blue    tested-not-significant = grey
  no data (gene present, not in table) = pale
  gene absent from strain genome       = single diagonal line

Row labels are the Oxford-style citation "First author et al., YEAR -
analysis type - strain"; a left-side bracket marks the nitrogen vs
phosphorus block.

Style (2026-09-09, researcher-requested): Arial throughout, no bold,
larger fonts, darker "tested-not-significant" grey.

Inputs: data/01_target_gene_experiment_matrix.csv
Outputs: figures/01_gene_experiment_heatmap.pptx

Usage (python-pptx is not a project dependency -- pulled just for this):
  uv run --with python-pptx analyses/2026-08-10-np_starvation_expression_walkthrough/5_analyze/scripts/07_heatmap_pptx.py
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

COLOR_UP = RGBColor(0xE3, 0x49, 0x48)
COLOR_DOWN = RGBColor(0x2A, 0x78, 0xD6)
COLOR_NS = RGBColor(0xA8, 0xA6, 0x9B)       # darker than the old #c3c2b7
COLOR_NODATA = RGBColor(0xE7, 0xE4, 0xD8)   # distinctly darker beige than the absent cell
COLOR_ABSENT_BG = RGBColor(0xFD, 0xFC, 0xF9)
COLOR_ABSENT_LINE = RGBColor(0x9A, 0x98, 0x8E)
GRIDLINE = RGBColor(0xD7, 0xD5, 0xC9)       # visible cell borders
INK = RGBColor(0x0B, 0x0B, 0x0B)
N_BLUE = RGBColor(0x2A, 0x78, 0xD6)
P_RED = RGBColor(0xE3, 0x49, 0x48)
WHITE = RGBColor(0xFF, 0xFF, 0xFF)
MUTED = RGBColor(0x89, 0x87, 0x81)
SURFACE = RGBColor(0xFC, 0xFC, 0xFB)

STATUS_TO_CAT = {
    "significant_up": "up",
    "significant_down": "down",
    "not_significant": "ns",
    "no_data_at_timepoint": "nodata",
    "no_locus_in_strain": "absent",
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

# geometry (inches)
SLIDE_W, SLIDE_H = 22.5, 8.6
BRACKET_X = 0.65
LEFT_LABELS = 5.4
TOP_TITLE = 0.66
TOP_COLLABELS = 1.85
GRID_LEFT = LEFT_LABELS + 0.1
GRID_TOP = TOP_TITLE + TOP_COLLABELS
CELL_H = 0.46


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


def _cell(slide, left, top, w, h, cat, name):
    shp = slide.shapes.add_shape(MSO_SHAPE.RECTANGLE, Inches(left), Inches(top), Inches(w), Inches(h))
    shp.name = name
    shp.shadow.inherit = False
    shp.line.color.rgb = GRIDLINE
    shp.line.width = Pt(1.0)
    if cat == "absent":
        shp.fill.solid()
        shp.fill.fore_color.rgb = COLOR_ABSENT_BG
        _line(slide, left, top + h, left + w, top, COLOR_ABSENT_LINE, 0.8)
        _line(slide, left, top + h / 2, left + w / 2, top, COLOR_ABSENT_LINE, 0.8)
        _line(slide, left + w / 2, top + h, left + w, top + h / 2, COLOR_ABSENT_LINE, 0.8)
    else:
        shp.fill.solid()
        shp.fill.fore_color.rgb = {
            "up": COLOR_UP, "down": COLOR_DOWN, "ns": COLOR_NS, "nodata": COLOR_NODATA,
        }[cat]
    return shp


def main() -> None:
    matrix = pd.read_csv(DATA_DIR / "01_target_gene_experiment_matrix.csv")
    matrix["cat"] = matrix["status"].map(STATUS_TO_CAT)

    has_signal = matrix.groupby("gene_name")["cat"].apply(lambda s: s.isin(["up", "down"]).any())
    kept = set(has_signal[has_signal].index)
    n_or_p = matrix.drop_duplicates("gene_name").set_index("gene_name")["n_or_p"]
    genes = sorted(kept, key=lambda g: (n_or_p[g] != "N", g))
    n_gene_count = sum(1 for g in genes if n_or_p[g] == "N")

    cell_map = {
        (r["experiment_id"], r["gene_name"]): r["cat"]
        for _, r in matrix.iterrows()
        if r["gene_name"] in kept
    }

    cw = (SLIDE_W - GRID_LEFT - 0.3) / len(genes)

    prs = Presentation()
    prs.slide_width, prs.slide_height = Inches(SLIDE_W), Inches(SLIDE_H)
    slide = prs.slides.add_slide(prs.slide_layouts[6])  # blank
    bg = slide.shapes.add_shape(MSO_SHAPE.RECTANGLE, 0, 0, prs.slide_width, prs.slide_height)
    bg.fill.solid()
    bg.fill.fore_color.rgb = SURFACE
    bg.line.fill.background()
    bg.shadow.inherit = False

    _text(slide, 0.4, 0.08, SLIDE_W - 0.8, 0.40,
          "Nitrogen- and phosphorus-acquisition gene response to nutrient starvation",
          22, align=PP_ALIGN.CENTER)
    _text(slide, 0.4, 0.48, SLIDE_W - 0.8, 0.26,
          "genes with a response in >=1 experiment; each experiment at its single chosen starvation timepoint",
          12, color=MUTED, align=PP_ALIGN.CENTER)

    # group headers
    _text(slide, GRID_LEFT, GRID_TOP - 1.5, cw * n_gene_count, 0.32,
          "N acquisition genes", 15, color=N_BLUE, align=PP_ALIGN.CENTER)
    _text(slide, GRID_LEFT + cw * n_gene_count, GRID_TOP - 1.5, cw * (len(genes) - n_gene_count), 0.32,
          "P acquisition genes", 15, color=P_RED, align=PP_ALIGN.CENTER)

    # column labels (rotated)
    for j, g in enumerate(genes):
        cx = GRID_LEFT + j * cw + cw / 2
        _text(slide, cx - 0.75, GRID_TOP - 1.18, 1.5, 0.24, g, 12,
              color=INK, align=PP_ALIGN.LEFT, rot=300)

    # row labels + cells
    for i, exp in enumerate(ROW_ORDER):
        top = GRID_TOP + i * CELL_H
        _text(slide, BRACKET_X + 0.35, top, LEFT_LABELS - BRACKET_X - 0.5, CELL_H,
              EXPERIMENT_LABELS[exp], 14, color=INK, align=PP_ALIGN.RIGHT)
        for j, g in enumerate(genes):
            cat = cell_map.get((exp, g), "nodata")
            _cell(slide, GRID_LEFT + j * cw, top, cw, CELL_H, cat,
                  f"cell|{EXPERIMENT_LABELS[exp]}|{g}|{cat}")

    grid_bottom = GRID_TOP + len(ROW_ORDER) * CELL_H
    grid_right = GRID_LEFT + len(genes) * cw

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

    # N/P dividers
    hline = slide.shapes.add_shape(MSO_SHAPE.RECTANGLE, Inches(GRID_LEFT),
                                   Inches(GRID_TOP + n_rows * CELL_H - 0.01),
                                   Inches(grid_right - GRID_LEFT), Inches(0.025))
    hline.fill.solid(); hline.fill.fore_color.rgb = INK; hline.line.fill.background()
    hline.shadow.inherit = False
    vline = slide.shapes.add_shape(MSO_SHAPE.RECTANGLE,
                                   Inches(GRID_LEFT + n_gene_count * cw - 0.01), Inches(GRID_TOP),
                                   Inches(0.025), Inches(grid_bottom - GRID_TOP))
    vline.fill.solid(); vline.fill.fore_color.rgb = INK; vline.line.fill.background()
    vline.shadow.inherit = False

    # legend
    ly = grid_bottom + 0.30
    legend = [("up", "Upregulated"), ("down", "Downregulated"), ("ns", "Tested, not significant"),
              ("nodata", "No data (gene present, not in table)"), ("absent", "Gene absent from strain genome")]
    lx = GRID_LEFT
    for cat, lab in legend:
        _cell(slide, lx, ly, 0.32, 0.32, cat, f"legend|{cat}")
        _text(slide, lx + 0.44, ly - 0.02, 4.0, 0.36, lab, 13)
        lx += 0.6 + 0.075 * len(lab) + 0.55

    _text(slide, GRID_LEFT, ly + 0.6, SLIDE_W - GRID_LEFT - 0.3, 0.24,
          "* Martiny MIT9313: significance-criterion / timepoint mismatch (24h used, not 48h) - see 5_analyze/notebook.md. "
          "Martiny 'tested, not significant' cells include genes absent from that publication's filtered table.",
          10.5, color=MUTED)

    out = FIG_DIR / "01_gene_experiment_heatmap.pptx"
    prs.save(out)
    print(f"Wrote {out}  ({len(ROW_ORDER)} rows x {len(genes)} genes = {len(ROW_ORDER) * len(genes)} editable cells)")


if __name__ == "__main__":
    main()
