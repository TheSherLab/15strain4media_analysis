"""
Step 5 -- Figure 1 as an editable PowerPoint.

Same data and layout as 03_heatmap_figure.py, but every cell is its own
rectangle shape and every label its own text box, so the researcher can
recolour individual cells, move labels, or relabel in PowerPoint. Reads
the same frozen matrix -- rerun after 01_extract_target_de.py changes.

Cell fill:
  up   = red     down = blue    tested-not-significant = grey
  no data (gene present, not in table) = pale
  gene absent from strain genome       = diagonal pattern fill

Inputs: data/01_target_gene_experiment_matrix.csv
Outputs: figures/01_gene_experiment_heatmap.pptx

Usage (python-pptx is not a project dependency -- pulled just for this):
  uv run --with python-pptx analyses/2026-08-10-np_starvation_expression_walkthrough/5_analyze/scripts/07_heatmap_pptx.py
"""

from pathlib import Path

import pandas as pd
from pptx import Presentation
from pptx.dml.color import RGBColor
from pptx.enum.shapes import MSO_SHAPE
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR
from pptx.oxml.ns import qn
from pptx.util import Inches, Pt

BASE = Path(__file__).resolve().parent.parent
DATA_DIR = BASE / "data"
FIG_DIR = BASE / "figures"

COLOR_UP = RGBColor(0xE3, 0x49, 0x48)
COLOR_DOWN = RGBColor(0x2A, 0x78, 0xD6)
COLOR_NS = RGBColor(0xC3, 0xC2, 0xB7)
COLOR_NODATA = RGBColor(0xF4, 0xF3, 0xEF)
COLOR_ABSENT_FG = RGBColor(0xB7, 0xB5, 0xAC)
COLOR_ABSENT_BG = RGBColor(0xFA, 0xF9, 0xF5)
INK = RGBColor(0x0B, 0x0B, 0x0B)
N_BLUE = RGBColor(0x2A, 0x78, 0xD6)
P_RED = RGBColor(0xE3, 0x49, 0x48)
WHITE = RGBColor(0xFF, 0xFF, 0xFF)
SURFACE = RGBColor(0xFC, 0xFC, 0xFB)

STATUS_TO_CAT = {
    "significant_up": "up",
    "significant_down": "down",
    "not_significant": "ns",
    "no_data_at_timepoint": "nodata",
    "no_locus_in_strain": "absent",
}

EXPERIMENT_LABELS = {
    "10.1101/2025.11.24.690089_growth_state_pro99lown_nutrient_starvation_med4_proteomics_axenic": "N: Weissberg — Proteomics (MED4)",
    "10.1101/2025.11.24.690089_growth_state_pro99lown_nutrient_starvation_med4_rnaseq_axenic": "N: Weissberg — RNA-seq (MED4)",
    "10.1038/ismej.2017.88_nitrogen_stress_ndepleted_pro99_medium_med4_rnaseq": "N: Read — RNA-seq (MED4)",
    "10.1038/msb4100087_nitrogen_nitrogen_deprivation_med4_med4_microarray": "N: Tolonen — Microarray (MED4)",
    "10.1038/msb4100087_nitrogen_nitrogen_deprivation_mit9313_mit9313_microarray": "N: Tolonen — Microarray (MIT9313)",
    "10.1111/1462-2920.13104_phosphorus_plimited_natl2a_rnaseq_uninfected": "P: Lin — RNA-seq (NATL2A), 59h",
    "10.1186/2046-9063-8-7_pi_limitation_mit9312_itraq": "P: Fuszard — Proteomics (MIT9312)",
    "10.1186/2046-9063-8-7_pi_limitation_natl2a_itraq": "P: Fuszard — Proteomics (NATL2A)",
    "10.1073/pnas.0601301103_phosphorus_phosphate_starvation_med4_microarray": "P: Martiny — Microarray (MED4)",
    "10.1073/pnas.0601301103_phosphorus_phosphate_starvation_mit9313_microarray": "P: Martiny — Microarray (MIT9313)*",
}
ROW_ORDER = list(EXPERIMENT_LABELS)

# geometry (inches)
SLIDE_W, SLIDE_H = 20.0, 8.0
LEFT_LABELS = 3.5
TOP_TITLE = 0.55
TOP_COLLABELS = 1.55
GRID_LEFT = LEFT_LABELS + 0.1
GRID_TOP = TOP_TITLE + TOP_COLLABELS
CELL_H = 0.40


def _set_pattern_fill(shape, pattern: str, fg: RGBColor, bg: RGBColor) -> None:
    """python-pptx has no high-level pattern-fill API; write the XML directly."""
    sp = shape.fill._xPr  # <a:spPr>
    for tag in ("a:noFill", "a:solidFill", "a:gradFill", "a:blipFill", "a:pattFill", "a:grpFill"):
        for el in sp.findall(qn(tag)):
            sp.remove(el)
    patt = sp.makeelement(qn("a:pattFill"), {"prst": pattern})
    fg_el = patt.makeelement(qn("a:fgClr"), {})
    fg_srgb = fg_el.makeelement(qn("a:srgbClr"), {"val": "%02X%02X%02X" % (fg[0], fg[1], fg[2])})
    fg_el.append(fg_srgb)
    bg_el = patt.makeelement(qn("a:bgClr"), {})
    bg_srgb = bg_el.makeelement(qn("a:srgbClr"), {"val": "%02X%02X%02X" % (bg[0], bg[1], bg[2])})
    bg_el.append(bg_srgb)
    patt.append(fg_el)
    patt.append(bg_el)
    ln = sp.find(qn("a:ln"))
    sp.insert(list(sp).index(ln) if ln is not None else len(sp), patt)


def _cell(slide, left, top, w, h, cat, name):
    shp = slide.shapes.add_shape(MSO_SHAPE.RECTANGLE, Inches(left), Inches(top), Inches(w), Inches(h))
    shp.name = name
    shp.shadow.inherit = False
    shp.line.color.rgb = WHITE
    shp.line.width = Pt(1.0)
    if cat == "absent":
        _set_pattern_fill(shp, "ltUpDiag", COLOR_ABSENT_FG, COLOR_ABSENT_BG)
    else:
        shp.fill.solid()
        shp.fill.fore_color.rgb = {
            "up": COLOR_UP, "down": COLOR_DOWN, "ns": COLOR_NS, "nodata": COLOR_NODATA,
        }[cat]
    return shp


def _text(slide, left, top, w, h, text, size, *, bold=False, color=INK,
          align=PP_ALIGN.LEFT, anchor=MSO_ANCHOR.MIDDLE, rot=0, font="Calibri"):
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
    r.font.bold = bold
    r.font.name = font
    r.font.color.rgb = color
    if rot:
        tb.rotation = rot
    return tb


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

    _text(slide, 0.4, 0.08, SLIDE_W - 0.8, 0.32,
          "Nitrogen- and phosphorus-acquisition gene response to nutrient starvation",
          16, bold=True, align=PP_ALIGN.CENTER)
    _text(slide, 0.4, 0.40, SLIDE_W - 0.8, 0.20,
          "genes with a response in ≥1 experiment; each experiment at its single chosen starvation timepoint",
          10, color=RGBColor(0x52, 0x51, 0x4E), align=PP_ALIGN.CENTER)

    # group headers
    _text(slide, GRID_LEFT, GRID_TOP - 1.5, cw * n_gene_count, 0.28,
          "N acquisition genes", 12, bold=True, color=N_BLUE, align=PP_ALIGN.CENTER)
    _text(slide, GRID_LEFT + cw * n_gene_count, GRID_TOP - 1.5, cw * (len(genes) - n_gene_count), 0.28,
          "P acquisition genes", 12, bold=True, color=P_RED, align=PP_ALIGN.CENTER)

    # column labels (rotated)
    for j, g in enumerate(genes):
        cx = GRID_LEFT + j * cw + cw / 2
        _text(slide, cx - 0.75, GRID_TOP - 1.18, 1.5, 0.24, g, 8,
              color=N_BLUE if n_or_p[g] == "N" else P_RED,
              align=PP_ALIGN.LEFT, rot=300, font="Consolas")

    # row labels + cells
    for i, exp in enumerate(ROW_ORDER):
        top = GRID_TOP + i * CELL_H
        _text(slide, 0.2, top, LEFT_LABELS - 0.3, CELL_H, EXPERIMENT_LABELS[exp], 10,
              color=RGBColor(0x52, 0x51, 0x4E), align=PP_ALIGN.RIGHT)
        for j, g in enumerate(genes):
            cat = cell_map.get((exp, g), "nodata")
            _cell(slide, GRID_LEFT + j * cw, top, cw, CELL_H, cat,
                  f"cell|{EXPERIMENT_LABELS[exp]}|{g}|{cat}")

    grid_bottom = GRID_TOP + len(ROW_ORDER) * CELL_H
    grid_right = GRID_LEFT + len(genes) * cw

    # N/P dividers
    n_rows = sum(1 for e in ROW_ORDER if EXPERIMENT_LABELS[e].startswith("N:"))
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
    ly = grid_bottom + 0.25
    legend = [("up", "Upregulated"), ("down", "Downregulated"), ("ns", "Tested, not significant"),
              ("nodata", "No data (gene present, not in table)"), ("absent", "Gene absent from strain genome")]
    lx = GRID_LEFT
    for cat, lab in legend:
        _cell(slide, lx, ly, 0.28, 0.28, cat, f"legend|{cat}")
        _text(slide, lx + 0.36, ly - 0.02, 3.1, 0.32, lab, 10)
        lx += 0.5 + 0.062 * len(lab) + 0.55

    _text(slide, GRID_LEFT, ly + 0.5, SLIDE_W - GRID_LEFT - 0.3, 0.24,
          "* Martiny MIT9313: significance-criterion / timepoint mismatch (24h used, not 48h) — see 5_analyze/notebook.md. "
          "Martiny 'tested, not significant' cells include genes absent from that publication's filtered table.",
          8.5, color=RGBColor(0x89, 0x87, 0x81))

    out = FIG_DIR / "01_gene_experiment_heatmap.pptx"
    prs.save(out)
    print(f"Wrote {out}  ({len(ROW_ORDER)} rows x {len(genes)} genes = {len(ROW_ORDER) * len(genes)} editable cells)")


if __name__ == "__main__":
    main()
