"""
Step 5 -- Figure 1 as an editable PowerPoint.

Same data and layout as 03_figures.py's heatmap (experiments x 41 COGs,
grouped by the researcher's "General annotation" category), but every cell
is its own rectangle + text box so the researcher can recolour, move, or
relabel anything in PowerPoint.

Cell fill: dominant direction (up red / down blue) blended toward grey by
the significant fraction; grey = tested but no gene significant; pale = no
data (COG present, not in this experiment's table); diagonal pattern = COG
absent from the strain genome.
Cell text: k / n / m  (significant / measured here / genome copies), always
all three.

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
from pptx.enum.shapes import MSO_SHAPE
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR
from pptx.oxml.ns import qn
from pptx.util import Inches, Pt

BASE = Path(__file__).resolve().parent.parent
DATA_DIR = BASE / "data"
FIG_DIR = BASE / "figures"
STEP2_DATA = BASE.parent / "2_kg_selection" / "data"

UP = (0xE3, 0x49, 0x48)
DOWN = (0x2A, 0x78, 0xD6)
NS = (0xC3, 0xC2, 0xB7)
NODATA = (0xF4, 0xF3, 0xEF)
ABSENT_FG = RGBColor(0xB7, 0xB5, 0xAC)
ABSENT_BG = RGBColor(0xFA, 0xF9, 0xF5)
INK = RGBColor(0x0B, 0x0B, 0x0B)
DIR_BLUE = RGBColor(0x2A, 0x78, 0xD6)
DIR_ORANGE = RGBColor(0xB5, 0x48, 0x0F)
WHITE = RGBColor(0xFF, 0xFF, 0xFF)
SECONDARY = RGBColor(0x52, 0x51, 0x4E)
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
    "10.1101/2025.11.24.690089_growth_state_pro99lown_nutrient_starvation_med4_proteomics_axenic": "N: Weissberg - Proteomics (MED4)",
    "10.1101/2025.11.24.690089_growth_state_pro99lown_nutrient_starvation_med4_rnaseq_axenic": "N: Weissberg - RNA-seq (MED4)",
    "10.1038/ismej.2017.88_nitrogen_stress_ndepleted_pro99_medium_med4_rnaseq": "N: Read - RNA-seq (MED4)",
    "10.1038/msb4100087_nitrogen_nitrogen_deprivation_med4_med4_microarray": "N: Tolonen - Microarray (MED4)",
    "10.1038/msb4100087_nitrogen_nitrogen_deprivation_mit9313_mit9313_microarray": "N: Tolonen - Microarray (MIT9313)",
    "10.1111/1462-2920.13104_phosphorus_plimited_natl2a_rnaseq_uninfected": "P: Lin - RNA-seq (NATL2A), 59h",
    "10.1186/2046-9063-8-7_pi_limitation_mit9312_itraq": "P: Fuszard - Proteomics (MIT9312)",
    "10.1186/2046-9063-8-7_pi_limitation_natl2a_itraq": "P: Fuszard - Proteomics (NATL2A)",
    "10.1073/pnas.0601301103_phosphorus_phosphate_starvation_med4_microarray": "P: Martiny - Microarray (MED4)",
    "10.1073/pnas.0601301103_phosphorus_phosphate_starvation_mit9313_microarray": "P: Martiny - Microarray (MIT9313)*",
}
ROW_ORDER = list(EXPERIMENT_LABELS)

SLIDE_W, SLIDE_H = 22.0, 9.2
LEFT_LABELS = 3.4
GRID_LEFT = LEFT_LABELS + 0.1
GRID_TOP = 2.85
CELL_H = 0.42


def _blend(base, frac):
    f = min(1.0, 0.25 + 0.75 * frac)
    return RGBColor(*(round(NS[i] + (base[i] - NS[i]) * f) for i in range(3)))


def _pattern_fill(shape, prst, fg, bg):
    sp = shape.fill._xPr
    for tag in ("a:noFill", "a:solidFill", "a:gradFill", "a:blipFill", "a:pattFill", "a:grpFill"):
        for el in sp.findall(qn(tag)):
            sp.remove(el)
    patt = sp.makeelement(qn("a:pattFill"), {"prst": prst})
    for tag, clr in (("a:fgClr", fg), ("a:bgClr", bg)):
        e = patt.makeelement(qn(tag), {})
        e.append(e.makeelement(qn("a:srgbClr"), {"val": "%02X%02X%02X" % (clr[0], clr[1], clr[2])}))
        patt.append(e)
    ln = sp.find(qn("a:ln"))
    sp.insert(list(sp).index(ln) if ln is not None else len(sp), patt)


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
    txt = f"{n_sig}/{n_tested}/{n_genome}"
    if n_sig == 0:
        return "ns", RGBColor(*NS), txt
    base = UP if n_up >= n_down else DOWN
    return ("up" if base is UP else "down"), _blend(base, n_sig / n_tested), txt


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

    _text(slide, 0.4, 0.10, SLIDE_W - 0.8, 0.34,
          "41 candidate starvation-sensitivity COGs vs. nutrient-starvation experiments",
          16, bold=True, align=PP_ALIGN.CENTER)
    _text(slide, 0.4, 0.46, SLIDE_W - 0.8, 0.44,
          "cell text k / n / m = significant / measured in this experiment / copies in that strain's genome  ·  "
          "fill = dominant direction, blended toward grey by k/n  ·  columns grouped by functional category, "
          "N / mixed strip = that COG's Direction tag  ·  * Martiny MIT9313 at 24h (no 48h row); table-absent COGs = not significant",
          9, color=MUTED, align=PP_ALIGN.CENTER)

    # group brackets + category labels, with a per-COG Direction strip beneath
    for a, b, cat in groups:
        x0 = GRID_LEFT + a * cw
        x1 = GRID_LEFT + (b + 1) * cw
        bar = slide.shapes.add_shape(MSO_SHAPE.RECTANGLE, Inches(x0 + 0.03), Inches(GRID_TOP - 0.52),
                                     Inches(x1 - x0 - 0.06), Inches(0.022))
        bar.fill.solid(); bar.fill.fore_color.rgb = INK
        bar.line.fill.background(); bar.shadow.inherit = False
        _text(slide, x0 - 0.4, GRID_TOP - 1.10, (x1 - x0) + 0.8, 0.5,
              CATEGORY_DISPLAY.get(cat, cat), 8, color=INK, align=PP_ALIGN.CENTER,
              anchor=MSO_ANCHOR.BOTTOM)

    # Direction strip: N / mixed per COG, black, just under the brackets
    for j, d in enumerate(col_dir):
        _text(slide, GRID_LEFT + j * cw, GRID_TOP - 0.32, cw, 0.24,
              "N" if d == "N" else "mixed", 7, color=INK, align=PP_ALIGN.CENTER)

    # column labels (rotated), all black
    for j, lab in enumerate(col_labels):
        cx = GRID_LEFT + j * cw + cw / 2
        _text(slide, cx - 1.0, GRID_TOP + len(ROW_ORDER) * CELL_H + 0.15, 2.0, 0.24, lab, 8,
              color=INK, align=PP_ALIGN.LEFT, rot=300, font="Consolas")

    # rows: label + 41 cells
    for i, exp in enumerate(ROW_ORDER):
        top = GRID_TOP + i * CELL_H
        _text(slide, 0.15, top, LEFT_LABELS - 0.25, CELL_H, EXPERIMENT_LABELS[exp], 9.5,
              color=SECONDARY, align=PP_ALIGN.RIGHT)
        for j, cog in enumerate(cog_order):
            cat, fill, txt = specs.get((exp, cog), ("absent", None, ""))
            left = GRID_LEFT + j * cw
            shp = slide.shapes.add_shape(MSO_SHAPE.RECTANGLE, Inches(left), Inches(top), Inches(cw), Inches(CELL_H))
            shp.name = f"cell|{EXPERIMENT_LABELS[exp]}|{cog}|{cat}"
            shp.shadow.inherit = False
            shp.line.color.rgb = WHITE
            shp.line.width = Pt(1.0)
            if fill is None:
                _pattern_fill(shp, "ltUpDiag", ABSENT_FG, ABSENT_BG)
            else:
                shp.fill.solid(); shp.fill.fore_color.rgb = fill
            if txt:
                _text(slide, left, top, cw, CELL_H, txt, 6, color=INK,
                      align=PP_ALIGN.CENTER, anchor=MSO_ANCHOR.MIDDLE)

    grid_bottom = GRID_TOP + len(ROW_ORDER) * CELL_H
    grid_right = GRID_LEFT + len(cog_order) * cw

    n_rows = sum(1 for e in ROW_ORDER if EXPERIMENT_LABELS[e].startswith("N:"))
    hline = slide.shapes.add_shape(MSO_SHAPE.RECTANGLE, Inches(GRID_LEFT),
                                   Inches(GRID_TOP + n_rows * CELL_H - 0.012),
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
              (RGBColor(*NS), "Tested, no gene significant"),
              (RGBColor(*NODATA), "No data (COG present, not in table)"),
              (None, "COG absent from strain genome")]
    lx = GRID_LEFT
    for clr, lab in legend:
        sw = slide.shapes.add_shape(MSO_SHAPE.RECTANGLE, Inches(lx), Inches(ly), Inches(0.28), Inches(0.28))
        sw.shadow.inherit = False; sw.line.color.rgb = MUTED; sw.line.width = Pt(0.75)
        if clr is None:
            _pattern_fill(sw, "ltUpDiag", ABSENT_FG, ABSENT_BG)
        else:
            sw.fill.solid(); sw.fill.fore_color.rgb = clr
        _text(slide, lx + 0.36, ly - 0.02, 3.4, 0.32, lab, 9.5)
        lx += 0.36 + 0.058 * len(lab) + 0.6

    out = FIG_DIR / "01_gene_experiment_heatmap.pptx"
    prs.save(out)
    print(f"Wrote {out}  ({len(ROW_ORDER)} rows x {len(cog_order)} COGs = "
          f"{len(ROW_ORDER) * len(cog_order)} editable cells)")


if __name__ == "__main__":
    main()
