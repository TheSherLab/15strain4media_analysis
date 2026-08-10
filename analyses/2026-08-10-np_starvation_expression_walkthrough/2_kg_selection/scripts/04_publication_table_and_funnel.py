"""
Step 2 -- publication funnel figure + publications table.

Purpose: the researcher asked to see how the KG's publication set narrows
down to the 6 publications actually used, and a table describing each of
those 6 (and the 2 that had N/P-treatment experiments but were fully
excluded by step-1 scope).

Inputs: KG queries (list_publications) -- see notebook.md for the exact
  calls and their counts. Funnel counts are hardcoded here because they
  come from 4 distinct queries run interactively (KG-total,
  Prochlorococcus-subset, N/P-treatment-subset) plus the step-1 scope
  filter already computed in 01_np_experiments.csv; recomputing the first
  three from scratch on every run would re-issue the same 3 KG calls for
  numbers that don't depend on anything downstream.
Outputs:
  data/04_publications_table.csv -- one row per publication (6 included +
    2 excluded), with authors, journal, what it tested, data type.
  figures/01_publication_funnel.png -- funnel bar chart.

Usage: uv run analyses/2026-08-10-np_starvation_expression_walkthrough/2_kg_selection/scripts/04_publication_table_and_funnel.py
"""

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

BASE = Path(__file__).resolve().parent.parent
DATA_DIR = BASE / "data"
FIG_DIR = BASE / "figures"

# Funnel stages, from list_publications queries run 2026-08-10 (see notebook.md).
FUNNEL = [
    ("All publications in KG", 43),
    ("...involving Prochlorococcus", 33),
    ("...with a nitrogen- or phosphorus-\ntreatment experiment", 8),
    ("...with >=1 experiment surviving\nstep-1 scope (final: 6)", 6),
]

# Ordinal sequential ramp (dataviz palette.md), light mode, funnel-stage use:
# step 250 -> 350 -> 450 -> 550, darkening as the set narrows.
STAGE_COLORS = ["#86b6ef", "#5598e7", "#2a78d6", "#1c5cab"]

INK = "#0b0b0b"
SECONDARY_INK = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
SURFACE = "#fcfcfb"

PUBLICATIONS = [
    # (doi, authors, journal, year, nutrient(s), what_it_tested, data_type, strains, status, reason_if_excluded)
    dict(doi="10.1101/2025.11.24.690089", authors="Weissberg, Aharonovich, Sher", journal="bioRxiv", year=2025,
         nutrient="N", what_it_tested="MED4 nutrient-starvation vs. exponential growth, axenic and in coculture with Alteromonas (90-day time course)",
         data_type="RNA-seq, Proteomics", strains="MED4", status="included",
         reason="axenic arm in scope; coculture arm excluded (confound)"),
    dict(doi="10.1038/ismej.2017.88", authors="Read, Berube, Biller, Neveux, Cubillos-Ruiz, Chisholm, Grzymski", journal="The ISME Journal", year=2017,
         nutrient="N", what_it_tested="MED4 N-depleted vs. N-replete Pro99 medium, transcription-start-site mapping",
         data_type="RNA-seq", strains="MED4", status="included", reason=""),
    dict(doi="10.1038/msb4100087", authors="Tolonen, Aach, Lindell, Johnson, Rector, Steen, Church, Chisholm", journal="Molecular Systems Biology", year=2006,
         nutrient="N", what_it_tested="MED4 and MIT9313 nitrogen deprivation time course vs. N-replete Pro99 (plus growth on alternate N sources, excluded here)",
         data_type="Microarray", strains="MED4, MIT9313", status="included",
         reason="deprivation experiments in scope; alt-N-source experiments excluded (substrate switch, not starvation)"),
    dict(doi="10.1111/1462-2920.13104", authors="Lin, Ding, Zeng", journal="Environmental Microbiology", year=2015,
         nutrient="P", what_it_tested="NATL2A P-limited vs. P-replete, with and without cyanophage P-SSM2 infection",
         data_type="RNA-seq", strains="NATL2A", status="included",
         reason="uninfected arm in scope; infected arm excluded (confound)"),
    dict(doi="10.1186/2046-9063-8-7", authors="Fuszard, Wright, Biggs", journal="Aquatic Biosystems", year=2012,
         nutrient="P", what_it_tested="MIT9312, NATL2A, SS120 phosphate-replete vs. phosphate-deplete (iTRAQ)",
         data_type="Proteomics", strains="MIT9312, NATL2A", status="included",
         reason="MIT9312/NATL2A in scope; SS120 excluded (not one of the researcher's 15 study strains)"),
    dict(doi="10.1073/pnas.0601301103", authors="Martiny, Coleman, Chisholm", journal="PNAS", year=2006,
         nutrient="P", what_it_tested="MED4 and MIT9313 phosphate starvation vs. P-replete Pro99, time course",
         data_type="Microarray", strains="MED4, MIT9313", status="included", reason=""),
    dict(doi="10.1128/mSystems.00008-17", authors="Dominguez-Martin, Gomez-Baena, Diez, Lopez-Grueso, Beynon, Garcia-Fernandez", journal="mSystems", year=2017,
         nutrient="N", what_it_tested="SS120 azaserine (GOGAT inhibitor) vs. untreated, as a chemical proxy for N-limitation",
         data_type="Proteomics", strains="SS120", status="excluded",
         reason="strain not in the researcher's 15 study strains; also a chemical-proxy design, not literal media starvation"),
    dict(doi="10.1128/msystems.01261-22", authors="Kujawinski, Braakman, Longnecker, Becker, Chisholm, Dooley, Kido Soule, Swarr, Halloran", journal="mSystems", year=2023,
         nutrient="P", what_it_tested="MIT9301 metabolite profiling, replete vs. P-limited (intra- and extracellular)",
         data_type="Metabolomics", strains="MIT9301", status="excluded",
         reason="metabolomics only -- no gene-level differential expression"),
]


def make_table() -> pd.DataFrame:
    df = pd.DataFrame(PUBLICATIONS)
    out_path = DATA_DIR / "04_publications_table.csv"
    df.to_csv(out_path, index=False)
    print(f"Wrote {len(df)} rows to {out_path}")
    return df


def make_funnel_figure() -> None:
    fig, ax = plt.subplots(figsize=(8.5, 4.2), dpi=300)
    fig.patch.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)

    labels = [s[0] for s in FUNNEL]
    values = [s[1] for s in FUNNEL]
    y_pos = range(len(FUNNEL))[::-1]  # first stage on top

    ax.barh(list(y_pos), values, color=STAGE_COLORS, height=0.6, zorder=3)

    for y, v in zip(y_pos, values):
        ax.text(v + 0.8, y, str(v), va="center", ha="left",
                color=INK, fontsize=12, fontweight="bold")

    ax.set_yticks(list(y_pos))
    ax.set_yticklabels(labels, color=SECONDARY_INK, fontsize=10)
    ax.set_xlabel("Number of publications", color=MUTED, fontsize=10)
    ax.set_xlim(0, 48)
    ax.tick_params(axis="x", colors=MUTED, labelsize=9)
    ax.tick_params(axis="y", length=0)

    ax.grid(axis="x", color=GRID, linewidth=0.8, zorder=0)
    for spine in ("top", "right", "left"):
        ax.spines[spine].set_visible(False)
    ax.spines["bottom"].set_color("#c3c2b7")

    ax.set_title(
        "From the KG's publication set to this analysis's 6 publications",
        color=INK, fontsize=11, fontweight="bold", loc="left", pad=12, wrap=True,
    )

    fig.tight_layout()
    out_path = FIG_DIR / "01_publication_funnel.png"
    fig.savefig(out_path, dpi=300, facecolor=SURFACE)
    print(f"Wrote {out_path}")


if __name__ == "__main__":
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    make_table()
    make_funnel_figure()
