"""
Step 5 QC -- enumerate the (COG x experiment) heatmap cells where the COG's
resolved genes disagree in direction (>=1 significant_up AND >=1
significant_down in the same cell).

Figure 1 colours each such cell by the dominant direction only, so the
disagreement is invisible unless the cell text spells it out. This script
lists every split cell with its genes, so the figure's "↑up↓down" k-text
can be checked against the underlying rows.

Input:  data/01_target_gene_experiment_matrix.csv
        ../../2_kg_selection/data/00_source_normixed_cogs.csv
Output: data/03_within_cell_direction_splits.csv

Usage:
  uv run analyses/2026-09-06-cog_starvation_sensitivity_validation/5_analyze/scripts/05_within_cell_direction_splits.py
"""

from pathlib import Path

import pandas as pd

BASE = Path(__file__).resolve().parent.parent
DATA_DIR = BASE / "data"
STEP2_DATA = BASE.parent / "2_kg_selection" / "data"

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
TESTED = {"significant_up", "significant_down", "not_significant"}


def main() -> None:
    m = pd.read_csv(DATA_DIR / "01_target_gene_experiment_matrix.csv")
    cogs = pd.read_csv(STEP2_DATA / "00_source_normixed_cogs.csv")
    name = dict(zip(cogs["COGs"], cogs["Protien"]))

    rows = []
    for (exp, cog), g in m.groupby(["experiment_id", "cog_number"]):
        st = g["status"].value_counts().to_dict()
        n_up = st.get("significant_up", 0)
        n_down = st.get("significant_down", 0)
        if n_up == 0 or n_down == 0:
            continue
        n_tested = sum(st.get(s, 0) for s in TESTED)
        n_genome = int(g["locus_tag"].notna().sum())
        sig = g[g["status"].isin(["significant_up", "significant_down"])]
        rows.append({
            "experiment": EXPERIMENT_LABELS.get(exp, exp),
            "experiment_id": exp,
            "cog_number": cog,
            "protein": name.get(cog),
            "direction_tag": g["direction"].iloc[0],
            "n_up": n_up,
            "n_down": n_down,
            "k_text": f"↑{n_up}↓{n_down}",
            "n_measured": n_tested,
            "n_genome_copies": n_genome,
            "dominant_colour": "up" if n_up >= n_down else "down",
            "up_loci": "; ".join(
                f"{r.locus_tag}({r.log2fc:+.2f})"
                for r in sig[sig.status == "significant_up"].itertuples()
            ),
            "down_loci": "; ".join(
                f"{r.locus_tag}({r.log2fc:+.2f})"
                for r in sig[sig.status == "significant_down"].itertuples()
            ),
        })

    out = pd.DataFrame(rows).sort_values(["experiment", "cog_number"]).reset_index(drop=True)
    out.to_csv(DATA_DIR / "03_within_cell_direction_splits.csv", index=False)
    print(f"Wrote data/03_within_cell_direction_splits.csv  ({len(out)} split cells)")
    print(out[["experiment", "cog_number", "protein", "direction_tag",
               "n_up", "n_down", "n_measured", "n_genome_copies", "dominant_colour"]].to_string(index=False))


if __name__ == "__main__":
    main()
