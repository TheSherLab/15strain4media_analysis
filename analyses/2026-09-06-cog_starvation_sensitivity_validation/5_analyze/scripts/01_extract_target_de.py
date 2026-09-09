"""
Step 5 (redo 2026-09-09) -- extract DE data for every gene carrying one of
the 41 COGs, across all 10 in-scope experiments (5 nitrogen + 5
phosphorus), each at its single chosen starvation timepoint.

Changes vs the 2026-09-06 version:
  - gene set = the full per-strain COG membership from step 2's rebuilt
    02_gene_locus_resolution.csv (was one Cyanorak ID per COG)
  - phosphorus back in scope (researcher-directed 2026-09-09)
  - timepoints from the walkthrough analysis's final table -- Lin at 59h,
    Martiny at 48h / 24h.
  - "gene present in strain but absent from the source table" reclassified
    from "no data" to "not_significant" (RECLASSIFY_ABSENT_AS_NOT_SIG) for
    the two Martiny genome-wide P microarrays AND for the Lin uninfected
    RNA-seq (2026-09-09, researcher-directed): the KG marks the Lin
    uninfected experiment `table_scope = significant_only` and it is a
    genome-wide RNA-seq, so a COG gene measured but not in Lin's reported
    set was tested and not significant -- same logic as Martiny. This
    overrides the walkthrough's 2026-09-08 "not applied to Lin" note.
    Still NOT applied to Fuszard (iTRAQ fixed detected-protein list).

One row per (COG, locus, experiment). A COG resolves to several loci per
strain; every locus is a row (matches how the copy-number test counted
them). States: significant_up / significant_down / not_significant /
no_data_at_timepoint / no_locus_in_strain.

Inputs: ../../2_kg_selection/data/01_np_experiments.csv,
        ../../2_kg_selection/data/02_gene_locus_resolution.csv
Outputs: data/01_target_gene_experiment_matrix.csv

Usage: uv run analyses/2026-09-06-cog_starvation_sensitivity_validation/5_analyze/scripts/01_extract_target_de.py
"""

from pathlib import Path

import pandas as pd
from multiomics_explorer import differential_expression_by_gene, GraphConnection

BASE = Path(__file__).resolve().parent.parent
DATA_DIR = BASE / "data"
STEP2_DATA = BASE.parent / "2_kg_selection" / "data"

# Chosen starvation timepoint per experiment -- reused verbatim from the
# walkthrough analysis 5_analyze/notebook.md (its final 2026-09-08 table).
# None = single-timepoint experiment.
STARVATION_TIMEPOINT = {
    "10.1101/2025.11.24.690089_growth_state_pro99lown_nutrient_starvation_med4_proteomics_axenic": "day 14",
    "10.1101/2025.11.24.690089_growth_state_pro99lown_nutrient_starvation_med4_rnaseq_axenic": None,
    "10.1038/ismej.2017.88_nitrogen_stress_ndepleted_pro99_medium_med4_rnaseq": "24h",
    "10.1038/msb4100087_nitrogen_nitrogen_deprivation_med4_med4_microarray": "12h",
    "10.1038/msb4100087_nitrogen_nitrogen_deprivation_mit9313_mit9313_microarray": "12h",
    "10.1111/1462-2920.13104_phosphorus_plimited_natl2a_rnaseq_uninfected": "59h",
    "10.1186/2046-9063-8-7_pi_limitation_mit9312_itraq": None,
    "10.1186/2046-9063-8-7_pi_limitation_natl2a_itraq": None,
    "10.1073/pnas.0601301103_phosphorus_phosphate_starvation_med4_microarray": "48h",
    "10.1073/pnas.0601301103_phosphorus_phosphate_starvation_mit9313_microarray": "24h",  # flagged, see notebook
}

# Genome-wide assays reported as a significance-filtered gene list
# (table_scope = filtered_subset / significant_only): a gene present in the
# strain but absent from the table was measured and did not pass.
RECLASSIFY_ABSENT_AS_NOT_SIG = {
    "10.1073/pnas.0601301103_phosphorus_phosphate_starvation_med4_microarray",
    "10.1073/pnas.0601301103_phosphorus_phosphate_starvation_mit9313_microarray",
    "10.1111/1462-2920.13104_phosphorus_plimited_natl2a_rnaseq_uninfected",  # Lin, significant_only (2026-09-09)
}


def main() -> None:
    experiments = pd.read_csv(STEP2_DATA / "01_np_experiments.csv")
    experiments = experiments[experiments["include"]].copy()
    assert len(experiments) == 10, len(experiments)
    experiments["organism_short"] = experiments["organism_name"].str.replace("Prochlorococcus ", "", regex=False)

    res = pd.read_csv(STEP2_DATA / "02_gene_locus_resolution.csv")
    res = res[res["status"] != "unmapped_in_kg"]
    cog_meta = res[["cog_number", "direction"]].drop_duplicates().set_index("cog_number")["direction"].to_dict()
    all_cogs = sorted(cog_meta)

    rows = []
    with GraphConnection() as conn:
        for _, exp in experiments.iterrows():
            exp_id = exp["experiment_id"]
            organism_full = exp["organism_name"]
            organism_short = exp["organism_short"]
            nutrient = exp["nutrient"]
            chosen_tp = STARVATION_TIMEPOINT[exp_id]

            strain_res = res[res["organism_name"] == organism_full]
            # (cog, locus) pairs present in this strain
            triples = list(strain_res[["cog_number", "locus_tag"]].itertuples(index=False, name=None))
            distinct_loci = sorted({lt for _, lt in triples})

            de_by_locus = {}
            if distinct_loci:
                result = differential_expression_by_gene(
                    organism=organism_short, locus_tags=distinct_loci,
                    experiment_ids=[exp_id], limit=None, conn=conn,
                )
                de = pd.DataFrame(result["results"])
                if chosen_tp is not None and len(de):
                    de = de[de["timepoint"] == chosen_tp]
                if len(de):
                    # a locus can have >1 row (shouldn't at one timepoint) -- keep the most significant
                    de = de.sort_values("padj").drop_duplicates("locus_tag", keep="first")
                    de_by_locus = {r["locus_tag"]: r for _, r in de.iterrows()}

            cogs_here = set(strain_res["cog_number"])
            for cog, locus in triples:
                direction = cog_meta[cog]
                if locus in de_by_locus:
                    r = de_by_locus[locus]
                    status, log2fc, padj = r["expression_status"], r["log2fc"], r["padj"]
                else:
                    status = ("not_significant" if exp_id in RECLASSIFY_ABSENT_AS_NOT_SIG
                              else "no_data_at_timepoint")
                    log2fc = padj = None
                rows.append({
                    "cog_number": cog, "direction": direction, "nutrient": nutrient,
                    "experiment_id": exp_id, "organism_name": organism_full,
                    "locus_tag": locus, "status": status, "log2fc": log2fc, "padj": padj,
                    "timepoint": chosen_tp,
                })

            for cog in all_cogs:
                if cog not in cogs_here:
                    rows.append({
                        "cog_number": cog, "direction": cog_meta[cog], "nutrient": nutrient,
                        "experiment_id": exp_id, "organism_name": organism_full,
                        "locus_tag": None, "status": "no_locus_in_strain",
                        "log2fc": None, "padj": None, "timepoint": chosen_tp,
                    })

            print(f"{exp_id} ({organism_short}, tp={chosen_tp}): {len(triples)} (cog,locus) pairs, "
                  f"{len(distinct_loci)} loci queried")

    df = pd.DataFrame(rows)
    out = DATA_DIR / "01_target_gene_experiment_matrix.csv"
    df.to_csv(out, index=False)
    print(f"\nWrote {len(df)} rows to {out}")
    print("\nstatus x nutrient:")
    print(pd.crosstab(df["status"], df["nutrient"]))


if __name__ == "__main__":
    main()
