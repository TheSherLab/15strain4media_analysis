"""
Step 3 -- document each in-scope experiment's actual significance criteria.

Purpose: the researcher asked, before any analysis runs, to know exactly
what p-value/threshold each experiment used to call a gene "significant",
and whether the source table reports every tested gene or only a filtered
subset (e.g. "only up/down regulated", "top 50%", "q<0.05 at one
timepoint"). table_scope / table_scope_detail / statistical_test (from
list_experiments) answer the "what's reported" half; this script adds the
empirical half -- pulling every DE row per experiment and finding the
padj/log2FC boundary that actually separates "significant" from
"not_significant" in the data, since expression_status uses each
publication's own threshold, not a uniform padj<0.05.

Inputs: ../../2_kg_selection/data/01_np_experiments.csv (the 10 in-scope
  experiments).
Outputs: data/01_significance_criteria.csv -- one row per experiment.

Usage: uv run analyses/2026-08-10-np_starvation_expression_walkthrough/3_analysis_framing/scripts/01_significance_criteria.py
"""

from pathlib import Path

import pandas as pd
from multiomics_explorer import differential_expression_by_gene, GraphConnection

BASE = Path(__file__).resolve().parent.parent
DATA_DIR = BASE / "data"
EXPERIMENTS_CSV = BASE.parent / "2_kg_selection" / "data" / "01_np_experiments.csv"


def main() -> None:
    experiments = pd.read_csv(EXPERIMENTS_CSV)
    included = experiments[experiments["include"]].copy()
    included["organism_short"] = included["organism_name"].str.replace("Prochlorococcus ", "", regex=False)

    rows = []
    with GraphConnection() as conn:
        for _, exp in included.iterrows():
            result = differential_expression_by_gene(
                organism=exp["organism_short"],
                experiment_ids=[exp["experiment_id"]],
                limit=None,
                conn=conn,
            )
            de = pd.DataFrame(result["results"])
            sig = de[de["expression_status"].isin(["significant_up", "significant_down"])]
            not_sig = de[de["expression_status"] == "not_significant"]

            row = {
                "nutrient": exp["nutrient"],
                "experiment_id": exp["experiment_id"],
                "organism": exp["organism_short"],
                "omics_type": exp["omics_type"],
                "statistical_test": None,  # filled from experiments csv below
                "table_scope": exp["table_scope"],
                "table_scope_detail": exp["table_scope_detail"],
                "n_rows_total": len(de),
                "n_significant": len(sig),
                "n_not_significant": len(not_sig),
                "has_not_significant_rows": len(not_sig) > 0,
                "max_padj_among_significant": sig["padj"].max() if len(sig) else None,
                "min_padj_among_not_significant": not_sig["padj"].min() if len(not_sig) else None,
                "min_abs_log2fc_among_significant": sig["log2fc"].abs().min() if len(sig) else None,
            }
            rows.append(row)
            print(f"{exp['experiment_id']}: {len(de)} rows, {len(sig)} significant, "
                  f"{len(not_sig)} not_significant, table_scope={exp['table_scope']}")

    df = pd.DataFrame(rows)
    # statistical_test isn't in the DE-row envelope; pull from list_experiments verbose data
    # already captured in the researcher-facing notebook table -- merge by experiment_id.
    stat_test_map = {
        "10.1101/2025.11.24.690089_growth_state_pro99lown_nutrient_starvation_med4_proteomics_axenic": "DESeq2",
        "10.1101/2025.11.24.690089_growth_state_pro99lown_nutrient_starvation_med4_rnaseq_axenic": "DESeq2",
        "10.1038/ismej.2017.88_nitrogen_stress_ndepleted_pro99_medium_med4_rnaseq": "Rockhopper",
        "10.1038/msb4100087_nitrogen_nitrogen_deprivation_med4_med4_microarray": "microarray_Goldenspike",
        "10.1038/msb4100087_nitrogen_nitrogen_deprivation_mit9313_mit9313_microarray": "microarray_Goldenspike",
        "10.1111/1462-2920.13104_phosphorus_plimited_natl2a_rnaseq_uninfected": "DESeq2",
        "10.1186/2046-9063-8-7_pi_limitation_mit9312_itraq": "iTRAQ_t-test",
        "10.1186/2046-9063-8-7_pi_limitation_natl2a_itraq": "iTRAQ_t-test",
        "10.1073/pnas.0601301103_phosphorus_phosphate_starvation_med4_microarray": "microarray_Cyber-T",
        "10.1073/pnas.0601301103_phosphorus_phosphate_starvation_mit9313_microarray": "microarray_Cyber-T",
    }
    df["statistical_test"] = df["experiment_id"].map(stat_test_map)

    out_path = DATA_DIR / "01_significance_criteria.csv"
    df.to_csv(out_path, index=False)
    print(f"\nWrote {len(df)} rows to {out_path}")
    print(df[["nutrient", "organism", "omics_type", "statistical_test", "table_scope",
              "has_not_significant_rows", "max_padj_among_significant",
              "min_padj_among_not_significant"]].to_string(index=False))


if __name__ == "__main__":
    main()
