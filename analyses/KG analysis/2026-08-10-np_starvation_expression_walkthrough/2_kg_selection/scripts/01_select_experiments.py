"""
Step 2 -- select N/P starvation experiments in scope.

Purpose: pull every Prochlorococcus nitrogen- and phosphorus-treatment
experiment from the KG, then apply the step-1 locked scope:
  - organism (strain) must be one of the researcher's 15 study strains
  - background factors must be axenic and not phage-infected (no coculture,
    no viral background)
  - must be a literal starvation/limitation contrast (not an
    alternate-N-source substrate switch, not a chemical-inhibitor proxy)
  - must carry gene-level differential expression (excludes
    metabolomics-only experiments)

Inputs: none (queries the KG directly).
Outputs: data/01_np_experiments.csv -- one row per experiment, with the
  include/exclude decision and the reason.

Usage: uv run analyses/2026-08-10-np_starvation_expression_walkthrough/2_kg_selection/scripts/01_select_experiments.py
"""

from pathlib import Path

import pandas as pd
from multiomics_explorer import list_experiments, GraphConnection

OUT_DIR = Path(__file__).resolve().parent.parent / "data"

STUDY_STRAINS = {
    "Prochlorococcus MED4", "Prochlorococcus MIT9312", "Prochlorococcus MIT9313",
    "Prochlorococcus MIT1327", "Prochlorococcus MIT0604", "Prochlorococcus NATL2A",
    "Prochlorococcus MIT9515", "Prochlorococcus MIT9215", "Prochlorococcus AS9601",
    "Prochlorococcus PAC1", "Prochlorococcus MIT9202", "Prochlorococcus SB",
    "Prochlorococcus MIT9301", "Prochlorococcus MIT1314", "Prochlorococcus NATL1A",
}

# Substring markers for alternate-N-source substrate-switch experiments
# (growth on cyanate/urea/nitrite as sole N source vs. replete -- not starvation).
ALT_N_SOURCE_MARKERS = ["as sole n source", "azaserine"]


def classify(row: dict) -> tuple[bool, str]:
    """Return (include, reason) for one experiment row."""
    organism = row["organism_name"]
    treatment = (row.get("treatment") or "").lower()
    experiment_name = (row.get("experiment_name") or "").lower()
    background = set(row.get("background_factors") or [])
    omics = row["omics_type"]

    if organism not in STUDY_STRAINS:
        return False, f"strain '{organism}' not in the researcher's 15 study strains"
    if omics == "METABOLOMICS":
        return False, "metabolomics-only: no gene-level differential expression"
    if "coculture" in background:
        return False, "coculture background factor (confound excluded per step 1)"
    if "viral" in background:
        return False, "phage-infection background factor (confound excluded per step 1)"
    if any(marker in treatment or marker in experiment_name for marker in ALT_N_SOURCE_MARKERS):
        return False, "alternate-N-source / chemical-proxy, not literal starvation (excluded per step 1)"
    return True, "axenic, uninfected, literal starvation/limitation vs. replete, in-scope strain"


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    rows = []
    with GraphConnection() as conn:
        for nutrient in ("nitrogen", "phosphorus"):
            result = list_experiments(
                organism="Prochlorococcus",
                treatment_type=[nutrient],
                verbose=True,
                conn=conn,
            )
            print(f"{nutrient}: {result['total_matching']} experiments returned by KG")
            for r in result["results"]:
                include, reason = classify(r)
                rows.append({
                    "nutrient": nutrient,
                    "experiment_id": r["experiment_id"],
                    "experiment_name": r["experiment_name"],
                    "organism_name": r["organism_name"],
                    "publication_doi": r["publication_doi"],
                    "publication_title": r.get("publication_title"),
                    "omics_type": r["omics_type"],
                    "background_factors": "|".join(r.get("background_factors") or []),
                    "treatment": r.get("treatment"),
                    "control": r.get("control"),
                    "table_scope": r.get("table_scope"),
                    "table_scope_detail": r.get("table_scope_detail"),
                    "is_time_course": r.get("is_time_course"),
                    "distinct_gene_count": r.get("distinct_gene_count"),
                    "significant_up": r.get("genes_by_status", {}).get("significant_up"),
                    "significant_down": r.get("genes_by_status", {}).get("significant_down"),
                    "include": include,
                    "reason": reason,
                })

    df = pd.DataFrame(rows)
    out_path = OUT_DIR / "01_np_experiments.csv"
    df.to_csv(out_path, index=False)

    print(f"\nWrote {len(df)} rows to {out_path}")
    print("\nFilter funnel:")
    print(df.groupby("nutrient")["include"].agg(["sum", "count"]))
    print("\nIncluded experiments:")
    included = df[df["include"]]
    print(included[["nutrient", "experiment_id", "organism_name", "omics_type"]].to_string(index=False))
    print("\nExcluded experiments (reason):")
    excluded = df[~df["include"]]
    print(excluded[["nutrient", "experiment_id", "organism_name", "reason"]].to_string(index=False))
    print("\nIncluded strains by nutrient:")
    print(included.groupby("nutrient")["organism_name"].unique())


if __name__ == "__main__":
    main()
