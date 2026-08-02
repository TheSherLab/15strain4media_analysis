"""Enumerate the axenic Prochlorococcus nitrogen and phosphorus starvation
experiments in the KG, per the scope locked in step 1 and refined during
step 2's co-define:

- Axenic only (no coculture-with-Alteromonas experiments).
- Nitrogen: strict starvation/deprivation contrasts only -- excludes the
  4 "growth on alternate N source (cyanate/urea/nitrite) vs N-replete"
  experiments, which test a different biological question (alternate
  N-source utilization, not starvation), per researcher decision.
- Phosphorus: all axenic P-limitation/deplete-vs-replete contrasts
  (no equivalent "alternate P source" complication exists in the KG for
  phosphorus treatments).
- SS120 (CCMP1375) excluded entirely -- not part of the researcher's
  study (step 3 co-define, 2026-08-02 redo).
- Gene-expression/proteomics data only (RNASEQ, PROTEOMICS, MICROARRAY)
  -- METABOLOMICS experiments are excluded: they measure metabolite
  abundance, not gene expression, and differential_expression_by_gene
  returns zero rows for them (discovered 2026-08-02 while building the
  step-3 background gene set; this drops MIT9301 entirely, since both of
  its phosphorus experiments are metabolomics-only).

Output: data/02_np_experiments.csv -- one row per experiment, with
organism, omics type, publication, treatment/control text.

Usage:
  uv run python 02_list_np_experiments.py
"""
from pathlib import Path

import pandas as pd
from multiomics_explorer import list_experiments

DATA_DIR = Path(__file__).resolve().parents[1] / "data"

# Experiment IDs corresponding to "growth on alternate N source vs
# N-replete" -- excluded per the researcher's strict-starvation-only
# decision (step 2 co-define, 2026-08-02).
ALTERNATE_N_SOURCE_EXCLUDED = {
    "10.1038/msb4100087_growth_medium_growth_on_cyanate_as_med4_microarray",
    "10.1038/msb4100087_growth_medium_growth_on_urea_as_med4_microarray",
    "10.1038/msb4100087_growth_medium_growth_on_nitrite_as_mit9313_microarray",
    "10.1038/msb4100087_growth_medium_growth_on_urea_as_mit9313_microarray",
}


def fetch(treatment_type: str) -> pd.DataFrame:
    result = list_experiments(
        treatment_type=[treatment_type],
        organism="Prochlorococcus",
        background_factors=["axenic"],
        verbose=True,
        limit=None,
    )
    rows = []
    for e in result["results"]:
        if e["experiment_id"] in ALTERNATE_N_SOURCE_EXCLUDED:
            continue
        if "SS120" in e["organism_name"] or "CCMP1375" in e["organism_name"]:
            continue
        if e["omics_type"] not in {"RNASEQ", "PROTEOMICS", "MICROARRAY"}:
            continue
        rows.append(
            {
                "n_or_p": "N" if treatment_type == "nitrogen" else "P",
                "experiment_id": e["experiment_id"],
                "organism_name": e["organism_name"],
                "omics_type": e["omics_type"],
                "background_factors": "|".join(e["background_factors"]),
                "treatment": e["treatment"],
                "control": e["control"],
                "publication_doi": e["publication_doi"],
                "publication_title": e["publication_title"],
                "is_time_course": e["is_time_course"],
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    n_df = fetch("nitrogen")
    p_df = fetch("phosphorus")
    out = pd.concat([n_df, p_df], ignore_index=True)

    out_path = DATA_DIR / "02_np_experiments.csv"
    out.to_csv(out_path, index=False)
    print(f"Wrote {len(out)} experiments to {out_path}")
    print()
    print(out.groupby(["n_or_p", "organism_name"]).size().to_string())
    print()
    print("Strains involved:", sorted(out["organism_name"].unique()))
    print()
    print("Publications:")
    print(out[["publication_doi", "publication_title"]].drop_duplicates().to_string(index=False))


if __name__ == "__main__":
    main()
