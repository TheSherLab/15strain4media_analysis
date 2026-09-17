"""
Step 4 -- driving example, real KG data. Reuses hit_rate() and
bootstrap_pvalue() unchanged from the prior analysis's
4_methods/np_response.py (already toy-verified there against
hand-computed synthetic data -- not re-verified here, since the functions
themselves are untouched; only their application to this new gene set is
new).

Driving example: COG1403 (McrA, N-tagged, HNH endonuclease family), the
COG with the most DE evidence (28 rows) of the 41. Pulls every DE row
across MED4 (2 in-genome-duplicate loci, PMM0031 and PMM1403) and MIT9313
in the 5 nitrogen experiments, hand-tallies significant/direction counts
independently of hit_rate()'s own computation, and confirms they match.

Inputs: ../../2_kg_selection/data/01_np_experiments.csv,
        ../../2_kg_selection/data/02_gene_locus_resolution.csv
Outputs: data/01_driving_example_mcra.csv

Usage: uv run analyses/2026-09-06-cog_starvation_sensitivity_validation/4_methods/scripts/01_driving_example.py
"""

import sys
from pathlib import Path

import pandas as pd
from multiomics_explorer import differential_expression_by_gene, GraphConnection

BASE = Path(__file__).resolve().parent.parent
DATA_DIR = BASE / "data"
STEP2_DATA = BASE.parent / "2_kg_selection" / "data"
PRIOR_METHODS_DIR = (
    BASE.parent.parent / "2026-08-10-np_starvation_expression_walkthrough" / "4_methods"
)

sys.path.insert(0, str(PRIOR_METHODS_DIR))
from np_response import hit_rate  # noqa: E402  (reused unchanged from the prior analysis)

COG = "COG1403"


def main() -> None:
    experiments = pd.read_csv(STEP2_DATA / "01_np_experiments.csv")
    # nitrogen only: 01_np_experiments.csv still carries the phosphorus rows
    # verified in step 2 (kept as the record of what was checked), but this
    # analysis is nitrogen-only from step 1's reopen decision.
    experiments = experiments[experiments["nutrient"] == "nitrogen"]

    resolution = pd.read_csv(STEP2_DATA / "02_gene_locus_resolution.csv")
    mcra = resolution[(resolution["cog_number"] == COG) & (resolution["status"] == "resolved")]
    print(f"{COG} (McrA) resolved loci (all 15 strains, before restricting to nitrogen-experiment strains):")
    print(mcra[["organism_name", "locus_tag"]].to_string(index=False))

    # GOTCHA (see gaps_and_friction.md): differential_expression_by_gene
    # treats experiment_ids=[] as "no filter", not "match nothing". Must
    # only query organisms that actually have a nitrogen experiment in
    # scope, never pass an empty exp_ids list.
    in_scope_organisms = set(experiments["organism_name"])
    mcra_in_scope = mcra[mcra["organism_name"].isin(in_scope_organisms)]
    print(f"\nRestricted to the {len(in_scope_organisms)} strain(s) with an in-scope "
          f"nitrogen experiment: {sorted(in_scope_organisms)}")

    rows = []
    with GraphConnection() as conn:
        for organism_full, group in mcra_in_scope.groupby("organism_name"):
            organism_short = organism_full.replace("Prochlorococcus ", "")
            exp_ids = experiments[experiments["organism_name"] == organism_full]["experiment_id"].tolist()
            assert exp_ids, f"{organism_full} has no in-scope experiment_ids -- would silently disable the filter"
            loci = group["locus_tag"].unique().tolist()
            result = differential_expression_by_gene(
                organism=organism_short, locus_tags=loci, experiment_ids=exp_ids,
                limit=None, conn=conn,
            )
            de = pd.DataFrame(result["results"])
            rows.append(de)
            print(f"\n{organism_short}: {len(de)} DE rows for loci {loci}")

    de_df = pd.concat(rows, ignore_index=True)
    out_path = DATA_DIR / "01_driving_example_mcra.csv"
    de_df.to_csv(out_path, index=False)

    # independent hand tally
    n_up = (de_df["expression_status"] == "significant_up").sum()
    n_down = (de_df["expression_status"] == "significant_down").sum()
    n_not_sig = (de_df["expression_status"] == "not_significant").sum()
    print(f"\nHand tally: {len(de_df)} rows, {n_up} up, {n_down} down, {n_not_sig} not significant "
          f"({100 * (n_up + n_down) / len(de_df):.1f}% significant)")

    result = hit_rate(de_df)
    print(f"hit_rate() result: {result}")

    assert result["n_up"] == n_up
    assert result["n_down"] == n_down
    assert result["n_not_significant"] == n_not_sig
    print("\nMATCH -- hit_rate() output matches the independent hand tally exactly.")


if __name__ == "__main__":
    main()
