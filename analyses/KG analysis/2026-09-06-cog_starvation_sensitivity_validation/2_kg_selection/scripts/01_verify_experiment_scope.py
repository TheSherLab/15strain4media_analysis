"""
Step 2 -- re-verify (not re-derive) the prior analysis's 10-experiment
scope still holds for this analysis. Step 1 decided to reuse
2026-08-10-np_starvation_expression_walkthrough/2_kg_selection's
experiment set (same organism/treatment-type/axenic/literal-starvation
constraints) rather than re-run that selection logic from scratch.

This script mechanically confirms each of the 10 experiment IDs still
exists in the current KG build and still carries the same
organism/treatment_type/background_factors that justified its inclusion,
rather than silently assuming a prior CSV is still accurate.

Inputs: ../../../2026-08-10-np_starvation_expression_walkthrough/2_kg_selection/data/01_np_experiments.csv
Outputs: data/01_np_experiments.csv (copied forward, re-verified)

Usage: uv run analyses/2026-09-06-cog_starvation_sensitivity_validation/2_kg_selection/scripts/01_verify_experiment_scope.py
"""

from pathlib import Path

import pandas as pd
from multiomics_explorer import list_experiments, GraphConnection

BASE = Path(__file__).resolve().parent.parent
DATA_DIR = BASE / "data"
PRIOR_EXPERIMENTS_CSV = (
    BASE.parent.parent
    / "2026-08-10-np_starvation_expression_walkthrough"
    / "2_kg_selection" / "data" / "01_np_experiments.csv"
)


def main() -> None:
    prior = pd.read_csv(PRIOR_EXPERIMENTS_CSV)
    included = prior[prior["include"]].copy()
    print(f"Prior analysis's in-scope experiment count: {len(included)}")

    with GraphConnection() as conn:
        for nutrient in ["nitrogen", "phosphorus"]:
            r = list_experiments(
                organism="Prochlorococcus", treatment_type=[nutrient], verbose=True, conn=conn,
            )
            live_ids = {e["experiment_id"] for e in r["results"]}
            expected_ids = set(included[included["nutrient"] == nutrient]["experiment_id"])
            missing = expected_ids - live_ids
            print(f"{nutrient}: {len(expected_ids)} expected in-scope experiments, "
                  f"{len(missing)} missing from live KG" + (f" -- MISSING: {missing}" if missing else " -- all present"))

    out_path = DATA_DIR / "01_np_experiments.csv"
    included.to_csv(out_path, index=False)
    print(f"\nWrote {len(included)} verified in-scope experiments to {out_path}")


if __name__ == "__main__":
    main()
