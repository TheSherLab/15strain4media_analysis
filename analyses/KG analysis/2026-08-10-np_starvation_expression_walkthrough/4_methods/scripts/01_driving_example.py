"""
Step 4 -- driving example: ntcA and glnA in MED4's 4 nitrogen experiments.

Purpose: verify hit_rate() against real KG data for two well-known
nitrogen-starvation marker genes, printing every individual row so the
result can be eyeballed and hand-tallied before np_response.py is trusted
on the full 60-gene analysis in step 5.

Usage: uv run analyses/2026-08-10-np_starvation_expression_walkthrough/4_methods/scripts/01_driving_example.py
"""

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from np_response import hit_rate  # noqa: E402

from multiomics_explorer import differential_expression_by_gene, GraphConnection

MED4_NITROGEN_EXPERIMENTS = [
    "10.1101/2025.11.24.690089_growth_state_pro99lown_nutrient_starvation_med4_proteomics_axenic",
    "10.1101/2025.11.24.690089_growth_state_pro99lown_nutrient_starvation_med4_rnaseq_axenic",
    "10.1038/ismej.2017.88_nitrogen_stress_ndepleted_pro99_medium_med4_rnaseq",
    "10.1038/msb4100087_nitrogen_nitrogen_deprivation_med4_med4_microarray",
]
GENES = {"ntcA": "PMM0246", "glnA": "PMM0920"}


def main() -> None:
    with GraphConnection() as conn:
        result = differential_expression_by_gene(
            organism="MED4",
            locus_tags=list(GENES.values()),
            experiment_ids=MED4_NITROGEN_EXPERIMENTS,
            limit=None,
            verbose=True,
            conn=conn,
        )
    de = pd.DataFrame(result["results"])
    print(f"Pulled {len(de)} DE rows for ntcA (PMM0246) + glnA (PMM0920) "
          f"across MED4's 4 nitrogen experiments\n")

    print("Every row (hand-tally these against the summary below):")
    cols = ["gene_name", "locus_tag", "experiment_name", "timepoint", "log2fc", "padj", "expression_status"]
    print(de[cols].sort_values(["gene_name", "experiment_name", "timepoint"]).to_string(index=False))

    print("\nPer-gene breakdown:")
    for name, tag in GENES.items():
        gene_rows = de[de["locus_tag"] == tag]
        r = hit_rate(gene_rows)
        print(f"  {name} ({tag}): {r['n_tests']} tests, {r['n_up']} up, {r['n_down']} down, "
              f"{r['n_not_significant']} not significant -> {r['pct_significant']:.1f}% significant")

    overall = hit_rate(de)
    print(f"\nPooled (both genes, hit_rate() output): {overall}")

    manual_n_tests = len(de)
    manual_n_up = int((de["expression_status"] == "significant_up").sum())
    manual_n_down = int((de["expression_status"] == "significant_down").sum())
    assert overall["n_tests"] == manual_n_tests
    assert overall["n_up"] == manual_n_up
    assert overall["n_down"] == manual_n_down
    print("\n[PASS] hit_rate() output matches an independent manual tally of the same rows.")


if __name__ == "__main__":
    main()
