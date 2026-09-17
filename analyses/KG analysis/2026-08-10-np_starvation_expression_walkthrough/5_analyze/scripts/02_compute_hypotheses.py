"""
Step 5 -- compute hypotheses 1 (matched), 2 (cross), and 3 (nitrogen-only
noise/bootstrap check) from the extracted gene x experiment matrix.

Inputs: data/01_target_gene_experiment_matrix.csv,
        ../../3_analysis_framing/data/02_positive_controls.csv,
        ../../3_analysis_framing/data/02_negative_background_pool.csv
Outputs: data/02_hypothesis_results.csv (summary table),
         data/02_bootstrap_null_distribution.csv (H3 raw null draws)

Usage: uv run analyses/2026-08-10-np_starvation_expression_walkthrough/5_analyze/scripts/02_compute_hypotheses.py
"""

import sys
import time
from pathlib import Path

import pandas as pd
from scipy.stats import fisher_exact

sys.path.insert(0, str((Path(__file__).resolve().parent.parent.parent / "4_methods")))
from np_response import hit_rate, bootstrap_pvalue  # noqa: E402

BASE = Path(__file__).resolve().parent.parent
DATA_DIR = BASE / "data"
FRAMING_DATA = BASE.parent / "3_analysis_framing" / "data"

TESTED_STATUSES = ["significant_up", "significant_down", "not_significant"]


def main() -> None:
    matrix = pd.read_csv(DATA_DIR / "01_target_gene_experiment_matrix.csv")
    tested = matrix[matrix["status"].isin(TESTED_STATUSES)].copy()
    tested["expression_status"] = tested["status"]  # np_response.hit_rate expects this column name

    results = []

    def add_result(label, subset):
        r = hit_rate(subset)
        r["label"] = label
        r["n_distinct_genes"] = subset["gene_name"].nunique()
        results.append(r)
        return r

    # H1: matched response, nitrogen and phosphorus separately
    matched_n = tested[(tested["match_type"] == "matched") & (tested["n_or_p"] == "N")]
    matched_p = tested[(tested["match_type"] == "matched") & (tested["n_or_p"] == "P")]
    add_result("H1: N-genes in N-experiments (matched)", matched_n)
    add_result("H1: P-genes in P-experiments (matched)", matched_p)

    # H2: cross response
    cross_n_in_p = tested[(tested["match_type"] == "cross") & (tested["n_or_p"] == "N")]
    cross_p_in_n = tested[(tested["match_type"] == "cross") & (tested["n_or_p"] == "P")]
    add_result("H2: N-genes in P-experiments (cross)", cross_n_in_p)
    add_result("H2: P-genes in N-experiments (cross)", cross_p_in_n)

    # Positive controls, for reference
    pos_controls = pd.read_csv(FRAMING_DATA / "02_positive_controls.csv")
    pos_n_genes = set(pos_controls[pos_controls["nutrient"] == "N"]["gene_name"])
    pos_p_genes = set(pos_controls[pos_controls["nutrient"] == "P"]["gene_name"])
    add_result("Positive controls: N genes in N-experiments", matched_n[matched_n["gene_name"].isin(pos_n_genes)])
    add_result("Positive controls: P genes in P-experiments", matched_p[matched_p["gene_name"].isin(pos_p_genes)])

    results_df = pd.DataFrame(results)[
        ["label", "n_distinct_genes", "n_tests", "n_up", "n_down", "n_not_significant",
         "pct_significant", "pct_up", "pct_down"]
    ]
    print("Hit-rate results:")
    print(results_df.to_string(index=False))

    # H3: nitrogen-only noise/bootstrap check + Fisher's exact, matched-N vs. background
    print("\n--- H3: nitrogen noise/bootstrap check ---")
    bg_pool = pd.read_csv(FRAMING_DATA / "02_negative_background_pool.csv")
    n_experiments = matched_n["experiment_id"].unique().tolist()
    genes_per_experiment = {
        exp_id: matched_n[matched_n["experiment_id"] == exp_id]["locus_tag"].nunique()
        for exp_id in n_experiments
    }
    print(f"Target genes per nitrogen experiment (for size-matched bootstrap draws): {genes_per_experiment}")

    background_pools = {exp_id: bg_pool[bg_pool["experiment_id"] == exp_id] for exp_id in n_experiments}
    observed = hit_rate(matched_n)

    t0 = time.time()
    boot = bootstrap_pvalue(
        observed_hit_rate=observed["pct_significant"],
        genes_per_experiment=genes_per_experiment,
        background_pools=background_pools,
        n_iterations=10_000,
        seed=42,
    )
    elapsed = time.time() - t0
    print(f"Bootstrap (10,000 iterations) took {elapsed:.2f}s wall-clock")
    print(f"Observed hit rate: {observed['pct_significant']:.1f}%")
    print(f"Null distribution: mean={boot['null_mean']:.1f}%, std={boot['null_std']:.1f}%")
    print(f"Empirical p-value: {boot['p_value']:.4f}")

    pd.DataFrame({"null_hit_rate": boot["null_distribution"]}).to_csv(
        DATA_DIR / "02_bootstrap_null_distribution.csv", index=False
    )

    # Fisher's exact: target-N significant/not vs. background significant/not (pooled across experiments)
    bg_sig = bg_pool["expression_status"].isin(["significant_up", "significant_down"]).sum()
    bg_not_sig = (bg_pool["expression_status"] == "not_significant").sum()
    target_sig = observed["n_significant"]
    target_not_sig = observed["n_not_significant"]
    table = [[target_sig, target_not_sig], [bg_sig, bg_not_sig]]
    odds_ratio, fisher_p = fisher_exact(table)
    print(f"\nFisher's exact test (target N matched vs. background): "
          f"table={table}, odds_ratio={odds_ratio:.2f}, p={fisher_p:.2e}")

    results_df.to_csv(DATA_DIR / "02_hypothesis_results.csv", index=False)
    print(f"\nWrote {DATA_DIR / '02_hypothesis_results.csv'}")

    with open(DATA_DIR / "02_h3_summary.txt", "w") as f:
        f.write(f"Observed hit rate (N matched): {observed['pct_significant']:.2f}%\n")
        f.write(f"Bootstrap null mean: {boot['null_mean']:.2f}%, std: {boot['null_std']:.2f}%\n")
        f.write(f"Bootstrap p-value (10,000 iterations, {elapsed:.2f}s): {boot['p_value']:.4f}\n")
        f.write(f"Fisher's exact: table={table}, odds_ratio={odds_ratio:.3f}, p={fisher_p:.3e}\n")
        f.write(f"genes_per_experiment: {genes_per_experiment}\n")
    print(f"Wrote {DATA_DIR / '02_h3_summary.txt'}")


if __name__ == "__main__":
    main()
