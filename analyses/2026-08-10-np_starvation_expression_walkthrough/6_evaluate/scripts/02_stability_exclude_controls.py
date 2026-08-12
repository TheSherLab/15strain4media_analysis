"""
Step 6 -- stability check: does the nitrogen matched-response finding (H1/H3)
depend on the 4 well-known positive-control genes (ntcA, glnA, amtB/amt1,
ureA), which individually score higher (80.0% significant) than the full
29-gene matched-N set (41.3%)?

Reruns hit_rate() and bootstrap_pvalue() on matched-N with the 4 positive
controls removed, using the same background pool and bootstrap machinery as
5_analyze/scripts/02_compute_hypotheses.py (same seed, same iteration count,
for direct comparability).

Inputs: ../../5_analyze/data/01_target_gene_experiment_matrix.csv
        ../../3_analysis_framing/data/02_positive_controls.csv
        ../../3_analysis_framing/data/02_negative_background_pool.csv
Outputs: data/04_stability_exclude_controls.txt

Usage: uv run analyses/2026-08-10-np_starvation_expression_walkthrough/6_evaluate/scripts/02_stability_exclude_controls.py
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
ANALYZE_DATA = BASE.parent / "5_analyze" / "data"
FRAMING_DATA = BASE.parent / "3_analysis_framing" / "data"

TESTED_STATUSES = ["significant_up", "significant_down", "not_significant"]


def main() -> None:
    matrix = pd.read_csv(ANALYZE_DATA / "01_target_gene_experiment_matrix.csv")
    tested = matrix[matrix["status"].isin(TESTED_STATUSES)].copy()
    tested["expression_status"] = tested["status"]

    pos_controls = pd.read_csv(FRAMING_DATA / "02_positive_controls.csv")
    pos_n_genes = set(pos_controls[pos_controls["nutrient"] == "N"]["gene_name"])

    matched_n = tested[(tested["match_type"] == "matched") & (tested["n_or_p"] == "N")]
    matched_n_no_controls = matched_n[~matched_n["gene_name"].isin(pos_n_genes)]

    full = hit_rate(matched_n)
    no_controls = hit_rate(matched_n_no_controls)

    print(f"Full matched-N set (29 genes): {full['n_tests']} tests, "
          f"{full['pct_significant']:.1f}% significant")
    print(f"Excluding 4 positive controls (25 genes): {no_controls['n_tests']} tests, "
          f"{no_controls['pct_significant']:.1f}% significant")

    # Bootstrap + Fisher on the reduced set, same machinery as step 5's H3.
    bg_pool = pd.read_csv(FRAMING_DATA / "02_negative_background_pool.csv")
    n_experiments = matched_n_no_controls["experiment_id"].unique().tolist()
    genes_per_experiment = {
        exp_id: matched_n_no_controls[matched_n_no_controls["experiment_id"] == exp_id]["locus_tag"].nunique()
        for exp_id in n_experiments
    }
    background_pools = {exp_id: bg_pool[bg_pool["experiment_id"] == exp_id] for exp_id in n_experiments}

    t0 = time.time()
    boot = bootstrap_pvalue(
        observed_hit_rate=no_controls["pct_significant"],
        genes_per_experiment=genes_per_experiment,
        background_pools=background_pools,
        n_iterations=10_000,
        seed=42,
    )
    elapsed = time.time() - t0
    print(f"\nBootstrap (10,000 iterations, no controls) took {elapsed:.2f}s")
    print(f"Null distribution: mean={boot['null_mean']:.1f}%, std={boot['null_std']:.1f}%")
    print(f"Empirical p-value: {boot['p_value']:.4f}")

    bg_sig = bg_pool["expression_status"].isin(["significant_up", "significant_down"]).sum()
    bg_not_sig = (bg_pool["expression_status"] == "not_significant").sum()
    table = [[no_controls["n_significant"], no_controls["n_not_significant"]], [bg_sig, bg_not_sig]]
    odds_ratio, fisher_p = fisher_exact(table)
    print(f"\nFisher's exact (no-controls target vs. background): "
          f"table={table}, odds_ratio={odds_ratio:.2f}, p={fisher_p:.2e}")

    with open(DATA_DIR / "04_stability_exclude_controls.txt", "w") as f:
        f.write(f"Full matched-N (29 genes, incl. 4 positive controls): "
                f"{full['n_tests']} tests, {full['n_up']} up, {full['n_down']} down, "
                f"{full['pct_significant']:.2f}% significant\n")
        f.write(f"Matched-N excluding 4 positive controls (25 genes): "
                f"{no_controls['n_tests']} tests, {no_controls['n_up']} up, {no_controls['n_down']} down, "
                f"{no_controls['pct_significant']:.2f}% significant\n")
        f.write(f"Bootstrap (10,000 iterations, no controls): null mean={boot['null_mean']:.2f}%, "
                f"std={boot['null_std']:.2f}%, p-value={boot['p_value']:.4f}\n")
        f.write(f"Fisher's exact (no controls): table={table}, odds_ratio={odds_ratio:.3f}, p={fisher_p:.3e}\n")
    print(f"\nWrote {DATA_DIR / '04_stability_exclude_controls.txt'}")


if __name__ == "__main__":
    main()
