"""
Step 4 -- toy-data verification of np_response.py before real data.

Hand-computed synthetic examples for hit_rate() and bootstrap_pvalue(),
per the repo's methodology (verify reusable utilities on toy data first).

Usage: uv run analyses/2026-08-10-np_starvation_expression_walkthrough/4_methods/scripts/00_toy_verification.py
"""

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from np_response import hit_rate, bootstrap_pvalue  # noqa: E402


def check(label: str, got, expected, tol=1e-9):
    ok = abs(got - expected) < tol if isinstance(expected, (int, float)) else got == expected
    status = "PASS" if ok else "FAIL"
    print(f"  [{status}] {label}: got {got}, expected {expected}")
    assert ok, f"{label} mismatch"


def test_hit_rate():
    print("hit_rate() -- 10-row synthetic table, hand-tallied: 3 up, 2 down, 5 not_significant")
    toy = pd.DataFrame({
        "expression_status": (
            ["significant_up"] * 3 + ["significant_down"] * 2 + ["not_significant"] * 5
        )
    })
    result = hit_rate(toy)
    check("n_tests", result["n_tests"], 10)
    check("n_up", result["n_up"], 3)
    check("n_down", result["n_down"], 2)
    check("n_significant", result["n_significant"], 5)
    check("pct_significant", result["pct_significant"], 50.0)
    check("pct_up", result["pct_up"], 30.0)
    check("pct_down", result["pct_down"], 20.0)

    print("hit_rate() -- empty table")
    empty_result = hit_rate(pd.DataFrame({"expression_status": []}))
    check("n_tests (empty)", empty_result["n_tests"], 0)


def test_bootstrap():
    print("\nbootstrap_pvalue() -- background pool where exactly 20% of genes are "
          "significant, sampling 5 genes 2000 times")
    # 100-gene background pool, first 20 flagged significant_up, rest not_significant --
    # every gene contributes exactly 1 row (1 timepoint), so pct_significant per
    # draw is just (# of the 20 "hot" genes drawn) / (# drawn) * 100.
    genes = [f"G{i:03d}" for i in range(100)]
    statuses = ["significant_up"] * 20 + ["not_significant"] * 80
    pool = pd.DataFrame({
        "locus_tag": genes, "experiment_id": "toy_exp", "expression_status": statuses,
    })

    # Observed: a "target set" of 5 genes, all significant (100% hit rate) --
    # should be a very unlikely draw from a 20%-hot pool, so p-value should be small.
    result_high = bootstrap_pvalue(
        observed_hit_rate=100.0,
        genes_per_experiment={"toy_exp": 5},
        background_pools={"toy_exp": pool},
        n_iterations=2000,
        seed=1,
    )
    print(f"  observed=100%, background hot rate=20% -> null_mean={result_high['null_mean']:.1f}"
          f" (expected ~20), p_value={result_high['p_value']:.4f} (expected small, <0.05)")
    check("null_mean near population rate (100 draws, tol 3pp)", round(result_high["null_mean"]), 20, tol=3)
    assert result_high["p_value"] < 0.05, "expected a small p-value for an extreme observed rate"
    print("  [PASS] p_value is small, as expected for an extreme observed hit rate")

    # Observed: exactly the population rate (20%) -- p-value should be roughly in the
    # middle of the distribution, not extreme.
    result_typical = bootstrap_pvalue(
        observed_hit_rate=20.0,
        genes_per_experiment={"toy_exp": 5},
        background_pools={"toy_exp": pool},
        n_iterations=2000,
        seed=2,
    )
    print(f"  observed=20% (== population rate) -> p_value={result_typical['p_value']:.4f} "
          f"(expected roughly 0.3-0.7, not extreme)")
    assert 0.2 < result_typical["p_value"] < 0.8, "expected a non-extreme p-value at the population rate"
    print("  [PASS] p_value is not extreme, as expected when observed == population rate")


if __name__ == "__main__":
    test_hit_rate()
    test_bootstrap()
    print("\nAll toy-data checks passed.")
