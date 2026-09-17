"""
Step 5 (redo 2026-09-09) -- hit rates for the 41-COG gene set.

  - nitrogen: pooled + N (13) / mixed (28) Direction subsets, each tested
    against the background pool via bootstrap + Fisher's exact
  - phosphorus: pooled + subsets, DESCRIPTIVE ONLY (no valid phosphorus
    background pool exists -- same reason as the walkthrough analysis;
    the phosphorus tables are pre-filtered / detection-limited, and the
    two Martiny microarrays only became usable via the reclassification)

`hit_rate()` / `bootstrap_pvalue()` reused unchanged from the walkthrough
analysis's 4_methods/np_response.py.

Inputs: data/01_target_gene_experiment_matrix.csv,
        ../../3_analysis_framing/data/02_negative_background_pool.csv
Outputs: data/02_hit_rate_results.csv, data/02_bootstrap_null_distribution.csv,
         data/02_h_summary.txt

Usage: uv run analyses/2026-09-06-cog_starvation_sensitivity_validation/5_analyze/scripts/02_compute_hit_rates.py
"""

import sys
from pathlib import Path

import pandas as pd
from scipy.stats import fisher_exact

BASE = Path(__file__).resolve().parent.parent
DATA_DIR = BASE / "data"
FRAMING_DATA = BASE.parent / "3_analysis_framing" / "data"
PRIOR_METHODS_DIR = BASE.parent.parent / "2026-08-10-np_starvation_expression_walkthrough" / "4_methods"
sys.path.insert(0, str(PRIOR_METHODS_DIR))
from np_response import hit_rate, bootstrap_pvalue  # noqa: E402

TESTED = ["significant_up", "significant_down", "not_significant"]


def main() -> None:
    m = pd.read_csv(DATA_DIR / "01_target_gene_experiment_matrix.csv").rename(
        columns={"status": "expression_status"})
    tested = m[m["expression_status"].isin(TESTED)]

    results = []
    for nutrient in ["nitrogen", "phosphorus"]:
        nut = tested[tested["nutrient"] == nutrient]
        for label, df in [
            (f"{nutrient} - all 41 pooled", nut),
            (f"{nutrient} - Direction N (13)", nut[nut["direction"] == "N"]),
            (f"{nutrient} - Direction mixed (28)", nut[nut["direction"] == "mixed"]),
        ]:
            r = hit_rate(df)
            r["group"] = label
            r["n_distinct_cogs"] = df["cog_number"].nunique()
            r["n_distinct_loci"] = df["locus_tag"].nunique()
            results.append(r)
            print(f"{label}: {r['n_significant']}/{r['n_tests']} = {r['pct_significant']:.1f}% "
                  f"({r['n_up']} up, {r['n_down']} down)")

    res_df = pd.DataFrame(results)[["group", "n_distinct_cogs", "n_distinct_loci", "n_tests",
                                    "n_significant", "n_up", "n_down", "n_not_significant",
                                    "pct_significant", "pct_up", "pct_down"]]
    res_df.to_csv(DATA_DIR / "02_hit_rate_results.csv", index=False)

    # --- nitrogen background comparison (pooled + subsets) ---
    bg = pd.read_csv(FRAMING_DATA / "02_negative_background_pool.csv")
    bg_sig = bg["expression_status"].isin(["significant_up", "significant_down"]).sum()
    bg_not = (bg["expression_status"] == "not_significant").sum()

    n_tested = tested[tested["nutrient"] == "nitrogen"]
    n_experiments = n_tested["experiment_id"].unique().tolist()
    bg_by_exp = {e: bg[bg["experiment_id"] == e] for e in n_experiments}

    summary = []
    null_dumps = {}
    for label, df in [
        ("nitrogen pooled (41)", n_tested),
        ("nitrogen Direction N (13)", n_tested[n_tested["direction"] == "N"]),
        ("nitrogen Direction mixed (28)", n_tested[n_tested["direction"] == "mixed"]),
    ]:
        obs = hit_rate(df)
        gpe = df.groupby("experiment_id")["locus_tag"].nunique().to_dict()
        boot = bootstrap_pvalue(obs["pct_significant"], gpe, bg_by_exp, n_iterations=10_000, seed=42)
        null_dumps[label] = boot["null_distribution"]
        odds, fp = fisher_exact([[obs["n_significant"], obs["n_not_significant"]], [bg_sig, bg_not]])
        summary.append(
            f"{label}: observed {obs['pct_significant']:.1f}% ({obs['n_significant']}/{obs['n_tests']}); "
            f"bootstrap null mean {boot['null_mean']:.1f}% (max {boot['null_distribution'].max():.1f}%), "
            f"p = {boot['p_value']:.4f}; Fisher OR {odds:.2f}, p = {fp:.3e}")
        print(summary[-1])

    pd.DataFrame(null_dumps).to_csv(DATA_DIR / "02_bootstrap_null_distribution.csv", index=False)

    p_pooled = hit_rate(tested[tested["nutrient"] == "phosphorus"])
    summary.append("")
    summary.append(f"phosphorus pooled (descriptive, no background test): "
                   f"{p_pooled['pct_significant']:.1f}% ({p_pooled['n_significant']}/{p_pooled['n_tests']}), "
                   f"{p_pooled['n_up']} up / {p_pooled['n_down']} down")
    summary.append(f"background pool overall: {100*bg_sig/(bg_sig+bg_not):.1f}% "
                   f"({bg_sig}/{bg_sig+bg_not})")
    (DATA_DIR / "02_h_summary.txt").write_text("\n".join(summary), encoding="utf-8")
    print("\n".join(summary[-3:]))


if __name__ == "__main__":
    main()
