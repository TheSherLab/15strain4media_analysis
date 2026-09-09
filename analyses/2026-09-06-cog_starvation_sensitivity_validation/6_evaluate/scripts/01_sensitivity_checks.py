"""
Step 6 -- robustness checks on the step-5 nitrogen result, triggered by
two features of the "use every gene per COG" decision:

  1. Big-COG dominance: 10 COGs hold ~70% of the gene rows (COG0477 alone
     is 25-34 MFS permeases per strain). Does dropping COGs with >8 genes
     in a strain change the pooled nitrogen hit rate / its background test?
  2. Per-COG (not per-gene) unit: recompute the nitrogen hit rate with
     each COG collapsed to one value per experiment (share of its genes
     significant, thresholded at >=50%), so a 30-gene COG counts once like
     a 1-gene COG -- closer to how the copy-number test treats a COG.

Also: per-COG significant-gene fractions (nitrogen), for the caveats.

Inputs: ../../5_analyze/data/01_target_gene_experiment_matrix.csv,
        ../../3_analysis_framing/data/02_negative_background_pool.csv
Outputs: data/01_sensitivity_summary.txt, data/02_per_cog_nitrogen_fractions.csv

Usage: uv run analyses/2026-09-06-cog_starvation_sensitivity_validation/6_evaluate/scripts/01_sensitivity_checks.py
"""

import sys
from pathlib import Path

import pandas as pd
from scipy.stats import fisher_exact

BASE = Path(__file__).resolve().parent.parent
DATA_DIR = BASE / "data"
ANALYZE_DATA = BASE.parent / "5_analyze" / "data"
FRAMING_DATA = BASE.parent / "3_analysis_framing" / "data"
sys.path.insert(0, str(BASE.parent.parent / "2026-08-10-np_starvation_expression_walkthrough" / "4_methods"))
from np_response import hit_rate, bootstrap_pvalue  # noqa: E402

TESTED = ["significant_up", "significant_down", "not_significant"]
BIG_COG_MAX_GENES = 8


def main() -> None:
    m = pd.read_csv(ANALYZE_DATA / "01_target_gene_experiment_matrix.csv").rename(
        columns={"status": "expression_status"})
    n = m[m["nutrient"] == "nitrogen"]
    n_tested = n[n["expression_status"].isin(TESTED)]

    bg = pd.read_csv(FRAMING_DATA / "02_negative_background_pool.csv")
    bg_sig = bg["expression_status"].isin(["significant_up", "significant_down"]).sum()
    bg_not = (bg["expression_status"] == "not_significant").sum()
    bg_by_exp = {e: bg[bg["experiment_id"] == e] for e in n_tested["experiment_id"].unique()}

    def test(df, label):
        obs = hit_rate(df)
        gpe = df.groupby("experiment_id")["locus_tag"].nunique().to_dict()
        boot = bootstrap_pvalue(obs["pct_significant"], gpe, bg_by_exp, n_iterations=10_000, seed=42)
        odds, fp = fisher_exact([[obs["n_significant"], obs["n_not_significant"]], [bg_sig, bg_not]])
        return (f"{label}: {obs['n_significant']}/{obs['n_tests']} = {obs['pct_significant']:.1f}%  "
                f"vs background 19.1%  |  bootstrap p={boot['p_value']:.3f}, "
                f"Fisher OR={odds:.2f} p={fp:.3f}")

    lines = ["Nitrogen sensitivity checks (background overall rate 19.1%)", ""]
    lines.append(test(n_tested, "all genes, all 41 COGs (step 5 headline)"))

    # 1. drop big COGs
    genes_per_cog_strain = (n_tested.groupby(["cog_number", "organism_name"])["locus_tag"]
                            .nunique().groupby("cog_number").max())
    big = sorted(genes_per_cog_strain[genes_per_cog_strain > BIG_COG_MAX_GENES].index)
    lines.append("")
    lines.append(f"COGs with >{BIG_COG_MAX_GENES} genes in a strain ({len(big)}): {big}")
    lines.append(test(n_tested[~n_tested["cog_number"].isin(big)],
                      f"excluding those {len(big)} big COGs"))

    # 2. per-COG unit: one value per (COG, experiment) = >=50% of genes significant
    lines.append("")
    per_ce = (n_tested.assign(sig=n_tested["expression_status"].isin(["significant_up", "significant_down"]))
              .groupby(["cog_number", "experiment_id"])["sig"].mean())
    per_ce_call = (per_ce >= 0.5)
    n_cog_tests = len(per_ce_call)
    n_cog_sig = int(per_ce_call.sum())
    lines.append(f"per-COG unit (COG 'responds' if >=50% of its genes significant in that experiment): "
                 f"{n_cog_sig}/{n_cog_tests} = {100*n_cog_sig/n_cog_tests:.1f}%  "
                 f"(background genes at the same >=50%-of-1 threshold would be ~19%)")

    (DATA_DIR / "01_sensitivity_summary.txt").write_text("\n".join(lines), encoding="utf-8")
    print("\n".join(lines))

    # per-COG nitrogen fractions for the caveats table
    frac = (n_tested.assign(sig=n_tested["expression_status"].isin(["significant_up", "significant_down"]),
                            up=n_tested["expression_status"].eq("significant_up"),
                            down=n_tested["expression_status"].eq("significant_down"))
            .groupby(["cog_number", "direction"])
            .agg(n_tests=("sig", "size"), n_sig=("sig", "sum"), n_up=("up", "sum"), n_down=("down", "sum")))
    frac["pct_sig"] = (100 * frac["n_sig"] / frac["n_tests"]).round(1)
    frac = frac.reset_index().sort_values("pct_sig", ascending=False)
    frac.to_csv(DATA_DIR / "02_per_cog_nitrogen_fractions.csv", index=False)
    print("\nTop nitrogen-responsive COGs:")
    print(frac.head(8).to_string(index=False))


if __name__ == "__main__":
    main()
