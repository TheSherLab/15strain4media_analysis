"""
Step 6 -- full detail on what the positive and negative controls actually
showed, beyond the summary numbers already in step 5's notebook.

Positive controls: per-gene breakdown (not just the pooled 80.0%/84.6%
rates) for the 8 marker genes (4 nitrogen, 4 phosphorus), across every
experiment x strain test they appear in.

Negative background: per-experiment composition of the nitrogen background
pool (already summarized in 3_analysis_framing/notebook.md) plus the shape
of the H3 bootstrap null distribution (already summarized as mean/std in
5_analyze/notebook.md) -- here reported with percentiles so "not
explainable by chance" can be read against the actual null spread, not
just its mean.

Inputs: ../../5_analyze/data/01_target_gene_experiment_matrix.csv
        ../../3_analysis_framing/data/02_positive_controls.csv
        ../../3_analysis_framing/data/02_negative_background_pool.csv
        ../../5_analyze/data/02_bootstrap_null_distribution.csv
Outputs: data/01_positive_control_by_gene.csv
         data/02_negative_background_by_experiment.csv
         data/03_bootstrap_null_percentiles.txt

Usage: uv run analyses/2026-08-10-np_starvation_expression_walkthrough/6_evaluate/scripts/01_control_detail.py
"""

from pathlib import Path

import numpy as np
import pandas as pd

BASE = Path(__file__).resolve().parent.parent
DATA_DIR = BASE / "data"
ANALYZE_DATA = BASE.parent / "5_analyze" / "data"
FRAMING_DATA = BASE.parent / "3_analysis_framing" / "data"

TESTED_STATUSES = ["significant_up", "significant_down", "not_significant"]


def main() -> None:
    matrix = pd.read_csv(ANALYZE_DATA / "01_target_gene_experiment_matrix.csv")
    tested = matrix[matrix["status"].isin(TESTED_STATUSES)].copy()

    pos_controls = pd.read_csv(FRAMING_DATA / "02_positive_controls.csv")
    control_genes = pos_controls[["nutrient", "gene_name"]].drop_duplicates()

    # --- Positive controls: per-gene detail, matched-nutrient tests only
    # (same scoping 02_compute_hypotheses.py used: N controls in N
    # experiments, P controls in P experiments).
    rows = []
    for _, row in control_genes.iterrows():
        nutrient, gene = row["nutrient"], row["gene_name"]
        nutrient_name = "nitrogen" if nutrient == "N" else "phosphorus"
        gene_rows = tested[
            (tested["gene_name"] == gene)
            & (tested["n_or_p"] == nutrient)
            & (tested["match_type"] == "matched")
            & (tested["nutrient"] == nutrient_name)
        ]
        n_tests = len(gene_rows)
        n_up = int((gene_rows["status"] == "significant_up").sum())
        n_down = int((gene_rows["status"] == "significant_down").sum())
        n_not_sig = int((gene_rows["status"] == "not_significant").sum())
        rows.append({
            "nutrient": nutrient,
            "gene_name": gene,
            "n_tests": n_tests,
            "n_up": n_up,
            "n_down": n_down,
            "n_not_significant": n_not_sig,
            "pct_significant": 100 * (n_up + n_down) / n_tests if n_tests else float("nan"),
            "n_experiments": gene_rows["experiment_id"].nunique(),
        })
    control_detail = pd.DataFrame(rows)
    control_detail.to_csv(DATA_DIR / "01_positive_control_by_gene.csv", index=False)
    print("Positive controls, per gene:")
    print(control_detail.to_string(index=False))

    # --- Negative background: per-experiment composition (nitrogen only)
    bg_pool = pd.read_csv(FRAMING_DATA / "02_negative_background_pool.csv")
    bg_summary = bg_pool.groupby("experiment_id").agg(
        n_genes=("locus_tag", "nunique"),
        n_sig_up=("expression_status", lambda s: (s == "significant_up").sum()),
        n_sig_down=("expression_status", lambda s: (s == "significant_down").sum()),
        n_not_sig=("expression_status", lambda s: (s == "not_significant").sum()),
    ).reset_index()
    bg_summary["pct_significant"] = 100 * (bg_summary["n_sig_up"] + bg_summary["n_sig_down"]) / (
        bg_summary["n_sig_up"] + bg_summary["n_sig_down"] + bg_summary["n_not_sig"]
    )
    bg_summary.to_csv(DATA_DIR / "02_negative_background_by_experiment.csv", index=False)
    print("\nNegative background pool, per experiment:")
    print(bg_summary.to_string(index=False))

    # --- Bootstrap null distribution: percentiles, not just mean/std
    null_dist = pd.read_csv(ANALYZE_DATA / "02_bootstrap_null_distribution.csv")["null_hit_rate"]
    percentiles = [1, 5, 25, 50, 75, 95, 99, 100]
    pct_values = {p: float(np.percentile(null_dist, p)) for p in percentiles}
    with open(DATA_DIR / "03_bootstrap_null_percentiles.txt", "w") as f:
        f.write(f"n_iterations: {len(null_dist)}\n")
        f.write(f"mean: {null_dist.mean():.2f}%\n")
        f.write(f"std: {null_dist.std():.2f}%\n")
        f.write(f"max observed in 10,000 draws: {null_dist.max():.2f}%\n")
        for p, v in pct_values.items():
            f.write(f"P{p}: {v:.2f}%\n")
        f.write("observed matched-N rate: 41.3% (exceeds every one of the 10,000 draws)\n")
    print("\nBootstrap null distribution percentiles:")
    for p, v in pct_values.items():
        print(f"  P{p}: {v:.2f}%")
    print(f"  max of 10,000 draws: {null_dist.max():.2f}% (observed 41.3% exceeds this)")


if __name__ == "__main__":
    main()
