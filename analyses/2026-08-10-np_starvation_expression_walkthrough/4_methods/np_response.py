"""
Methodology for the N/P starvation gene-expression response question.

Operates on already-extracted DE dataframes (no KG calls in this module --
extraction happens in the calling scripts via differential_expression_by_gene).
A DE dataframe is one row per (gene x experiment x timepoint), with at least
an `expression_status` column (values: "significant_up", "significant_down",
"not_significant").

Two functions:

- hit_rate(de_df): fraction of rows that are significant, split by direction.
  Answers hypotheses 1 and 2 (matched / cross response) directly -- run it
  once on the target gene set's DE rows for the matching-nutrient
  experiments (H1), and once for the non-matching-nutrient experiments (H2).

- bootstrap_pvalue(...): empirical p-value for whether the target gene set's
  hit rate exceeds what a random, size-matched gene set from the same
  experiments would achieve, built by resampling from a background pool.
  Nitrogen only (see 3_analysis_framing/notebook.md for why phosphorus has
  no valid background).

Toy-data verification: see 4_methods/notebook.md for a hand-computed check
of both functions against a synthetic 10-row example before either is run
on real data.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

SIGNIFICANT_STATUSES = {"significant_up", "significant_down"}


def hit_rate(de_df: pd.DataFrame) -> dict:
    """Fraction of (gene x experiment x timepoint) rows that are significant.

    Parameters
    ----------
    de_df : DataFrame with an `expression_status` column.

    Returns
    -------
    dict with n_tests, n_significant, n_up, n_down, n_not_significant,
    pct_significant, pct_up, pct_down (percentages of n_tests).
    """
    n_tests = len(de_df)
    if n_tests == 0:
        return {
            "n_tests": 0, "n_significant": 0, "n_up": 0, "n_down": 0,
            "n_not_significant": 0, "pct_significant": float("nan"),
            "pct_up": float("nan"), "pct_down": float("nan"),
        }
    n_up = int((de_df["expression_status"] == "significant_up").sum())
    n_down = int((de_df["expression_status"] == "significant_down").sum())
    n_sig = n_up + n_down
    n_not_sig = int((de_df["expression_status"] == "not_significant").sum())
    return {
        "n_tests": n_tests,
        "n_significant": n_sig,
        "n_up": n_up,
        "n_down": n_down,
        "n_not_significant": n_not_sig,
        "pct_significant": 100 * n_sig / n_tests,
        "pct_up": 100 * n_up / n_tests,
        "pct_down": 100 * n_down / n_tests,
    }


def bootstrap_pvalue(
    observed_hit_rate: float,
    genes_per_experiment: dict[str, int],
    background_pools: dict[str, pd.DataFrame],
    n_iterations: int = 10_000,
    seed: int = 42,
) -> dict:
    """Empirical p-value: how often does a random, size-matched gene sample
    from the background pool achieve a hit rate >= the observed one?

    Parameters
    ----------
    observed_hit_rate : the target gene set's actual pct_significant (0-100).
    genes_per_experiment : {experiment_id: n_target_genes_tested_in_this_experiment}
        -- how many target genes had evidence in each experiment, so each
        bootstrap draw samples the same number of genes per experiment as
        the real target set did (not a flat pooled count).
    background_pools : {experiment_id: DataFrame} -- each experiment's
        background gene pool (distinct locus_tag rows), from
        3_analysis_framing/data/02_negative_background_pool.csv grouped by
        experiment_id. Sampling is by distinct locus_tag, then all DE rows
        (all timepoints) for the sampled genes are pooled into that
        iteration's hit-rate calculation.
    n_iterations : number of bootstrap draws.
    seed : RNG seed for reproducibility.

    Returns
    -------
    dict with observed_hit_rate, null_mean, null_std, p_value (fraction of
    iterations with null hit rate >= observed), and the raw null_distribution
    array.
    """
    rng = np.random.default_rng(seed)
    null_hit_rates = np.empty(n_iterations)

    # Precompute each experiment's distinct genes once.
    pool_genes = {
        exp_id: pool["locus_tag"].unique() for exp_id, pool in background_pools.items()
    }

    for i in range(n_iterations):
        sampled_rows = []
        for exp_id, n_genes in genes_per_experiment.items():
            genes = pool_genes[exp_id]
            n_draw = min(n_genes, len(genes))
            chosen = rng.choice(genes, size=n_draw, replace=False)
            exp_pool = background_pools[exp_id]
            sampled_rows.append(exp_pool[exp_pool["locus_tag"].isin(chosen)])
        iteration_df = pd.concat(sampled_rows, ignore_index=True)
        null_hit_rates[i] = hit_rate(iteration_df)["pct_significant"]

    p_value = float(np.mean(null_hit_rates >= observed_hit_rate))
    return {
        "observed_hit_rate": observed_hit_rate,
        "null_mean": float(np.mean(null_hit_rates)),
        "null_std": float(np.std(null_hit_rates)),
        "n_iterations": n_iterations,
        "p_value": p_value,
        "null_distribution": null_hit_rates,
    }
