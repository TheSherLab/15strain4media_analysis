"""Collapse the (gene x experiment x timepoint) rows in
06_nitrogen_volcano_data.csv down to one row per (gene, platform) --
so a gene measured 19 times (like ntcA, across 2 strains x multiple
timepoints) counts once per platform, not 19 times, when comparing
genes.

For each (gene_name_raw, omics_type): n measurements, n significant,
fraction significant, and the median log2FC across ALL measurements
(not just significant ones -- avoids cherry-picking only the biggest
hits, which would bias the summary toward overstating each gene's
typical effect).

Input:  data/06_nitrogen_volcano_data.csv
Output: data/09_gene_level_summary.csv (one row per gene x platform)

Usage:
  uv run python analyses/2026-08-02-np_gene_expression_full_analysis/5_analyze/scripts/09_pull_gene_level_summary.py
"""
from pathlib import Path

import pandas as pd

STEP_DIR = Path(__file__).resolve().parents[1]
DATA_DIR = STEP_DIR / "data"


def main() -> None:
    df = pd.read_csv(DATA_DIR / "06_nitrogen_volcano_data.csv")

    grouped = df.groupby(["gene_name_raw", "n_or_p", "omics_type"])
    summary = grouped.agg(
        n_measurements=("log2fc", "size"),
        n_significant=("padj", lambda s: (s < 0.05).sum()),
        median_log2fc=("log2fc", "median"),
    ).reset_index()
    summary["fraction_significant"] = summary["n_significant"] / summary["n_measurements"]

    out_path = DATA_DIR / "09_gene_level_summary.csv"
    summary.to_csv(out_path, index=False)
    print(f"Wrote {len(summary)} (gene x platform) rows to {out_path}")
    print(f"\nDistinct genes: {df['gene_name_raw'].nunique()} "
          f"(was {len(df)} raw measurement rows -> {len(summary)} gene x platform rows)")
    print("\nRows per platform:")
    print(summary.groupby("omics_type").size().to_string())


if __name__ == "__main__":
    main()
