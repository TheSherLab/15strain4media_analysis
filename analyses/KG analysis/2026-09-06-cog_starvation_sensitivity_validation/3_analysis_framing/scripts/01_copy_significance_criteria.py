"""
Step 3 -- reuse (not re-derive) the significance-criteria table the prior
walkthrough analysis already built empirically for these same 5 nitrogen
experiments (what padj/log2FC threshold each platform actually uses,
whether the source table reports every gene or a filtered subset).

Inputs: ../../../2026-08-10-np_starvation_expression_walkthrough/3_analysis_framing/data/01_significance_criteria.csv
Outputs: data/01_significance_criteria.csv -- nitrogen rows only

Usage: uv run analyses/2026-09-06-cog_starvation_sensitivity_validation/3_analysis_framing/scripts/01_copy_significance_criteria.py
"""

from pathlib import Path

import pandas as pd

BASE = Path(__file__).resolve().parent.parent
DATA_DIR = BASE / "data"
PRIOR_CSV = (
    BASE.parent.parent
    / "2026-08-10-np_starvation_expression_walkthrough"
    / "3_analysis_framing" / "data" / "01_significance_criteria.csv"
)


def main() -> None:
    df = pd.read_csv(PRIOR_CSV)
    n_df = df[df["nutrient"] == "nitrogen"].reset_index(drop=True)
    assert len(n_df) == 5, f"expected 5 nitrogen rows, got {len(n_df)}"

    out_path = DATA_DIR / "01_significance_criteria.csv"
    n_df.to_csv(out_path, index=False)
    print(f"Wrote {len(n_df)} nitrogen-experiment significance-criteria rows to {out_path}")
    print(n_df[["experiment_id", "statistical_test", "table_scope"]].to_string(index=False))


if __name__ == "__main__":
    main()
