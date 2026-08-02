"""Fisher's exact tests on the two valid enrichment-style comparisons in
this analysis -- both against the nitrogen background gene set, the
only unbiased background this KG build has:

1. H1_N vs background_N: do nitrogen target genes (matched) hit more
   than background, in nitrogen experiments? (also serves as the H3
   noise check for nitrogen)
2. H2_Pgenes_Nexp vs background_N: do phosphorus-annotated genes
   (cross-nutrient) hit more than background, in nitrogen experiments?

No phosphorus-background test is possible: step 5 deliberately did not
compute a phosphorus background set, since it would be drawn from the
same significance-pre-filtered tables as every phosphorus target gene
(see 5_analyze/notebook.md), so H1_P and H2_Ngenes_Pexp have no valid
background to test against -- their rates are reported descriptively
only, not statistically tested here.

2x2 tables (from data/01_hit_rate_summary.csv):

                        significant   not significant
    target N genes (H1)      116            287        = 403
    P genes in N exp (H2)     62            333        = 395
    background N genes       236           1360        = 1596

Uses scipy.stats.fisher_exact directly (not the package's fisher_ora,
which is built for many-term GO/pathway ORA with per-cluster
backgrounds -- overkill for two one-off 2x2 comparisons sharing one
background).

Input:  ../5_analyze/data/01_hit_rate_summary.csv
Output: data/01_fisher_test_nitrogen.csv

Usage:
  uv run python analyses/2026-08-02-np_gene_expression_full_analysis/6_evaluate/scripts/01_fisher_test.py
"""
from pathlib import Path

import pandas as pd
from scipy.stats import fisher_exact

STEP_DIR = Path(__file__).resolve().parents[1]
DATA_DIR = STEP_DIR / "data"
ANALYZE_DATA = STEP_DIR.parent / "5_analyze" / "data"


def run_test(name: str, target: pd.Series, background: pd.Series) -> dict:
    target_sig = int(target["n_up"] + target["n_down"])
    target_ns = int(target["n_not_significant"])
    bg_sig = int(background["n_up"] + background["n_down"])
    bg_ns = int(background["n_not_significant"])

    table = [[target_sig, target_ns], [bg_sig, bg_ns]]
    odds_ratio, p_value = fisher_exact(table, alternative="two-sided")

    target_rate = target_sig / (target_sig + target_ns)
    bg_rate = bg_sig / (bg_sig + bg_ns)
    fold_enrichment = target_rate / bg_rate

    print(f"--- {name} ---")
    print(f"                     significant   not_significant   total")
    print(f"  target             {target_sig:>11}   {target_ns:>15}   {target_sig + target_ns}")
    print(f"  background N genes {bg_sig:>11}   {bg_ns:>15}   {bg_sig + bg_ns}")
    print(f"  target rate: {target_rate:.4f}  background rate: {bg_rate:.4f}  "
          f"fold enrichment: {fold_enrichment:.2f}x  odds ratio: {odds_ratio:.3f}  "
          f"p (two-sided): {p_value:.3e}")
    print()

    return {
        "comparison": name,
        "target_significant": target_sig,
        "target_not_significant": target_ns,
        "background_significant": bg_sig,
        "background_not_significant": bg_ns,
        "target_rate": target_rate,
        "background_rate": bg_rate,
        "fold_enrichment": fold_enrichment,
        "odds_ratio": odds_ratio,
        "p_value": p_value,
    }


def main() -> None:
    summary = pd.read_csv(ANALYZE_DATA / "01_hit_rate_summary.csv").set_index("group")
    bgn = summary.loc["background_N"]

    rows = [
        run_test("H1_N vs background_N (matched N genes, also serves as H3)", summary.loc["H1_N"], bgn),
        run_test("H2_Pgenes_Nexp vs background_N (cross: P genes under N starvation)", summary.loc["H2_Pgenes_Nexp"], bgn),
    ]

    out = pd.DataFrame(rows)
    out_path = DATA_DIR / "01_fisher_test_nitrogen.csv"
    out.to_csv(out_path, index=False)
    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
