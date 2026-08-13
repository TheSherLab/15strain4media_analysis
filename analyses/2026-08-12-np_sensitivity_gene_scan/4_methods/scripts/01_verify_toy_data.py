"""
Step 4 -- verify the methods module's test functions on hand-built toy data
with known, planted patterns before trusting them on the real ~5,000+
ortholog-group table (step 5).

Toy strains are synthetic (not the real 15 study strains) so each case can
be a clean, hand-reasoned textbook example, decoupled from the real
sensitivity/ecotype confound (that confound is real biology, demonstrated
instead on real data in 02_driving_example.py). Every case below asserts
against a hand-derived expectation; the script fails loudly (AssertionError)
if the module's behavior doesn't match.

Usage:
  uv run python analyses/2026-08-12-np_sensitivity_gene_scan/4_methods/scripts/01_verify_toy_data.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from np_sensitivity_scan import (  # noqa: E402
    benjamini_hochberg,
    classify_significance,
    copy_number_test,
    is_copy_number_variable,
    normalize_copy_number,
    presence_absence_test,
)

GROUP_A9 = [f"a{i}" for i in range(1, 10)]  # 9 toy strains, mirrors the real 9-vs-6 sensitivity split size
GROUP_B6 = [f"b{i}" for i in range(1, 7)]  # 6 toy strains
GENE_COUNT = {s: 2000 for s in GROUP_A9 + GROUP_B6}  # flat gene count so normalization doesn't distort these cases


def case_1_clear_presence_signal() -> None:
    """8 of 9 group-A strains carry it, only 1 of 6 group-B strains does -- should read
    as a strong, significant enrichment in A."""
    copy_counts = {s: 1 for s in GROUP_A9[:8]}
    copy_counts[GROUP_A9[8]] = 0
    copy_counts.update({s: 0 for s in GROUP_B6[:5]})
    copy_counts[GROUP_B6[5]] = 1

    result = presence_absence_test(copy_counts, GROUP_A9, GROUP_B6, "A", "B")
    print("case 1 (clear presence signal):", result)
    assert result["p_value"] < 0.05, "expected a significant p-value for an 8/9 vs 1/6 split"
    assert result["direction"] == "enriched in A"
    assert result["n_carry_a"] == 8 and result["n_carry_b"] == 1


def case_2_no_presence_signal() -> None:
    """Roughly equal carry-rates in both groups -- should read as not significant."""
    copy_counts = {s: 1 for s in GROUP_A9[:4]}
    copy_counts.update({s: 0 for s in GROUP_A9[4:]})
    copy_counts.update({s: 1 for s in GROUP_B6[:3]})
    copy_counts.update({s: 0 for s in GROUP_B6[3:]})

    result = presence_absence_test(copy_counts, GROUP_A9, GROUP_B6, "A", "B")
    print("case 2 (no presence signal):", result)
    assert result["p_value"] > 0.05, "expected a non-significant p-value for a near-even split"


def case_3_universal_presence_no_error() -> None:
    """Present in every strain -- degenerate 2x2 (one column all-zero); must not error,
    and must correctly report no distinguishing signal (p == 1.0)."""
    copy_counts = {s: 1 for s in GROUP_A9 + GROUP_B6}
    result = presence_absence_test(copy_counts, GROUP_A9, GROUP_B6, "A", "B")
    print("case 3 (universal presence, edge case):", result)
    assert result["p_value"] == 1.0
    assert result["direction"] == "equal"


def case_4_clear_copy_number_signal() -> None:
    """Group A strains all carry 2 copies, group B strains all carry 1 -- should read as
    a strong, significant, higher-in-A copy-number signal."""
    copy_counts = {s: 2 for s in GROUP_A9}
    copy_counts.update({s: 1 for s in GROUP_B6})

    assert is_copy_number_variable(copy_counts, GROUP_A9 + GROUP_B6)
    result = copy_number_test(copy_counts, GENE_COUNT, GROUP_A9, GROUP_B6, "A", "B")
    print("case 4 (clear copy-number signal):", result)
    assert result is not None
    assert result["p_value"] < 0.05, "expected a significant p-value for a clean 2-copy vs 1-copy split"
    assert result["direction"] == "higher in A"
    assert result["median_normalized_a"] > result["median_normalized_b"]


def case_5_copy_number_not_variable_is_skipped() -> None:
    """Every strain carries exactly 1 copy -- no variation to test; the function must
    return None rather than run a meaningless test."""
    copy_counts = {s: 1 for s in GROUP_A9 + GROUP_B6}
    assert not is_copy_number_variable(copy_counts, GROUP_A9 + GROUP_B6)
    result = copy_number_test(copy_counts, GENE_COUNT, GROUP_A9, GROUP_B6, "A", "B")
    print("case 5 (copy number not variable):", result)
    assert result is None


def case_6_variability_counts_absence_as_zero() -> None:
    """Every carrier has exactly 1 copy, but 3 of the 15 strains don't carry it at all
    (copy count 0). Co-defined rule: absence counts as 0 when judging variability, so
    this MUST be treated as variable even though every carrier's own count is identical."""
    copy_counts = {s: 1 for s in GROUP_A9 + GROUP_B6}
    copy_counts[GROUP_A9[0]] = 0
    copy_counts[GROUP_A9[1]] = 0
    copy_counts[GROUP_B6[0]] = 0
    assert is_copy_number_variable(copy_counts, GROUP_A9 + GROUP_B6)
    result = copy_number_test(copy_counts, GENE_COUNT, GROUP_A9, GROUP_B6, "A", "B")
    print("case 6 (variability via absence, not just among carriers):", result)
    assert result is not None


def case_7_fdr_correction_wiring() -> None:
    """Sanity-check the BH wiring: adjusted p-values are >= raw p-values (BH's core
    property), same length/order as input, and small raw p-values stay small."""
    raw = [0.001, 0.01, 0.02, 0.20, 0.50, 0.80, 0.90]
    adjusted = benjamini_hochberg(raw)
    print("case 7 (FDR correction):", list(zip(raw, adjusted)))
    assert len(adjusted) == len(raw)
    assert all(a >= r for a, r in zip(adjusted, raw)), "BH-adjusted p-values must never be smaller than raw"
    assert adjusted[0] < 0.05, "the smallest raw p-value should survive correction here"


def case_8_classification() -> None:
    """Sensitivity-significant + ecotype-not -> candidate; both significant -> flagged but
    not a candidate; ecotype-only -> flagged but not a candidate; neither -> not significant;
    a missing (None) p-value (skipped copy-number test) -> no classification, not an error."""
    assert classify_significance(0.01, 0.50) == "candidate"
    assert classify_significance(0.01, 0.01) == "significant_both"
    assert classify_significance(0.50, 0.01) == "ecotype_only"
    assert classify_significance(0.50, 0.50) == "not_significant"
    assert classify_significance(None, 0.50) is None
    assert classify_significance(0.01, None) is None
    print("case 8 (classification): all 6 sub-cases passed")


def main() -> None:
    case_1_clear_presence_signal()
    case_2_no_presence_signal()
    case_3_universal_presence_no_error()
    case_4_clear_copy_number_signal()
    case_5_copy_number_not_variable_is_skipped()
    case_6_variability_counts_absence_as_zero()
    case_7_fdr_correction_wiring()
    case_8_classification()
    print("\nAll toy-data checks passed.")


if __name__ == "__main__":
    main()
