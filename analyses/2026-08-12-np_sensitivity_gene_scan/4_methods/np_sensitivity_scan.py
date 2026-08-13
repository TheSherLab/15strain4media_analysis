"""
Ad-hoc methods module for the genome-wide N/P sensitivity gene scan (step 4).

Implements, per step 3's framing (../3_analysis_framing/notebook.md):
  - presence/absence association test (Fisher's exact)
  - copy-number association test (Mann-Whitney U, genome-size-normalized)
  - Benjamini-Hochberg FDR correction
  - direction and significance-classification logic

Every test function is grouping-agnostic: callers pass the two strain-name
lists to compare (the sensitivity split or the ecotype split), so step 5
calls each function twice per ortholog group -- once per grouping -- through
the same code path. FDR correction is a separate function applied afterward,
once per (test type x grouping) combination, across the full batch of
p-values from that combination -- not something a single-group test call can
do on its own (co-defined with the researcher: 4 separate corrections, not
pooled, since presence/absence and copy-number cover different, only
partly-overlapping sets of ortholog groups).
"""
from __future__ import annotations

from statistics import median

from scipy.stats import fisher_exact, mannwhitneyu
from statsmodels.stats.multitest import multipletests


def normalize_copy_number(copy_count: int, strain_gene_count: int) -> float:
    """Copies per 1,000 genes, so a strain's larger genome doesn't mechanically
    inflate its copy count relative to smaller-genome strains."""
    return copy_count / strain_gene_count * 1000


def presence_absence_test(
    copy_counts: dict[str, int],
    group_a: list[str],
    group_b: list[str],
    label_a: str,
    label_b: str,
) -> dict:
    """Fisher's exact test on carries-vs-not, group_a vs group_b.

    copy_counts: strain name -> copy count in this ortholog group (0 if absent).
    Every strain in group_a/group_b must have an entry in copy_counts.
    """
    carries_a = sum(1 for s in group_a if copy_counts[s] > 0)
    not_a = len(group_a) - carries_a
    carries_b = sum(1 for s in group_b if copy_counts[s] > 0)
    not_b = len(group_b) - carries_b

    _, p_value = fisher_exact([[carries_a, not_a], [carries_b, not_b]], alternative="two-sided")

    rate_a = carries_a / len(group_a)
    rate_b = carries_b / len(group_b)
    if rate_a > rate_b:
        direction = f"enriched in {label_a}"
    elif rate_b > rate_a:
        direction = f"enriched in {label_b}"
    else:
        direction = "equal"

    return {
        "test": "presence_absence",
        "p_value": p_value,
        "direction": direction,
        "n_carry_a": carries_a,
        "n_total_a": len(group_a),
        "rate_a": rate_a,
        "n_carry_b": carries_b,
        "n_total_b": len(group_b),
        "rate_b": rate_b,
        "label_a": label_a,
        "label_b": label_b,
    }


def is_copy_number_variable(copy_counts: dict[str, int], strains: list[str]) -> bool:
    """True if raw copy count differs across at least two of the given strains.
    Strains that don't carry the group count as 0 (co-defined with the researcher:
    'variable' is judged across all strains, absence included, not just among carriers)."""
    return len({copy_counts[s] for s in strains}) > 1


def copy_number_test(
    copy_counts: dict[str, int],
    strain_gene_counts: dict[str, int],
    group_a: list[str],
    group_b: list[str],
    label_a: str,
    label_b: str,
) -> dict | None:
    """Mann-Whitney U test on genome-size-normalized copy number, group_a vs group_b.

    Returns None if copy count doesn't vary across all of group_a + group_b
    (per step 3: this test runs "where variable"; see is_copy_number_variable).
    """
    all_strains = group_a + group_b
    if not is_copy_number_variable(copy_counts, all_strains):
        return None

    norm_a = [normalize_copy_number(copy_counts[s], strain_gene_counts[s]) for s in group_a]
    norm_b = [normalize_copy_number(copy_counts[s], strain_gene_counts[s]) for s in group_b]

    _, p_value = mannwhitneyu(norm_a, norm_b, alternative="two-sided")

    median_a = median(norm_a)
    median_b = median(norm_b)
    if median_a > median_b:
        direction = f"higher in {label_a}"
    elif median_b > median_a:
        direction = f"higher in {label_b}"
    else:
        direction = "equal"

    return {
        "test": "copy_number",
        "p_value": p_value,
        "direction": direction,
        "median_normalized_a": median_a,
        "median_normalized_b": median_b,
        "label_a": label_a,
        "label_b": label_b,
    }


def benjamini_hochberg(p_values: list[float]) -> list[float]:
    """BH FDR-adjusted p-values, same order and length as the input list.
    Called once per (test type x grouping) combination across the full batch
    of ortholog groups tested by that combination -- not per-group."""
    _, adjusted, _, _ = multipletests(p_values, method="fdr_bh")
    return list(adjusted)


def classify_significance(sens_padj: float | None, eco_padj: float | None, alpha: float = 0.05) -> str | None:
    """Per step 3: 'candidate' means sensitivity-significant and ecotype-not-significant.
    Genes significant on both, or ecotype-only, stay labeled but aren't candidates --
    nothing is filtered out of the output table. Returns None if either adjusted
    p-value is unavailable (e.g. the copy-number test was skipped as not variable)."""
    if sens_padj is None or eco_padj is None:
        return None
    sens_sig = sens_padj < alpha
    eco_sig = eco_padj < alpha
    if sens_sig and not eco_sig:
        return "candidate"
    if sens_sig and eco_sig:
        return "significant_both"
    if eco_sig and not sens_sig:
        return "ecotype_only"
    return "not_significant"
