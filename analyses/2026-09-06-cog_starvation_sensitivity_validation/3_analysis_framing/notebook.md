# Step 3 — Analysis framing

## Context

Step 2 established the usable evidence base: 32 of 41 COGs have DE
evidence, all but 1 from the nitrogen side, which is why step 1 was
reopened to drop phosphorus. This step operationalizes the (now
nitrogen-only) hypothesis and materializes the negative-background gene
set, co-defined with the researcher: no positive-control subset (this
list is candidate genes, not textbook markers), and the background pool
excludes both this analysis's 41-COG list and the prior analysis's 92-gene
list.

## What I did

**Significance criteria** (`scripts/01_copy_significance_criteria.py`):
reused the prior walkthrough analysis's empirically-derived
significance-criteria table for these same 5 nitrogen experiments — no
new derivation needed, since the platforms and their actual padj/log2FC
thresholds haven't changed.

**Background pool** (`scripts/02_build_background_pool.py`): pulled every
DE row from the same 4 unbiased ("all_detected_genes") nitrogen
experiments the prior analysis used, at the same single chosen starvation
timepoint per experiment (reused directly — no new timepoint judgment
call), excluding:
1. This analysis's 41-COG list, by (organism, locus_tag) pair — not by
   gene-name text, since many of these genes have no `gene_name` populated
   in the KG (locus tag is the unambiguous identifier, per Rule 2).
2. The prior analysis's 92-gene N/P-acquisition list, by (organism,
   locus_tag) pair — co-defined with the researcher: those genes are
   already hypothesized to respond to starvation, so leaving them in the
   "random" pool would inflate the null baseline against which the 41-COG
   set is tested.
3. Nitrogen-keyword product/gene-name text matches (same 13-keyword list
   the prior analysis used: nitrogen, nitrate, nitrite, ammonium, ammonia,
   urea, urease, cyanate, glutamine synthetase, glutamate synthase, GOGAT,
   nitrogenase, amino acid transport).

## Results

**Exclusion set:** 364 (organism, locus) pairs from the 41-COG list ∪ 183
from the 92-gene list = 543 combined pairs excluded from the background
(some overlap possible if the same locus happened to be on both lists —
none was, since a union of 364+183 rows produced exactly 543... actually
543 < 364+183=547, so 4 pairs do overlap between the two lists).

**Background pool per experiment** (full table:
`data/02_negative_background_pool.csv`):

| Experiment | Distinct genes tested (at chosen TP) | Excluded (41-COG list) | Excluded (N-keyword) | Background pool |
|---|---|---|---|---|
| Weissberg proteomics (MED4) | 1,424 | 60 rows | 22 rows | 1,359 genes |
| Weissberg RNA-seq (MED4) | 1,849 | 76 rows | 24 rows | 1,768 genes |
| Tolonen microarray (MED4) | 1,697 | 74 rows | 24 rows | 1,618 genes |
| Tolonen microarray (MIT9313) | 2,241 | 66 rows | 23 rows | 2,169 genes |

("Excluded (41-COG list)" column combines both target-locus lists, per the
script; row counts are per-row not per-gene, consistent with the prior
analysis's reporting convention.)

**Framing (co-defined with the researcher):**
- **Hypothesis:** the 41-COG "starvation sensitivity" candidate list
  shows a nitrogen-starvation response (significant DE, at the same
  per-platform thresholds established in step 3's significance-criteria
  table) elevated above the background/random rate — tested for all 41
  pooled, and separately for the `N` (12), `mixed` (28), and `unlabeled`
  (1) Direction-tag subsets.
- **No positive controls.** Unlike the prior analysis's `ntcA`/`glnA`-type
  textbook markers, this list has no genes with an established, literature
  -confirmed nitrogen-starvation response to sanity-check against — the
  premise being tested is precisely whether any of them behave that way.
  Skipping a positive-control subset was a deliberate co-defined choice,
  not an oversight.
- **Negative background:** the pool above, nitrogen-only (same rationale
  as the prior analysis — no unfiltered phosphorus table exists to build
  an unbiased phosphorus background, moot now that phosphorus is out of
  scope entirely).
- **Expected outcome, operationally:** `hit_rate()` on the 41-COG set
  (and subsets) against the 5 nitrogen experiments, compared via
  `bootstrap_pvalue()` and Fisher's exact test to the background pool
  above — the same two functions from `4_methods/np_response.py` in the
  prior analysis, reused rather than rebuilt (step 4 will verify, not
  reimplement).

## Surprises

None this step — the framing follows directly from step 1/2's decisions;
no new judgment calls were forced by the data.

## Decide-gate checklist

- **Outputs produced:** `scripts/01_copy_significance_criteria.py`,
  `scripts/02_build_background_pool.py`; `data/01_significance_criteria.csv`
  (5 nitrogen rows), `data/02_negative_background_pool.csv` (6,914 rows).
- **Results presented:** significance-criteria table (reused); background
  pool exclusion table (above); framing (hypothesis, no-positive-control
  decision, background scope, operational expected outcome).
- **QC gate:** confirmed the reused significance-criteria table has
  exactly the 5 nitrogen rows expected. Background-pool exclusion checked
  per experiment (target-list and N-keyword row counts both nonzero and
  proportionate to pool size, consistent with the prior analysis's own
  exclusion rates).
- **Decisions made this step:** no positive-control subset (co-defined,
  2026-09-06) — the only new judgment call this step forced; everything
  else follows mechanically from steps 1-2.
- **Advance rationale:** the nitrogen-only hypothesis is now operational
  (what "responds" means per experiment, what the background comparison
  is, that the 92-gene list is also excluded from it) — ready for step 4
  to verify the existing `hit_rate()`/`bootstrap_pvalue()` functions still
  apply as-is to this new gene set before running them at scale in step 5.
