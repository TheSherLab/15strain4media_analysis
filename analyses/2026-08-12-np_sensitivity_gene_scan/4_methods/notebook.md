# Step 4 — Methods

## Context

Step 3 locked the test design: Fisher's exact for presence/absence, Mann-Whitney U for
genome-size-normalized copy number, both run against the sensitivity split (9 N-sensitive vs. 6
mixed) and the ecotype split (10 HL vs. 5 LL), with Benjamini-Hochberg FDR correction and a
`candidate` tag for sensitivity-significant-but-not-ecotype-significant genes. This step turns
that design into an ad-hoc Python module, verified on hand-built toy data before running on the
real ~5,000+ ortholog groups.

## What I did

**Co-defined two judgment calls with the researcher before writing code:**
1. **FDR pooling.** Presence/absence and copy-number are different tests over different (only
   partly overlapping) sets of ortholog groups — copy number is only tested where it varies at
   all. Agreed: 4 separate BH corrections, one per (test type × grouping) combination
   (presence/absence×sensitivity, presence/absence×ecotype, copy-number×sensitivity,
   copy-number×ecotype), not 2 pooled-by-grouping corrections.
2. **"Variable" copy number.** Agreed: build each group's copy-count vector across all 15
   strains, counting strains that don't carry the group as 0, and run the copy-number test if
   that vector has more than one distinct value (i.e. absence counts toward variability, not
   just variation among carriers).

**Built the module** (`np_sensitivity_scan.py`): `normalize_copy_number`, `presence_absence_test`
(scipy `fisher_exact`), `is_copy_number_variable`, `copy_number_test` (scipy `mannwhitneyu`),
`benjamini_hochberg` (statsmodels `multipletests`, method `fdr_bh`), `classify_significance`.
Every test function is grouping-agnostic — callers pass in whichever two strain lists and labels
to compare, so step 5 will call the same functions once for sensitivity and once for ecotype.

**Verified on toy data** (`scripts/01_verify_toy_data.py`): 8 hand-built cases using synthetic
strain names (decoupled from the real sensitivity/ecotype confound, so each case has a clean,
hand-reasoned expected outcome) — a clear presence signal, a null presence case, a degenerate
universal-presence edge case (checks no divide-by-zero/error), a clear copy-number signal, a
not-variable copy-number case (confirms the test is skipped, not run meaninglessly), a
variable-via-absence case (confirms the co-defined variability rule), FDR-correction wiring
(adjusted ≥ raw, correct length), and all 6 classification outcomes including the `None`/skipped
case. All 8 passed (see script output in Results).

**Ran the driving example** (`scripts/02_driving_example.py`): both tests, both groupings, on one
real ortholog group from step 2's data (`2_kg_selection/data/03_gene_group_assignment.csv`) and
step 3's strain-group table (`3_analysis_framing/data/01_strain_group_membership.csv`). Group
selection rule, fixed before looking at any result: the first ortholog group (sorted by
`final_group_id`) with real presence/absence variation (carried by more than 1 but fewer than all
15 strains) *and* real copy-number variation (at least one strain with >1 copy) — chosen so the
demo exercises both test functions, not a trivial universal/single-copy case. FDR correction and
classification are **not** run in this example — both need the full batch of p-values from a
given (test type × grouping) combination, which only exists once step 5 runs the loop across all
groups; this example reports raw p-values and direction only.

## Results

**Toy verification** (`scripts/01_verify_toy_data.py` output):

```
case 1 (clear presence signal): p=0.0110, direction=enriched in A, 8/9 vs 1/6
case 2 (no presence signal): p=1.0000, 4/9 vs 3/6
case 3 (universal presence, edge case): p=1.0000, direction=equal, 9/9 vs 6/6
case 4 (clear copy-number signal): p=0.00024, direction=higher in A, median 1.0 vs 0.5
case 5 (copy number not variable): None (test skipped)
case 6 (variability via absence, not just among carriers): test ran, p=0.865, direction=equal
case 7 (FDR correction): raw [0.001,0.01,0.02,0.20,0.50,0.80,0.90] -> adjusted
  [0.007,0.035,0.047,0.35,0.70,0.90,0.90] — all adjusted >= raw
case 8 (classification): all 6 sub-cases passed
All toy-data checks passed.
```

**Driving example** — ortholog group `cyanorak:CK_00000001` (gene name `rpoD8`, RNA polymerase
sigma factor type II), carried by 14 of 15 strains (absent only in MIT1314), with 2 copies in
MIT1327 and MIT9313, 1 copy elsewhere (full table: `data/02_driving_example_readout.csv`):

| Test | Grouping | p-value | Direction |
|---|---|---|---|
| Presence/absence | sensitivity (9 vs 6) | 1.000 | enriched in mixed (8/9 vs 6/6) |
| Presence/absence | ecotype (10 vs 5) | 1.000 | enriched in LL (9/10 vs 5/5) |
| Copy number | sensitivity | 0.328 | higher in mixed |
| Copy number | ecotype | 0.859 | higher in HL |

This particular gene shows no signal on either grouping or either test (expected — it's an
arbitrarily, mechanically selected group, not a result of interest) and confirms the pipeline
runs end to end on real data without error, including a real near-universal-presence case
(14/15) landing correctly at p=1.0 rather than erroring.

## Surprises

- None. The toy cases and the driving example both behaved as hand-reasoned/expected.

## Decisions

**2026-08-13 — 4 separate FDR corrections, one per (test type × grouping) combination.**
Co-defined with the researcher (see Context/What I did). Rationale: presence/absence and
copy-number test different, only partly-overlapping sets of ortholog groups, so pooling their
p-values into one correction would mix incomparable test families.

**2026-08-13 — Copy-number "variable" means the per-strain copy-count vector (absent = 0) across
all strains being compared has more than one distinct value.** Co-defined with the researcher.
Rationale: simpler, literal reading of "variable copy count across strains"; treats absence as
part of the copy-number picture rather than a separate axis.

## Decide-gate checklist

- **Outputs produced:** `np_sensitivity_scan.py` (methods module);
  `scripts/01_verify_toy_data.py`, `scripts/02_driving_example.py`;
  `data/02_driving_example_readout.csv`.
- **Results presented:** toy-verification output (8/8 passed), driving-example readout table
  (all above).
- **QC gate:** toy script asserts all 8 hand-derived expectations and exits without error;
  driving-example script asserts the sensitivity split resolves to 9/6 and the ecotype split to
  10/5 (matching step 3's locked grouping) before running any test; degenerate 2×2 (universal
  presence) case confirmed not to error and to report p=1.0 as expected, both in toy data and in
  the real driving example (14/15 case).
- **Decisions made this step:** FDR-pooling choice (4 separate corrections); copy-number
  variability rule (absence counts as 0) — both above, dated 2026-08-13.
- **Advance rationale:** the module implements step 3's design exactly, is verified against
  hand-reasoned toy expectations covering the normal cases and the two edge cases (universal
  presence, non-variable copy number), and runs cleanly end-to-end on one real ortholog group.
  Ready for step 5 to loop it across all tested ortholog groups and apply FDR correction across
  the full batch.
