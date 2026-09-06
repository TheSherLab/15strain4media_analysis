# Step 4 — Methods

## Context

Step 3 operationalized the nitrogen-only hypothesis and materialized the
background pool. This analysis needs no new method code — `hit_rate()`
and `bootstrap_pvalue()` from the prior walkthrough's
`4_methods/np_response.py` (already toy-verified there against
hand-computed synthetic data) apply unchanged to any DE dataframe shaped
the same way. This step's job is a driving example on real data from the
41-COG list, verified independently, before running at scale in step 5 —
not re-verifying the functions themselves, which are untouched.

## What I did

**Driving example** (`scripts/01_driving_example.py`): picked `COG1403`
(McrA, `N`-tagged), the COG with the most DE evidence (28 rows) of the 41.
Pulled every DE row across its resolved loci in MED4 (2 in-genome-duplicate
loci, `PMM0031`/`PMM1403`) and MIT9313 (`PMT0036`) — the only 2 strains
with an in-scope nitrogen experiment — and hand-tallied significant/
direction counts independently of `hit_rate()`'s own computation.

**Bug caught in the process**: the first version queried per organism
across all 15 strains (not per actual in-scope experiment), passing an
empty `experiment_ids` list for the 13 strains with no in-scope nitrogen
experiment. `differential_expression_by_gene(experiment_ids=[])` turned
out to behave as *no filter at all*, not "match nothing" — it returned an
unrelated 1-row hit for AS9601 from an out-of-scope experiment, inflating
the row count before the hand-tally cross-check caught the mismatch. Fixed
by restricting to the 2 organisms with a real in-scope experiment_ids list
before querying, with an assertion guarding against ever passing an empty
list again. Logged in `gaps_and_friction.md` as a standing caution for
future extraction scripts.

## Results

**Driving example — every row, McrA (COG1403), MED4 + MIT9313, the 5
nitrogen experiments** (full data: `data/01_driving_example_mcra.csv`):

23 total rows (17 from MED4's 2 duplicate loci across 3 experiments, 6
from MIT9313's 1 locus across its 1 experiment): 1 significant_up, 0
significant_down, 22 not_significant (4.3% significant).

`hit_rate()` output: `{'n_tests': 23, 'n_significant': 1, 'n_up': 1,
'n_down': 0, 'n_not_significant': 22, 'pct_significant': 4.35, 'pct_up':
4.35, 'pct_down': 0.0}` — **matches the independent hand tally exactly.**

Unlike the original analysis's `ntcA`/`glnA` driving example (which scored
overwhelmingly significant, as expected for known nitrogen markers),
McrA's low hit rate here (4.3%) is unsurprising and not concerning: McrA
is a restriction-modification/anti-viral defense gene with no established
nitrogen-starvation role — this driving example only verifies the
*function* computes correctly, not that this particular gene is expected
to respond.

## Surprises

- The `experiment_ids=[]` no-op-filter gotcha (above) — worth flagging
  because it's silent (no error, no warning) and would have gone unnoticed
  without the independent hand-tally check.

## Decide-gate checklist

- **Outputs produced:** `scripts/01_driving_example.py`;
  `data/01_driving_example_mcra.csv`.
- **Results presented:** McrA driving-example row counts and `hit_rate()`
  match (above).
- **QC gate:** `hit_rate()` output matches an independent manual tally of
  the real DE rows exactly. Caught and fixed the `experiment_ids=[]`
  no-filter gotcha before it could corrupt step 5's full run script.
- **Decisions made this step:** none (the bug fix was correctness work,
  not a judgment call).
- **Advance rationale:** the reused methods functions are verified against
  real data from this analysis's own gene list, and the empty-experiment-
  ids gotcha is caught and fixed here rather than silently corrupting
  step 5's full extraction — ready for step 5 to run at scale.
