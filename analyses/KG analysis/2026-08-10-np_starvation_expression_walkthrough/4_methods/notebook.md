# Step 4 — Methods

## Context

Step 3 operationalized the three hypotheses and materialized the control
gene sets. This step builds the two functions needed to compute them —
`hit_rate()` (hypotheses 1 and 2) and `bootstrap_pvalue()` (hypothesis 3,
nitrogen only) — and verifies both before they run on real analysis data
in step 5.

## What I did

Wrote `np_response.py` (module, not a script — imported by step 5) with
two functions:

- `hit_rate(de_df)`: given a DataFrame of DE rows (one row per gene x
  experiment x timepoint), returns the count and percentage significant,
  split by up/down direction.
- `bootstrap_pvalue(...)`: draws random, size-matched gene samples from
  the nitrogen background pool (per experiment, matching how many target
  genes had evidence in that experiment) many times, computes the hit
  rate for each draw, and returns an empirical p-value — the fraction of
  random draws that reach or exceed the observed target hit rate.

**Toy verification** (`scripts/00_toy_verification.py`, before touching
real data): `hit_rate()` checked against a 10-row synthetic table
hand-tallied to 3 up / 2 down / 5 not-significant (50% significant, 30%
up, 20% down — matched exactly). `bootstrap_pvalue()` checked against a
100-gene synthetic background where exactly 20 genes (20%) are flagged
significant: sampling 5 genes 2,000 times, the null distribution's mean
landed at 20.5% (expected ~20%); an observed rate of 100% (all 5 sampled
genes significant) got p=0.0005 (correctly rare against a 20% background);
an observed rate of exactly 20% (matching the population) got p=0.68
(correctly unremarkable). All checks passed.

**Driving example** (`scripts/01_driving_example.py`, real KG data): pulled
every DE row for `ntcA` and `glnA` across MED4's 4 nitrogen experiments
(26 rows total), printed each one, and hand-tallied significant/direction
counts independently of `hit_rate()`'s own computation — the two matched
exactly (see Results). Both genes are known nitrogen-starvation markers,
so this also served as a sanity check on the extraction itself: both come
back overwhelmingly upregulated, as expected.

**Runtime check** (the researcher asked how long 10,000 bootstrap
iterations would take): the bootstrap operates entirely on the
already-extracted background-pool CSV from step 3 — no KG queries inside
the loop. 2,000 iterations in the toy verification ran in well under a
second of actual compute time (`user 0m0.015s` for the whole toy-check
script, most of the 5.8s wall-clock being Python/uv startup overhead, not
computation). 10,000 iterations on the real background pool is expected
to be similarly fast; step 5 will report the actual time.

## Results

**Driving example — every row, ntcA and glnA, MED4, 4 nitrogen
experiments:**

| Gene | Experiment | Timepoint | log2FC | padj | Status |
|---|---|---|---|---|---|
| glnA | Read RNASEQ | 3h | 1.50 | 1.0 | not_significant |
| glnA | Read RNASEQ | 12h | 2.85 | 0.0 | significant_up |
| glnA | Read RNASEQ | 24h | 2.45 | 0.0 | significant_up |
| glnA | Tolonen MICROARRAY | 0h | 0.51 | 1.0 | not_significant |
| glnA | Tolonen MICROARRAY | 3h | 3.76 | 0.01 | significant_up |
| glnA | Tolonen MICROARRAY | 6h | 5.06 | 0.01 | significant_up |
| glnA | Tolonen MICROARRAY | 12h | 4.14 | 0.01 | significant_up |
| glnA | Tolonen MICROARRAY | 24h | 3.92 | 0.01 | significant_up |
| glnA | Tolonen MICROARRAY | 48h | 3.80 | 0.01 | significant_up |
| glnA | Weissberg PROTEOMICS | day 14 | 2.12 | 8.9e-05 | significant_up |
| glnA | Weissberg PROTEOMICS | day 31 | 1.38 | 4.5e-04 | significant_up |
| glnA | Weissberg PROTEOMICS | day 89 | 2.53 | 3.5e-06 | significant_up |
| glnA | Weissberg RNASEQ | day 14 | 3.36 | 7.4e-06 | significant_up |
| ntcA | Read RNASEQ | 3h | 0.57 | 1.0 | not_significant |
| ntcA | Read RNASEQ | 12h | 1.71 | 0.26 | not_significant |
| ntcA | Read RNASEQ | 24h | 2.17 | 0.0 | significant_up |
| ntcA | Tolonen MICROARRAY | 0h | -2.17 | 1.0 | not_significant |
| ntcA | Tolonen MICROARRAY | 3h | 3.04 | 1.0 | not_significant |
| ntcA | Tolonen MICROARRAY | 6h | 4.45 | 0.01 | significant_up |
| ntcA | Tolonen MICROARRAY | 12h | 3.02 | 0.01 | significant_up |
| ntcA | Tolonen MICROARRAY | 24h | 1.68 | 0.01 | significant_up |
| ntcA | Tolonen MICROARRAY | 48h | 1.23 | 1.0 | not_significant |
| ntcA | Weissberg PROTEOMICS | day 14 | 1.87 | 2.2e-04 | significant_up |
| ntcA | Weissberg PROTEOMICS | day 31 | 3.59 | 2.5e-07 | significant_up |
| ntcA | Weissberg PROTEOMICS | day 89 | 2.68 | 2.5e-06 | significant_up |
| ntcA | Weissberg RNASEQ | day 14 | 3.35 | 1.0e-05 | significant_up |

**Summary:** ntcA 13 tests, 8 up, 0 down, 5 not-significant (61.5%
significant); glnA 13 tests, 11 up, 0 down, 2 not-significant (84.6%
significant); pooled 26 tests, 19 significant (all up), 73.1%. Both
biologically expected (both are canonical nitrogen-starvation-induced
genes) and exactly matched an independent manual tally of the same rows.

## Surprises

None this step — both toy and real-data verification matched expectations
on the first attempt.

## Decide-gate checklist

- **Outputs produced:** `np_response.py` (methods module);
  `scripts/00_toy_verification.py`, `scripts/01_driving_example.py`.
- **Results presented:** toy-check pass/fail table; ntcA/glnA full row
  table and summary (above).
- **QC gate:** `hit_rate()` matches independent manual tally on both toy
  data and real KG data (ntcA/glnA); `bootstrap_pvalue()` matches expected
  behavior on synthetic data (extreme observed rate -> small p-value;
  population-matching observed rate -> non-extreme p-value).
- **Decisions made this step:** none (toy verification and driving example
  ran as planned; no forks required).
- **Advance rationale:** both methods functions are verified against
  hand-computable ground truth before touching the full 60-gene analysis;
  ready for step 5 to run them at scale and produce the actual result
  tables and figures.
