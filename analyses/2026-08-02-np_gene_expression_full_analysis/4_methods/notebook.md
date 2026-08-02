# Step 4 — Methods

## Context

Step 3 defined three hypotheses and the target/positive-control/
negative-control gene sets, each needing the same underlying
computation: "for this set of genes, in this set of experiments, what
fraction of tests come back significantly differentially expressed, and
in which direction?" This step builds that computation once as a shared
method and verifies it against known biology before applying it in
step 5.

## Co-define

Proposed to the researcher: a single `hit_rate()` function, applied six
ways (H1 N/P, H2 N/P, positive control, negative control), with a
worked-example verification on two canonical nitrogen genes before
trusting it on the full gene sets. Agreed before building.

## What I did

- `4_methods/hit_rate.py` — the method module. Takes a list of
  (strain, locus_tag) pairs and a list of experiment IDs; calls
  `differential_expression_by_gene(significant_only=False, verbose=True)`
  batched per strain; tallies `expression_status` per gene into
  `n_up` / `n_down` / `n_not_significant` / `n_tests`, and computes
  `rate_significant`, `rate_up`, `rate_down` per gene and in aggregate.
  Genes with no row for a given experiment are excluded from that
  experiment's denominator (not measured there, not counted as
  "not significant").
- `scripts/01_worked_example.py` — traced `ntcA` (PMM0246) and `glnA`
  (PMM0920), MED4's two nitrogen positive-control genes, through
  `hit_rate()` against MED4's 4 nitrogen experiments, and independently
  hand-tallied the same raw rows to cross-check the function's counts.

## Results

**Worked example — `hit_rate()` vs. hand tally, MED4 nitrogen
experiments:**

| locus_tag | gene | n_tests | n_up | n_down | n_not_significant | rate_significant |
|---|---|---|---|---|---|---|
| PMM0246 | ntcA | 13 | 8 | 0 | 5 | 0.615 |
| PMM0920 | glnA | 13 | 11 | 0 | 2 | 0.846 |

The hand tally from the raw `differential_expression_by_gene` rows
(printed per-row in the script's output) matches `hit_rate()`'s counts
exactly for both genes (13 rows each, correct up/down/not-significant
split, `sum check: True`).

**Sanity check against known biology:** both genes are exclusively
`significant_up` when significant (0 `significant_down` rows for
either), consistent with `ntcA` (master nitrogen-starvation regulator)
and `glnA` (nitrogen-assimilation enzyme) being induced, not repressed,
under nitrogen starvation — the expected direction. Neither gene is
significant in 100% of tests (13/13); both have 2-5 not-significant
timepoints/experiments out of 13, which is expected — not every
timepoint of a time-course captures peak induction.

## Surprises

None this step — the method behaved as expected on the first pass, and
the worked-example cross-check matched exactly.

## Decisions

None this step.

## Decide-gate checklist

- **Outputs produced** — `hit_rate.py` (method module);
  `scripts/01_worked_example.py` (toy verification, no data file output
  — results are the printed comparison table above, captured in this
  notebook).
- **Results presented** — worked-example table shown inline above,
  matching script output shown to the researcher in chat.
- **QC gate** — `hit_rate()`'s per-gene counts match an independent
  hand tally of the same raw rows exactly (0 discrepancies, 2/2 genes);
  both positive-control genes score `significant_up`-only in their
  significant hits, matching known nitrogen-starvation-induction
  biology.
- **Decisions made this step** — none.
- **Advance rationale** — the shared method is verified correct on a
  known-biology worked example; ready for step 5 (analyze: apply
  `hit_rate()` across all 6 gene-set x experiment-set combinations for
  H1, H2, H3, and the positive-control sanity check).
