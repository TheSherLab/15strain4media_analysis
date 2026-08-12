# Step 6 — Evaluate

## Context

Step 5 computed all three hypotheses on the full 61-gene, 10-experiment
dataset. This step assesses those results against the step-3 framing,
digs into what the positive and negative controls actually showed gene by
gene (not just as pooled rates), runs one stability check the pooled
numbers raised, harvests caveats, and writes the paper's Discussion.

## What I did

**Positive-control detail** (`scripts/01_control_detail.py`): pulled every
(gene x experiment) test for the 8 positive-control genes — 4 nitrogen
(`ntcA`, `glnA`, `amtB/amt1`, `ureA`), 4 phosphorus (`pstS`, `phoA`,
`phoB`, `phoR`) — from the same matched-nutrient scope step 5 used, broken
out per gene rather than pooled into the step-5 summary numbers (80.0%
nitrogen, 84.6% phosphorus).

**Negative-background detail** (same script): broke the nitrogen
background pool's "% significant" rate out per experiment (it was only
reported pooled/mean before), and pulled percentiles (not just mean/std)
of the H3 bootstrap's 10,000-draw null distribution.

**Stability check** (`scripts/02_stability_exclude_controls.py`,
co-defined with the researcher before running): the 4 nitrogen positive
controls individually score much higher (80.0%) than the full 29-gene
matched-N set (41.3%) — reran `hit_rate()` and `bootstrap_pvalue()` on
matched-N with the 4 controls excluded, same background pools, same seed
and iteration count as step 5's H3, to check whether the matched-nitrogen
signal depends on those 4 well-established genes.

## Results

**Positive controls, per gene** (matched-nutrient tests only; full table:
`data/01_positive_control_by_gene.csv`):

| Nutrient | Gene | Tests | Up | Down | Not sig. | % significant |
|---|---|---|---|---|---|---|
| N | ntcA | 5 | 5 | 0 | 0 | 100.0% |
| N | glnA | 5 | 5 | 0 | 0 | 100.0% |
| N | ureA | 5 | 4 | 0 | 1 | 80.0% |
| N | amtB/amt1 | 5 | 2 | 0 | 3 | 40.0% |
| P | pstS | 5 | 4 | 1 | 0 | 100.0% |
| P | phoB | 3 | 3 | 0 | 0 | 100.0% |
| P | phoA | 3 | 2 | 0 | 1 | 66.7% |
| P | phoR | 2 | 1 | 0 | 1 | 50.0% |

Every nitrogen control that is ever significant is significant upregulated
— never downregulated, including `amtB/amt1`, the weakest of the four.
`pstS` is the only phosphorus control with a downregulated hit (1 of 5
tests); the other 3 phosphorus controls are up-only when significant.
`phoR` has only 2 tests in the dataset (smallest evidence base of the 8).

**Negative background pool, per experiment** (nitrogen only; full table:
`data/02_negative_background_by_experiment.csv`):

| Experiment | Background genes | % significant |
|---|---|---|
| Tolonen microarray (MIT9313) | 2,196 | 3.0% |
| Weissberg proteomics (MED4) | 1,387 | 9.1% |
| Tolonen microarray (MED4) | 1,656 | 17.6% |
| Weissberg RNA-seq (MED4) | 1,808 | 47.2% |

The pooled mean reported in step 5 (18.8%) sits inside this range, but the
range itself is wide — a nearly 16-fold spread from lowest to highest.
The Weissberg RNA-seq background rate (47.2%) is on its own higher than
the observed matched-nitrogen target rate (41.3%).

**Bootstrap null distribution, with percentiles** (full data:
`../5_analyze/data/02_bootstrap_null_distribution.csv`; percentiles:
`data/03_bootstrap_null_percentiles.txt`):

| Percentile | Null hit rate |
|---|---|
| P1 | 10.5% |
| P25 | 16.3% |
| P50 (median) | 18.6% |
| P75 | 20.9% |
| P95 | 25.6% |
| P99 | 27.9% |
| P100 (max of 10,000 draws) | 32.6% |

The observed matched-nitrogen rate (41.3%) exceeds every one of the
10,000 random size-matched draws — not just the mean.

**Stability check — matched-N with positive controls excluded** (full
output: `data/04_stability_exclude_controls.txt`):

| | Genes | Tests | % significant | Bootstrap null mean | Bootstrap p | Fisher's OR | Fisher's p |
|---|---|---|---|---|---|---|---|
| Full matched-N | 29 | 104 | 41.3% | 18.8% | <0.0001 | 3.01 | 2.15e-07 |
| Excluding 4 positive controls | 25 | 84 | 32.1% | 18.7% | 0.0012 | 2.02 | 4.70e-03 |

The signal weakens without the 4 well-known markers but remains
statistically distinguishable from background by both tests — the
matched-nitrogen response is not solely carried by `ntcA`/`glnA`/`ureA`/
`amtB`.

## Evaluation against the step-3 framing

**H1 (matched response) — held for both nutrients.** Nitrogen-annotated
genes respond to nitrogen starvation (41.3% significant, 100% of that
upregulated); phosphorus-annotated genes respond to phosphorus starvation
(75.7% significant, mostly upregulated). `[interpretation]` The higher
phosphorus rate is not necessarily a stronger biological response: all 5
phosphorus source tables are pre-filtered toward genes the source
publication already found interesting (`significant_only`/
`filtered_subset`, step 3), while the 4 unbiased nitrogen tables report
every gene tested — so the phosphorus percentage is measured against a
smaller, pre-curated denominator, which can mechanically inflate the rate
independent of biology. This is exactly why H3 was restricted to nitrogen
only.

**H2 (cross-nutrient response) — not supported as a general stress
response, though only demonstrated in one direction.** Phosphorus genes
under nitrogen starvation: 20.2% significant, split roughly evenly
between up (9) and down (13) — both a lower rate and a different
direction-pattern than the matched-nitrogen result (41.3%, uniformly up).
`[interpretation]` Lower rate with no consistent direction is what a
nutrient-specific response predicts. Nitrogen genes under phosphorus
starvation: 1 test total, too sparse to interpret — a direct consequence
of the phosphorus tables' narrow pre-filtered coverage (step 3). The
"no general cross-response" conclusion therefore rests on real evidence
in only one of the two cross-nutrient directions.

**H3 (nitrogen noise/bootstrap check) — held, and robust to the positive
controls' influence.** 41.3% exceeds the full 10,000-draw bootstrap null
distribution's maximum (32.6%), p<0.0001; Fisher's exact agrees (OR=3.01,
p=2.15e-07). The stability check above shows this is not an artifact of
the 4 positive-control genes: the other 25 genes alone still clear
background at 32.1% vs. a null mean of 18.7% (bootstrap p=0.0012,
Fisher's OR=2.02, p=4.70e-03).

No formal preregistration was written for this analysis (step 3's framing
stated hypotheses in prose, not point predictions with thresholds) — this
evaluation compares the three hypotheses' direction and the H3
significance test to what step 3 said would count as support, not to a
locked numeric prediction.

## Surprises

- The nitrogen background "noise" rate is not one number — it ranges from
  3.0% (Tolonen MIT9313 microarray) to 47.2% (Weissberg RNA-seq), a
  ~16-fold spread across the 4 experiments it's built from. The pooled
  18.8% mean reported in step 5 is real but hides this heterogeneity.
  `[interpretation]` The Weissberg RNA-seq experiment's high baseline
  rate likely reflects a broad transcriptional response to nitrogen
  starvation affecting many genes genome-wide, or a comparatively liberal
  DESeq2 call in that dataset specifically — not "noise" in the everyday
  sense of the word. The bootstrap is not misled by this because it draws
  each iteration's sample from the same per-experiment pools the real
  target genes were tested in (not one pooled rate), but the term
  "background/noise rate" oversimplifies what is really a mix of
  per-platform baseline rates.
- The observed matched-nitrogen rate (41.3%) exceeds not just the
  bootstrap's mean but its maximum across all 10,000 draws (32.6%) —
  worth stating plainly since "p<0.0001" alone doesn't convey how far
  outside the null distribution the observed value sits.

## Decisions

**2026-08-12 — Ran a stability check excluding the 4 nitrogen positive
controls from H1/H3, co-defined with the researcher.** Motivated by the
controls' much higher individual rate (80.0%) than the full matched-N set
(41.3%) raising the question of whether the significant background-beating
result was mostly driven by 4 already-known marker genes. Confirmed it is
not: the remaining 25 genes still clear background significantly, just
less strongly.

## Decide-gate checklist

- **Outputs produced:** `scripts/01_control_detail.py`,
  `scripts/02_stability_exclude_controls.py`;
  `data/01_positive_control_by_gene.csv`,
  `data/02_negative_background_by_experiment.csv`,
  `data/03_bootstrap_null_percentiles.txt`,
  `data/04_stability_exclude_controls.txt`.
- **Results presented:** positive-control per-gene table, negative-
  background per-experiment table, bootstrap-null percentile table,
  stability-check table (all above).
- **QC gate:** per-gene positive-control totals cross-checked against
  step 5's pooled numbers (N: 16 up/0 down/4 not-sig of 20 tests = 80.0%,
  matches; P: 10 up/1 down/2 not-sig of 13 tests = 84.6%, matches).
  Stability-check bootstrap reused the same seed (42) and iteration count
  (10,000) as step 5's H3 for direct comparability.
- **Decisions made this step:** stability check excluding positive
  controls, co-defined and run (2026-08-12, above).
- **Advance rationale:** all three hypotheses are evaluated against the
  step-3 framing with the full positive/negative-control detail
  presented, the one result-triggered stability check is done, and
  caveats are harvested — ready to finalize `paper.md`'s Discussion.
