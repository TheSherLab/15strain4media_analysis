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

> **Re-evaluated 2026-09-08.** Step 2 corrected 6 gene loci (4× `phoE`,
> `unkP2`, `urtA`; `ptrA`-NATL2A added), step 5 moved Lin to 59h and
> reclassified Martiny table-absent genes. This step's scripts (`01`,
> `02`) were re-run. All tables and evaluation text below are post-change;
> the prior values are shown alongside where they moved. `gaps_and_friction.md`
> (2026-09-08) has the methodology record.

**Positive controls, per gene** (matched-nutrient tests only; full table:
`data/01_positive_control_by_gene.csv`):

| Nutrient | Gene | Tests | Up | Down | Not sig. | % significant |
|---|---|---|---|---|---|---|
| N | ntcA | 5 | 5 | 0 | 0 | 100.0% |
| N | glnA | 5 | 5 | 0 | 0 | 100.0% |
| N | ureA | 5 | 4 | 0 | 1 | 80.0% |
| N | amtB/amt1 | 5 | 2 | 0 | 3 | 40.0% |
| P | pstS | 5 | 4 | 1 | 0 | 100.0% |
| P | phoA | 3 | 3 | 0 | 0 | 100.0% |
| P | phoB | 3 | 3 | 0 | 0 | 100.0% |
| P | phoR | 3 | 2 | 0 | 1 | 66.7% |

Pooled: N controls 16/20 = 80.0% (unchanged); P controls 13/14 = 92.9%
(was 84.6%, 11/13). `phoA` moved 66.7% → 100% and `phoR` 50% → 66.7%,
both because Lin now sits at 59h where the whole pho regulon is engaged
(`phoA` and `phoR` are `not_significant` at Lin 46h, `significant_up` at
59h). Every nitrogen control that is ever significant is significant
upregulated. `pstS` is still the only phosphorus control with a
downregulated hit (1 of 5).

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
the observed matched-nitrogen target rate (44.2%). (The background pool is
unchanged by every 2026-09-08 fix — it is built from the nitrogen
`all_detected_genes` experiments and excludes target genes by name, so
the `urtA` locus swap and the phosphorus-side changes do not touch it.)

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

The observed matched-nitrogen rate (44.2%) exceeds every one of the
10,000 random size-matched draws — not just the mean.

**Stability check — matched-N with positive controls excluded** (full
output: `data/04_stability_exclude_controls.txt`; post-2026-09-08):

| | Genes | Tests | % significant | Bootstrap null mean | Bootstrap p | Fisher's OR | Fisher's p |
|---|---|---|---|---|---|---|---|
| Full matched-N | 29 | 104 | 44.2% | 18.8% | <0.0001 | 3.39 | 4.98e-09 |
| Excluding 4 positive controls | 25 | 84 | 35.7% | 18.7% | <0.0001 | 2.37 | 3.77e-04 |

(Prior: full 41.3%, OR 3.01; no-controls 32.1%, OR 2.02, bootstrap
p=0.0012.) The signal weakens without the 4 well-known markers but
remains statistically distinguishable from background by both tests, and
slightly firmer than before — the matched-nitrogen response is not solely
carried by `ntcA`/`glnA`/`ureA`/`amtB`, and the real `urtA` locus (step 2
fix) contributes to the non-control set.

## Evaluation against the step-3 framing

**H1 (matched response) — held for both nutrients.** Nitrogen-annotated
genes respond to nitrogen starvation (44.2% significant, 100% of that
upregulated); phosphorus-annotated genes respond to phosphorus starvation
(60.6% significant, mostly upregulated). `[interpretation]` The phosphorus
rate (60.6%, down from a pre-2026-09-08 75.7%) is now measured against a
fairer denominator: the Martiny microarray's genome-present-but-table-
absent genes are counted as tested-not-significant, so its contribution is
no longer all-significant. It is still above the nitrogen rate, and still
not directly comparable to it — the other 3 phosphorus tables (Lin,
Fuszard ×2) remain pre-filtered/detection-limited, which is why H3 stays
nitrogen-only. The `phoE` fix (step 2) added a canonical pho-regulon
responder that had been pointing at the wrong porin.

**H2 (cross-nutrient response) — not supported as a general stress
response; now demonstrated in both directions.** Phosphorus genes under
nitrogen starvation: 20.7% significant, split roughly evenly between up
(10) and down (13). Nitrogen genes under phosphorus starvation: **48
tests, 4.2% significant** (was 1 test) — the Martiny reclassification
(step 5) makes the microarray's nitrogen-gene cross-tests visible, and
they are almost all flat. `[interpretation]` A low rate with no
consistent direction in *both* cross directions is what a nutrient-
specific response predicts and a general-stress response would not. The
"no general cross-response" conclusion no longer rests on a data vacuum in
the N-genes-under-P direction.

**H3 (nitrogen noise/bootstrap check) — held, and robust to the positive
controls' influence.** 44.2% exceeds the full 10,000-draw bootstrap null
distribution's maximum (32.6%), p<0.0001; Fisher's exact agrees (OR=3.39,
p=4.98e-09) — slightly firmer than before (OR 3.01) because the real
`urtA` locus is more N-responsive than the `urtE` locus step 2 had been
using. The stability check above shows this is not an artifact of the 4
positive-control genes: the other 25 genes alone still clear background at
35.7% vs. a null mean of 18.7% (bootstrap p<0.0001,
Fisher's OR=2.37, p=3.77e-04).

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
- The observed matched-nitrogen rate (44.2%) exceeds not just the
  bootstrap's mean but its maximum across all 10,000 draws (32.6%) —
  worth stating plainly since "p<0.0001" alone doesn't convey how far
  outside the null distribution the observed value sits.

## Decisions

**2026-08-12 — Ran a stability check excluding the 4 nitrogen positive
controls from H1/H3, co-defined with the researcher.** Motivated by the
controls' much higher individual rate (80.0%) than the full matched-N set
raising the question of whether the significant background-beating result
was mostly driven by 4 already-known marker genes. Confirmed it is not:
the remaining 25 genes still clear background significantly, just less
strongly.

**2026-09-08 — Re-ran step 6 (`01`, `02`) after the step-2 gene-identity
reopen and the step-5 Lin/Martiny changes.** Not a new judgment call —
step 6 consumes steps 2 and 5, so it re-runs when they change. Numbers
moved (P positive controls 84.6% → 92.9%; H3 Fisher OR 3.01 → 3.39;
stability-check no-controls 32.1% → 35.7%); every conclusion holds and the
nitrogen result is firmer. `02_stability_exclude_controls.py` still uses
seed 42 / 10,000 iterations. The `03`/`04` post-closure scripts (2026-09-06
phosphorus-coverage checks) are unaffected and not re-run.

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
  matches; P: 12 up/1 down/1 not-sig of 14 tests = 92.9%, matches
  post-2026-09-08). Stability-check bootstrap reused the same seed (42)
  and iteration count (10,000) as step 5's H3 for direct comparability.
- **Decisions made this step:** stability check excluding positive
  controls, co-defined and run (2026-08-12, above).
- **Advance rationale:** all three hypotheses are evaluated against the
  step-3 framing with the full positive/negative-control detail
  presented, the one result-triggered stability check is done, and
  caveats are harvested — ready to finalize `paper.md`'s Discussion.

## Post-closure addendum (2026-09-06)

After this step was closed and committed, the researcher asked why
specific target genes (cynD, pipX, urtD in MED4; phoE in MIT9313) show
"no data" for particular experiments in Figure 1, and whether the
phosphorus experiments' narrow coverage was a property of the source
tables or a KG gap. Re-verified directly against the live KG (see
`scripts/03_post_closure_missing_gene_checks.py`,
`scripts/04_post_closure_phosphorus_coverage.py`) rather than trusting
this notebook's own prior summary text. Full write-up:
`gaps_and_friction.md` (2026-09-06 entry). Two distinct mechanisms
confirmed:

- **Read et al. 2017 MED4 RNA-seq's "top 50% by expression" filter is
  genome-wide, not per-timepoint** — cynD, pipX, urtD are absent from all
  3 of the experiment's timepoints (3h, 12h, 24h), not just the 24h
  timepoint this analysis chose as the representative starvation point.
- **phoE (locus `PMT_2631`) in MIT9313 has zero DE rows in any MIT9313
  experiment in the entire KG** — a platform/probe coverage gap, distinct
  from Read's declared filter, since Tolonen's MIT9313 table (the only
  other in-scope MIT9313 experiment) is a broad, largely unfiltered
  sample (76% of the genome).
- **Quantified genome coverage for all 5 phosphorus experiments**
  (`data/04_phosphorus_genome_coverage.csv`): 1.5-6.0% of each strain's
  genome, regardless of method (RNA-seq, proteomics, microarray) —
  versus 76-86% for the nitrogen microarray tables on the same platform
  type. Confirms the phosphorus side's narrow coverage is a property of
  each source publication's own filtered supplementary table, not a
  technology limitation.

No results or conclusions change — this quantifies and confirms caveats
already present in this notebook and `paper.md`'s Discussion, rather than
overturning them. `paper.md`'s phosphorus-background caveat updated with
the exact percentages.
