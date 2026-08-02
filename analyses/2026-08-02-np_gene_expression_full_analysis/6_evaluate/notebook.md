# Step 6 — Evaluate

## Context

Step 5 produced significant-hit rates for all target/control gene sets
against nitrogen and phosphorus starvation. This step assesses those
numbers against the three hypotheses locked in step 3, adds the one
formal statistical test the rates need (Fisher's exact, nitrogen side
only), harvests caveats, and finalizes `paper.md`.

## Co-define

Proposed to the researcher: run a Fisher's exact test on the one valid
comparison (nitrogen target genes vs. nitrogen background), assess
H1/H2/H3 in prose, harvest caveats from everything logged across steps,
finalize paper.md's Discussion and References. Agreed before proceeding.
Mid-step, the researcher asked exactly what the test covered and
pointed out H2's nitrogen-side cross-nutrient comparison
(`H2_Pgenes_Nexp` vs. `background_N`) is also a valid, unbiased
comparison that hadn't been statistically tested — added as a second
test (see What I did).

## What I did

- `01_fisher_test.py` — two Fisher's exact tests, both against
  `background_N` (the only unbiased background this KG build has):
  1. `H1_N` vs. `background_N` — matched nitrogen genes (also serves as
     the H3 noise check for nitrogen).
  2. `H2_Pgenes_Nexp` vs. `background_N` — phosphorus-annotated genes
     under nitrogen starvation (cross-nutrient).
  No phosphorus-background test is possible — see step 5's decision not
  to compute a phosphorus background set at all; `H1_P` and
  `H2_Ngenes_Pexp` are reported descriptively only, never statistically
  tested.

## Results

**Fisher's exact tests (both vs. nitrogen background):**

| comparison | target significant | target total | background significant | background total | fold enrichment | odds ratio | p-value |
|---|---|---|---|---|---|---|---|
| H1_N (matched N genes) | 116 | 403 | 236 | 1596 | 1.95x | 2.33 | 3.2e-10 |
| H2_Pgenes_Nexp (P genes, cross) | 62 | 395 | 236 | 1596 | 1.06x | 1.07 | 0.637 |

H1_N: target rate 28.8% vs. background rate 14.8% — strongly
significant (per this methodology's significance language guide).
H2_Pgenes_Nexp: target rate 15.7% vs. the same 14.8% background —
statistically indistinguishable from background (p = 0.64). Full
detail: `data/01_fisher_test_nitrogen.csv`.

**Caveat on this p-value:** the 403 vs. 1596 "tests" are (gene x
experiment x timepoint) rows, not independent biological samples — a
single gene contributes multiple correlated rows across timepoints and
experiments. The p-value should be read as strong evidence a real
difference exists, not taken as a literal sample-size-403 precision
estimate. The 1.95x fold enrichment and the direction split (below) are
the more robust take-aways than the p-value's exact magnitude.

**H1 (matched, directional) — nitrogen: supported.** Target genes are
significantly enriched for a hit (1.95x, p = 3.2e-10) and, more
tellingly, are almost entirely upregulated when significant (27.5 of
28.8 percentage points, i.e. 96% of hits are "up"). The 4 positive
controls go further: 57.9% significant, 100% of that up, 0% down. This
matches H1's prediction exactly: nitrogen-annotated genes go up under
nitrogen starvation, more than chance, and specifically in the
predicted direction.

**H1 (matched, directional) — phosphorus: directionally consistent,
not independently confirmed.** Phosphorus target genes show the same
qualitative pattern (39.3% up, 1.4% down — a 28:1 up:down ratio,
comparable to nitrogen's) and positive controls are similarly
up-skewed (47.6% up, 1.6% down). No statistical test is reported here
(see caveats) because there is no valid phosphorus background rate to
test against. `[interpretation]`: the up:down asymmetry is still mildly
informative on its own — the phosphorus papers' inclusion filters
select on fold-change magnitude in either direction (e.g. ">1.6 or
<0.6", a symmetric cutoff), so if the high phosphorus hit rate were
purely a table_scope artifact with no real directional biology, a
roughly even up/down split among the pre-selected genes would be a
more likely outcome than the 28:1 skew actually observed. This is
suggestive, not proof — a properly powered comparison would need an
unbiased phosphorus dataset this KG build does not contain.

**H2 (cross-nutrient) — phosphorus genes under nitrogen starvation:
not supported.** 15.7% significant across a real sample (395 tests) —
tested directly against nitrogen background (14.8%): fold enrichment
1.06x, p = 0.64, statistically indistinguishable from background. The
direction is also roughly balanced (6.8% up, 8.9% down), unlike H1's
strong directional signature. This is a clean negative result, not
just "weaker evidence" — nitrogen starvation does not measurably
perturb phosphorus-acquisition genes beyond what background genes show.

**H2 (cross-nutrient) — nitrogen genes under phosphorus starvation:
untestable.** Only 6 total tests (2 of 93 gene-instances had any
phosphorus-experiment data at all) — this KG build cannot answer this
side of H2 given the phosphorus experiments' narrow, pre-filtered gene
coverage.

**H3 (noise / method-limitation check) — nitrogen: refuted (the
signal is not just noise).** The nitrogen target-gene response is both
statistically distinguishable from background (p = 3.2e-10) and
directionally distinct: target genes are upregulated when significant
96% of the time, while the background set itself leans the *other*
way (33% up / 67% down among its significant hits — `[interpretation]`,
plausibly a general transcriptional-shutdown response among genes not
specifically tied to nitrogen acquisition, though this KG build wasn't
designed to test that directly). A pure significance-calling artifact
would not be expected to produce this directional asymmetry.

**H3 — phosphorus: untestable, not merely "weaker."** No phosphorus
background set was computed (step 5 decision) because it would share
the same pre-filtering as the target genes — there is no dataset in
this KG build that could answer whether the phosphorus signal is noise
or biology.

## Caveats (harvested)

- **Phosphorus data is uniformly pre-filtered to already-significant
  genes** (`table_scope` `significant_only`/`filtered_subset` in all 5
  phosphorus experiments) — inflates absolute phosphorus rates and
  blocks a phosphorus noise check entirely. This is the single largest
  limitation of this analysis; see `2_kg_selection/data/02_np_experiments.csv`
  for the cutoff recorded per experiment.
- **MIT9301 has zero usable data** in this KG build (its only
  phosphorus experiments are metabolomics-only; it has no nitrogen
  experiments at all) — dropped entirely, and `ptxB/phnD2` loses its
  only evidence as a consequence.
- **SS120 (CCMP1375) excluded** — not part of the researcher's strain
  scope, confirmed 2026-08-02.
- **4 "growth on alternate N source" experiments excluded** from the
  nitrogen scope (cyanate/urea/nitrite vs. N-replete) because they test
  a different question (alternate-source utilization, not starvation).
  This likely undercounts the nitrogen-side response for `cynA/B/D/S`,
  `ureA-G`, and `nirA` specifically, since those experiments are
  arguably the most on-target evidence for exactly those genes'
  designed function.
- **H2's nitrogen-genes-under-phosphorus-starvation side is
  untestable** with only 6 tests, a direct consequence of the
  phosphorus table_scope narrowness, not a separate data gap.
- **Uneven strain representation.** MIT9312 and NATL2A's phosphorus
  evidence comes from a targeted iTRAQ proteomics panel (19-25 genes
  tested), far narrower than MED4/MIT9313's microarray/RNA-seq coverage
  (1,400+ genes) — hit rates pool across strains, so strains with
  broader coverage dominate the aggregate numbers.
- **Positive-control genes were selected by researcher/Claude domain
  knowledge** (`[interpretation]`, canonical literature markers), not
  independently validated within this KG — a reasonable sanity check,
  but not itself evidence for H1.
- **Starvation-induction method** (centrifugation vs. dilution/transfer)
  is not recorded in this KG build for any of the 10 experiments —
  flagged, not resolved (researcher decision: not needed for this
  analysis).
- **This analysis pools across timepoints and experiments** within
  each group; it does not test per-strain or per-timepoint consistency
  separately (the prior triage analysis's gene x strain matrices are a
  related resource for that finer-grained view, not redone here).

## Decisions

None this step (evaluation only; no new forks).

## Decide-gate checklist

- **Outputs produced** — `scripts/01_fisher_test.py` →
  `data/01_fisher_test_nitrogen.csv`.
- **Results presented** — Fisher's exact test result and per-hypothesis
  assessment shown inline above, matching what was shown to the
  researcher in chat.
- **QC gate** — confirmed both 2x2 tables' row/column sums reconcile
  with step 5's `01_hit_rate_summary.csv` before running the tests;
  confirmed no phosphorus-side statistical test is reported, consistent
  with step 5's decision not to compute a phosphorus background set.
- **Decisions made this step** — none.
- **Advance rationale** — all three hypotheses have been assessed
  against the data with their appropriate caveats; `paper.md` is
  finalized (Discussion + References, below); the analysis is complete.
