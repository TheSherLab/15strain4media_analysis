# Step 5 — Analyze

## Context

Step 4 built and verified `hit_rate()`. This step applies it to the
group x experiment-set combinations defined in step 3, to produce the
numbers that step 6 will assess against H1, H2, and H3.

## Co-define

Proposed to the researcher: run `hit_rate()` on all 8 combinations,
produce a per-gene table, an aggregate summary table, and a comparison
figure. Agreed before running. (One of the 8 — `background_P` — was
removed after the first pass; see Decisions.)

## What I did

- `01_compute_hit_rates.py` — ran `hit_rate()` for: H1_N, H1_P,
  H2_Ngenes_Pexp, H2_Pgenes_Nexp, positive_control_N,
  positive_control_P, background_N. Output: 346 per-gene rows, 7
  group-level aggregate rows. `background_P` (H3 for phosphorus) was
  deliberately not computed — see Decisions.
- Added `table_scope` and `table_scope_detail` (the cutoff, when the
  source publication states one) columns to step 2's
  `02_np_experiments.csv`, at the researcher's request, so the
  significance-pre-filtering caveat below is traceable per experiment
  rather than only asserted in prose (see `2_kg_selection/notebook.md`'s
  updated experiment table).
- `02_plot_hit_rate_comparison.py` — grouped horizontal bar chart of all
  8 groups' up/down/not-significant rates, colored with the validated
  blue/red diverging pair (`validate_palette.js`, both checks pass) plus
  a neutral gray for "not significant," direct value labels as the
  required mitigation for the gray segment's low contrast.

While building this step's tables, a data-scope issue surfaced that
changes how the phosphorus-side numbers must be read (see Surprises) —
raised with and confirmed by the researcher before finalizing the
write-up below.

## Results

**Group summary** (`n_genes_in_set` = number of (strain, locus_tag)
instances eligible for that group, i.e. gene x strain, not distinct
gene names; `n_genes_tested` = how many of those actually had >=1 row
in the relevant experiments):

| group | n_genes_in_set | n_genes_tested | n_tests | n_up | n_down | n_not_significant | rate_up | rate_down | rate_significant |
|---|---|---|---|---|---|---|---|---|---|
| H1_N (N genes x N starvation, matched) | 93 | 46 | 403 | 111 | 5 | 287 | 0.275 | 0.012 | 0.288 |
| H1_P (P genes x P starvation, matched)* | 85 | 31 | 145 | 57 | 2 | 86 | 0.393 | 0.014 | 0.407 |
| H2_Ngenes_Pexp (N genes x P starvation, cross)* | 93 | 2 | 6 | 0 | 3 | 3 | 0.000 | 0.500 | 0.500 |
| H2_Pgenes_Nexp (P genes x N starvation, cross) | 85 | 45 | 395 | 27 | 35 | 333 | 0.068 | 0.089 | 0.157 |
| positive_control_N | 16 | 8 | 76 | 44 | 0 | 32 | 0.579 | 0.000 | 0.579 |
| positive_control_P* | 18 | 14 | 63 | 30 | 1 | 32 | 0.476 | 0.016 | 0.492 |
| background_N (H3) | 303 | 200 | 1596 | 78 | 158 | 1360 | 0.049 | 0.099 | 0.148 |

(* = phosphorus-experiment-side numbers; all pre-filtered by table_scope,
see Surprises — not directly comparable to nitrogen-side numbers.
`background_P`, H3's phosphorus counterpart, was not computed — see
Decisions.)

Full detail in `data/01_hit_rate_per_gene.csv` and
`data/01_hit_rate_summary.csv`. Figure: `figures/01_hit_rate_comparison.png`.

**Nitrogen side (valid H1/H3 comparison):** N target genes are
significant in 28.8% of tests (27.5% up, 1.2% down) versus background
genes at 14.8% (4.9% up, 9.9% down) in the same nitrogen experiments —
roughly double the overall hit rate, and a qualitatively different
direction profile: target genes are almost entirely upregulated when
significant, while background genes lean slightly toward downregulated.
The 4 N-side positive-control genes (`ntcA`, `glnA`, `amtB/amt1`,
`ureA`) score higher still: 57.9% significant, all upregulated, 0 down.
Top individual N-gene rates (>=3 tests): `glnA` (MED4) 85% up-rate,
`urtA` (MIT9313) 83%, `urtA` (MED4) 69%, `cynA` (MED4) 69%.

**Phosphorus side (numbers reported, no noise-check comparison):** P
target genes are significant in 40.7% of tests (39.3% up, 1.4% down);
the P-side positive controls score similarly (49.2%, 47.6% up). No
phosphorus background/H3 rate is reported (see Decisions) — with every
phosphorus experiment's gene table pre-filtered to already-significant
genes, a background sample drawn from the same tables would not be a
random sample of the genome, so it cannot serve as a noise check the
way the nitrogen background can. Top individual P-gene rates (>=3
tests): `gap3` (NATL2A) 100% up-rate, `phoB` (NATL2A) 83%,
`gap3`/`phoB`/`phoA` (MED4) 60% each, `pstS` (NATL2A) 57%.

**Cross-nutrient (H2):** P genes tested in nitrogen-starvation
experiments show 15.7% significant (6.8% up, 8.9% down, roughly
balanced direction) across 395 tests — a real-sized sample, sitting
between the N background rate (14.8%) and the N target-gene rate
(28.8%), closer to background. N genes tested in phosphorus-starvation
experiments returned only 6 total tests across 2 of 93 gene-instances —
too small to interpret (see Surprises).

## Surprises

- **Phosphorus experiments are all significance-pre-filtered
  (`table_scope`).** Checked via `list_experiments(experiment_ids=...)`:
  all 5 phosphorus experiments have `table_scope` of `significant_only`
  or `filtered_subset` (34-177 genes each — only genes that already
  passed a significance or fold-change cutoff in the source
  publication's supplementary table are present in the KG at all; genes
  never significant are absent, not recorded as "not significant"). By
  contrast, 4 of 5 nitrogen experiments are `all_detected_genes`
  (1,400-2,200 genes each, no significance threshold on inclusion); the
  5th nitrogen experiment has a milder, non-significance filter ("top
  50% of genes by expression level"). Practical effect: every
  phosphorus-side rate in this step (H1_P, H2_Ngenes_Pexp,
  positive_control_P) is inflated relative to what it would be against
  a fully-measured gene population, and — critically — a phosphorus
  background/H3 set cannot serve as a valid noise check, since it would
  be drawn from these same pre-filtered tables. Researcher confirmed:
  report the phosphorus target/positive-control numbers with this
  caveat attached, but do not compute a phosphorus background set at
  all (see Decisions) — only the nitrogen-side H3 comparison is valid.
- **H2_Ngenes_Pexp has almost no data (n=6).** Only 2 of the ~23-25
  nitrogen-annotated gene-instances per relevant strain had any row at
  all in the phosphorus experiments -- a direct consequence of the same
  table_scope issue: phosphorus tables are so narrow (34-177 genes) that
  most nitrogen-annotated genes simply never appear in them, regardless
  of whether they'd have been significant. This group's 50% rate is an
  artifact of a 6-row sample and should not be read as evidence for or
  against H2 on the nitrogen-genes-under-phosphorus-starvation side.

## Decisions

- **2026-08-02** — Report all phosphorus target/positive-control
  numbers (not drop phosphorus from the analysis), with the
  table_scope caveat stated prominently.
- **2026-08-02** — `background_P` (H3's phosphorus side) is **not
  computed at all**, superseding an earlier draft of this step that
  computed and reported it as "invalid." Researcher's rationale: since
  the phosphorus background pool would be drawn from the same
  pre-filtered tables as the target genes, reporting a number and
  labeling it invalid is worse than not reporting it — the number
  itself has no interpretable meaning here, not just a weaker one.
  Instead, the underlying `table_scope`/cutoff evidence is now carried
  as data (see What I did) so the caveat is verifiable, not just
  asserted.

## Addendum (2026-08-02) — log2FC magnitude figures

Everything above (the hit-rate figure and tables) used only the
categorical `expression_status` call (up/down/not-significant), never
the actual fold-change magnitude. The researcher asked what metric
underlies that call (it's `log2fc` + `padj` per row, direction from the
sign of log2fc, significance from padj against the study's own
threshold) and asked for figures showing magnitude, split by omics
platform (RNA-seq/proteomics/microarray are not directly comparable in
magnitude per this methodology's statistical-rigor rule).

**What I did:**
- `03_pull_log2fc_data.py` — re-queried the same 5 groups (H1_N, H1_P,
  positive_control_N, positive_control_P, background_N; background_P
  still excluded) keeping `log2fc`, `padj`, and `omics_type` per row.
  Row counts match `01_compute_hit_rates.py` exactly (403/145/76/63/1596).
- `04_plot_log2fc_distribution.py` → `figures/02_log2fc_distribution.png`
  — box+strip plot of log2FC among significant hits only, per group,
  colored by direction (blue=up, red=down).
- `05_plot_log2fc_by_platform.py` → `figures/03_log2fc_by_platform.png`
  — same, faceted into 3 panels (RNASEQ/PROTEOMICS/MICROARRAY) instead
  of pooled, per the cross-platform-magnitude caveat.

**Results:** the direction pattern established in step 6 is visibly
sharper once magnitude is shown, and holds within every platform, not
just in aggregate. Nitrogen target genes and positive controls cluster
tightly positive (median log2FC 2.3-2.4, narrow spread, almost no
negative points). The nitrogen background set is the mirror image —
its significant hits are predominantly *negative* (median magnitude
1.7, but visibly down-skewed rather than up-skewed), which is the same
asymmetry noted qualitatively in step 6 now visible as an actual
fold-change pattern, not just a rate. Phosphorus target genes and
positive controls also skew positive, with a wider spread and higher
ceiling (median log2FC 2.6-3.1) — including 3 known outlier genes
(`PMM0707`, `PMM0708`=`phoA`, `PMM1416`, all MED4 phosphate-starvation
microarray, log2FC 30-162) flagged as likely artifacts by the prior
triage analysis and not filtered here either; the figures clip the
axis and note the count rather than hide or silently include them.

**Surprise:** the 3 known outliers from the prior analysis
(`2026-07-13-n_p_genes_in_vivo_experiments`) reappeared here
independently — good cross-analysis consistency, and confirms they're
a property of the underlying KG data (MED4 phosphate microarray study,
`10.1073/pnas.0601301103`), not an artifact of either analysis's own
processing.

**QC gate (addendum):** row counts from `03_pull_log2fc_data.py`
reconcile exactly with `01_hit_rate_summary.csv`'s `n_tests` per group;
manually cross-checked the 3 outlier loci against the prior triage
analysis's own flagged-outlier note (same loci, same experiment, same
magnitude range).

This addendum is additive: it does not change any number in
`01_hit_rate_summary.csv` or `01_hit_rate_per_gene.csv`, which step 6
already evaluated and locked, so step 6's conclusions stand unchanged.
It gets its own commit (step 5 was already closed and step 6 already
built on it) rather than amending the prior step-5 commit.

## Addendum 2 (2026-08-02) — nitrogen volcano plot, all 53 genes

Researcher asked: since the nitrogen result is established as real (not
background noise), show a volcano plot of all 53 usable genes'
response to nitrogen starvation, colored by each gene's *original*
N/P-acquisition annotation from the researcher's source spreadsheet
(not by what this analysis found).

**What I did:**
- `06_pull_nitrogen_volcano_data.py` — pulled DE data for all 53 genes
  (both N- and P-annotated) restricted to the 5 nitrogen experiments;
  798 rows = 403 (H1_N) + 395 (H2_Pgenes_Nexp), reconciling exactly with
  step 5's earlier totals.
- `07_plot_nitrogen_volcano.py` → `figures/04_nitrogen_volcano.png` —
  volcano plot (log2FC vs. -log10 padj), faceted by platform (padj
  granularity differs sharply: RNASEQ/PROTEOMICS are continuous,
  MICROARRAY is coarsely discretized to ~0.01/1.0 in this KG build, so
  pooling platforms would be misleading here too). padj=0 rows (16,
  all RNASEQ) are floored to one order of magnitude below the smallest
  observed nonzero padj (1.4e-18) so they plot as most-significant, not
  tied at an arbitrary constant. Top 3 genes per panel by |log2FC|
  among significant hits are labeled.

**Results:** the RNASEQ and PROTEOMICS panels visibly show blue
(N-annotated) points concentrated in the upper-right (significant,
upregulated) — the same pattern established quantitatively in step 6,
now visible gene-by-gene. Three P-annotated genes stand out as notable
exceptions worth flagging individually (not washed out by the
aggregate H2 "no effect" result, which describes the P-gene set on
average, not every gene in it):
- **`pstS`** is strongly significantly *down* in MED4 under nitrogen
  starvation across all 3 platforms it's measured in (RNASEQ log2FC
  -3.1 to -4.5, padj as low as 1.4e-17; MICROARRAY log2FC -1.2 to -2.9)
  — but strongly significantly *up* in MIT9313 (MICROARRAY log2FC
  1.3-2.1). A strain-specific divergence in the same gene, opposite
  directions.
- **`phnD`** is consistently, strongly down in MED4 across all 3
  platforms (log2FC -2.1 to -3.2, padj down to 1.4e-8).
- **`PMM722`** is strongly up in MED4 (RNASEQ log2FC 4.3, padj~0;
  PROTEOMICS log2FC 1.9, padj 0.005).

These are `[interpretation]`-flagged candidates for follow-up, not a
revision of H2's aggregate conclusion — H2 tested and did not find a
P-gene-set-wide effect, and individual strongly-responding genes within
a set that shows no *average* effect is expected, not contradictory.

**QC gate (addendum 2):** row count (798) reconciles with H1_N (403) +
H2_Pgenes_Nexp (395) from `01_hit_rate_summary.csv`; spot-checked that
`cynA`/`glnA`/`ntcA` plot as N-annotated (blue) and `pstS`/`phnD`/
`PMM722` as P-annotated (orange) against `2_kg_selection/data/03_gene_loci.csv`'s
`n_or_p` column, confirming the color-mapping code isn't reversed.

## Decide-gate checklist

- **Outputs produced** — `scripts/01_compute_hit_rates.py` →
  `data/01_hit_rate_per_gene.csv` (346 rows), `data/01_hit_rate_summary.csv`
  (7 rows); `scripts/02_plot_hit_rate_comparison.py` →
  `figures/01_hit_rate_comparison.png`. Also updated
  `2_kg_selection/scripts/02_list_np_experiments.py` and re-ran it to add
  `table_scope`/`table_scope_detail` columns (redo 3 of step 2, additive
  only). Addendum: `scripts/03_pull_log2fc_data.py` →
  `data/03_log2fc_raw.csv` (2283 rows); `scripts/04_plot_log2fc_distribution.py`
  → `figures/02_log2fc_distribution.png`; `scripts/05_plot_log2fc_by_platform.py`
  → `figures/03_log2fc_by_platform.png`. Addendum 2:
  `scripts/06_pull_nitrogen_volcano_data.py` → `data/06_nitrogen_volcano_data.csv`
  (798 rows); `scripts/07_plot_nitrogen_volcano.py` →
  `figures/04_nitrogen_volcano.png`.
- **Results presented** — group summary table and per-gene highlights
  shown inline above, matching what was shown to the researcher in chat;
  figure included.
- **QC gate** — confirmed `hit_rate()`'s underlying data source
  (`table_scope`) for every one of the 10 experiments before trusting
  the aggregate rates; confirmed the H2_Ngenes_Pexp small-sample issue
  traces to the same root cause rather than a separate bug; palette
  validated via `validate_palette.js` (blue/red pair: all checks pass).
- **Decisions made this step** — see Decisions above.
- **Advance rationale** — all 3 hypotheses now have numbers to assess
  against, with the phosphorus-side and H2-nitrogen-genes caveats
  explicit; ready for step 6 (evaluate: assess against the framing,
  harvest caveats, finalize the paper).
