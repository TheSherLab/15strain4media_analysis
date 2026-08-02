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

## Decide-gate checklist

- **Outputs produced** — `scripts/01_compute_hit_rates.py` →
  `data/01_hit_rate_per_gene.csv` (346 rows), `data/01_hit_rate_summary.csv`
  (7 rows); `scripts/02_plot_hit_rate_comparison.py` →
  `figures/01_hit_rate_comparison.png`. Also updated
  `2_kg_selection/scripts/02_list_np_experiments.py` and re-ran it to add
  `table_scope`/`table_scope_detail` columns (redo 3 of step 2, additive
  only).
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
