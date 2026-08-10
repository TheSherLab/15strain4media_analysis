# Step 3 — Analysis framing

## Context

Step 2 established the usable experiment set (10 experiments, 4 strains,
6 publications) and gene universe (60 of 92 genes with actual DE evidence).
This step operationalizes the three hypotheses from step 1, documents
exactly what "significant" means in each experiment, and materializes the
positive and negative control gene sets.

## What I did

**Significance criteria** (`scripts/01_significance_criteria.py`): for
each of the 10 in-scope experiments, pulled every DE row
(`differential_expression_by_gene`) and computed the empirical boundary
between `significant_up`/`significant_down` and `not_significant` rows —
not just reading the KG's `table_scope`/`table_scope_detail` text, but
checking what padj and log2FC values the data itself shows at that
boundary. This directly answers the researcher's question: what p-value
was actually used, and does the source table report every tested gene or
only a filtered subset.

**Co-defined with the researcher** (three judgment calls):
1. Hypothesis 3 (noise/bootstrap check) runs on nitrogen only. All 5
   phosphorus experiments are pre-filtered (`significant_only` or
   `filtered_subset`) — none has an unbiased "every gene we tested"
   population, so a phosphorus background/bootstrap would secretly compare
   against an already-cherry-picked set. Phosphorus is still analyzed for
   hypotheses 1 and 2 (matched and cross response), just not the noise
   check.
2. MIT9312 stays in the phosphorus strain pool (2 genes' worth of real
   evidence; no scientific downside to keeping it).
3. Positive controls: `ntcA`, `glnA`, `amtB/amt1`, `ureA` (nitrogen);
   `pstS`, `phoA`, `phoB`, `phoR` (phosphorus) — all 8 confirmed present
   in the step-2 evidence set. Negative/background: the full population of
   genes actually tested in each of the 4 nitrogen `all_detected_genes`
   experiments (no phosphorus equivalent exists — see judgment call 1),
   excluding the 92 target genes and any gene whose product/name contains
   a nitrogen-metabolism keyword. No pre-set background size — the full
   filtered population is kept; the bootstrap (steps 4-5) draws
   size-matched random samples from it.

**Built the control sets** (`scripts/02_build_controls.py`): resolved the
8 positive-control genes to locus tags per strain; built the nitrogen
background pool by pulling all DE rows from the 4 unbiased nitrogen
experiments and excluding target-list genes and nitrogen-keyword matches
in gene name/product text (keywords: nitrogen, nitrate, nitrite,
ammonium, ammonia, urea, urease, cyanate, glutamine synthetase, glutamate
synthase, GOGAT, nitrogenase, amino acid transport).

## Results

**Significance criteria per experiment** (full table:
`data/01_significance_criteria.csv`):

| Nutrient | Strain | Omics | Test | Table contains | Effective threshold |
|---|---|---|---|---|---|
| N | MED4 | Proteomics | DESeq2 | all tested genes | padj<0.05 AND \|log2FC\|>~1.0 |
| N | MED4 | RNA-seq | DESeq2 | all tested genes | padj<0.05 AND \|log2FC\|>~1.0 |
| N | MED4 | RNA-seq | Rockhopper | top 50% by expression | padj<0.05, \|log2FC\|>~1.0 |
| N | MED4 | Microarray | Goldenspike | all tested genes | padj<0.01, no fold-change floor |
| N | MIT9313 | Microarray | Goldenspike | all tested genes | padj<0.01, no fold-change floor |
| P | NATL2A | RNA-seq | DESeq2 | genes sig. at >=1 timepoint | padj<0.05 AND \|log2FC\|>~1.0 |
| P | MIT9312 | Proteomics | iTRAQ t-test | **fold-change cutoff only, no p-value** | \|fold-change\|>1.6 (no padj reported) |
| P | NATL2A | Proteomics | iTRAQ t-test | **fold-change cutoff only, no p-value** | \|fold-change\|>1.6 (no padj reported) |
| P | MED4 | Microarray | Cyber-T | q<0.05 at t=48h, all timepoints shown | padj<0.05 (at the qualifying timepoint) |
| P | MIT9313 | Microarray | Cyber-T | q<0.05 at t=48h, all timepoints shown | padj<0.05 (at the qualifying timepoint) |

Two things worth flagging plainly: (1) the DESeq2- and Rockhopper-based
experiments use a **dual criterion** (padj AND a fold-change floor around
2-fold, i.e. \|log2FC\|>1) — a gene can have a very low padj and still be
called `not_significant` if its fold-change is small; this is visible in
the data (e.g. MED4 nitrogen proteomics has a `not_significant` row with
padj=3.6e-05, far below the 0.05 line, because its fold-change didn't
clear the floor). (2) The two Fuszard iTRAQ proteomics experiments
(MIT9312, NATL2A phosphorus) report **no p-value at all** for their
significance call — despite running a t-test, the KG's significance flag
for these rows is purely the paper's fold-change cutoff (>1.6-fold up or
<0.6-fold down). This means "significant" means something different across
platforms in this analysis, not just a different number — a caveat that
carries into how results get compared in step 5.

**Positive controls** (`data/02_positive_controls.csv`): 4 nitrogen genes
(`ntcA`, `glnA`, `amtB/amt1`, `ureA`) and 4 phosphorus genes (`pstS`,
`phoA`, `phoB`, `phoR`), 31 (gene x strain) rows total across the 4
in-scope strains.

**Negative/background pool** (`data/02_negative_background_pool.csv`,
nitrogen only): built from the 4 unbiased nitrogen experiments, no
pre-set size —

| Experiment | Distinct genes tested | Excluded (target list) | Excluded (N-keyword) | Background pool |
|---|---|---|---|---|
| Weissberg proteomics (MED4) | 1,424 | 96 rows | 66 rows | 1,387 genes |
| Weissberg RNA-seq (MED4) | 1,849 | 35 rows | 24 rows | 1,808 genes |
| Tolonen microarray (MED4) | 1,697 | 210 rows | 144 rows | 1,656 genes |
| Tolonen microarray (MIT9313) | 2,241 | 228 rows | 138 rows | 2,196 genes |

(Exclusion row counts are per-row, not per-gene — a gene tested at
multiple timepoints contributes multiple rows.)

## Surprises

- The dual significance criterion (padj + fold-change floor) on the
  DESeq2/Rockhopper experiments wasn't visible from `table_scope` alone —
  it only showed up by checking the actual padj values at the
  significant/not-significant boundary. Without this check we would have
  assumed a single padj<0.05 cutoff and misreported it.
- The two Fuszard iTRAQ experiments report no padj whatsoever — every
  `padj` cell is null. Their "significant" call is 100% fold-change-based.
  This wasn't visible from `table_scope_detail` text alone either
  (`"Only proteins reaching the paper's fold-change cutoff..."` reads like
  a filter description, not "there is no p-value in this table at all"
  until the data is actually checked).

## Decide-gate checklist

- **Outputs produced:** `scripts/01_significance_criteria.py`,
  `scripts/02_build_controls.py`; `data/01_significance_criteria.csv`,
  `data/02_positive_controls.csv`, `data/02_negative_background_pool.csv`
  (4.6MB / 29k rows, gitignored as a large KG-reproducible extract —
  regenerate with `uv run scripts/02_build_controls.py`).
- **Results presented:** significance-criteria table; positive-control
  list; negative-background-pool table (all above).
- **QC gate:** significance boundary checked empirically (not assumed
  from `table_scope` text) for all 10 experiments; all 8 positive-control
  genes confirmed resolved in >=1 in-scope strain; background pool
  exclusion counts checked (target-list and N-keyword rows both
  nonzero and plausible relative to pool size).
- **Decisions made this step:** H3 restricted to nitrogen only; MIT9312
  kept; positive controls set; negative-background construction rule set
  with no pre-set size (all dated 2026-08-10, co-defined with researcher).
- **Advance rationale:** hypotheses are now operational (what "responds"
  means per experiment, what counts as a valid noise-check comparison),
  controls are materialized as concrete gene/locus-tag lists, and the
  significance-criteria table answers the researcher's standing question
  about what threshold each platform actually used. Ready for step 4 to
  build the method that computes hit rates and the bootstrap comparison.
