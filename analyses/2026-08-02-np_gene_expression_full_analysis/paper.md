# N/P acquisition gene expression response — full analysis

## Question

Do the 54 nitrogen/phosphorus-acquisition genes with confirmed in vivo
expression evidence (identified in the prior triage analysis,
`analyses/2026-07-13-n_p_genes_in_vivo_experiments/`) show a significant
expression change under nutrient starvation in axenic Prochlorococcus —
both matched (nitrogen-annotated genes under nitrogen starvation,
phosphorus-annotated genes under phosphorus starvation) and cross
(nitrogen-annotated genes under phosphorus starvation, phosphorus-annotated
genes under nitrogen starvation)?

## Background

The 54 genes analyzed here (29 nitrogen-annotated, 25 phosphorus-annotated)
were identified as having in vivo expression evidence in the prior triage
analysis (`analyses/2026-07-13-n_p_genes_in_vivo_experiments/`), which
screened all 91 genes in the researcher's N/P-acquisition gene spreadsheet
against the KG for the mere presence of expression data.

The KG holds 10 axenic Prochlorococcus nitrogen/phosphorus starvation
experiments used here, across 6 publications: 5 nitrogen (strict
starvation/deprivation contrasts — 4 experiments testing growth on an
alternate N source such as cyanate/urea/nitrite were excluded as a
different biological question) and 5 phosphorus (limitation/deplete-vs-
replete contrasts; metabolomics-only experiments excluded, as they carry
no gene-level differential expression), spanning 4 Prochlorococcus
strains in the researcher's study: MED4, MIT9312, MIT9313, and NATL2A
(SS120/CCMP1375 excluded — not part of this study; MIT9301 excluded —
its only phosphorus data was metabolomics-only). 53 of the 54 target
genes resolved to at least one locus tag among these 4 strains
(`ptxB/phnD2`'s only evidence was in the now-excluded MIT9301). See
`2_kg_selection/notebook.md` for the full experiment table and
resolution detail.

## Methods

Three hypotheses are tested: (H1) nitrogen-annotated genes are
upregulated under nitrogen starvation and phosphorus-annotated genes
under phosphorus starvation (matched, directional); (H2) genes are also
differentially expressed under the non-matching nutrient's starvation —
nitrogen-annotated genes under phosphorus starvation and vice versa
(cross, non-directional); (H3) a noise/method-limitation check comparing
the target genes' significant-hit rate against a same-experiment
background set. Positive controls are 8 canonical nutrient-starvation
marker genes drawn from the target set itself (`pstS`, `phoA`, `phoB`,
`phoR` for phosphorus; `ntcA`, `glnA`, `amtB/amt1`, `ureA` for nitrogen).
The negative control is 303 genes (up to 100 per strain) sampled from
each strain's actually-tested gene pool in these experiments, excluding
the target genes and any gene with N/P-metabolism-related product/name
text. See `3_analysis_framing/notebook.md` for full framing detail and
`2_kg_selection/notebook.md` for the underlying experiment/gene
selection.

All three hypotheses are computed with one shared function,
`4_methods/hit_rate.py`: for a gene set and an experiment set, the
fraction of (gene x experiment x timepoint) tests with a significant
DE call, split by direction. Verified against a hand tally on `ntcA`/
`glnA` (MED4, nitrogen experiments) before use — see
`4_methods/notebook.md`.

## Results

**Nitrogen side (unbiased gene coverage, valid H1/H3 comparison).**
Nitrogen-annotated target genes are significantly differentially
expressed in 28.8% of (gene x experiment x timepoint) tests under
nitrogen starvation (27.5% up, 1.2% down), versus 14.8% for a
same-experiment background gene set (4.9% up, 9.9% down) — roughly
double the overall rate, with a qualitatively different direction
profile (target genes almost entirely up; background genes lean down).
The 4 canonical nitrogen positive-control genes (`ntcA`, `glnA`,
`amtB/amt1`, `ureA`) score higher still: 57.9% significant, 100% of
that upregulated.

**Phosphorus side (numbers reported, no noise check).** All 5
phosphorus experiments in the KG are pre-filtered by the source
publications to already-significant genes (`table_scope`
`significant_only`/`filtered_subset`, with the cutoff recorded per
experiment in `2_kg_selection/data/02_np_experiments.csv`; nitrogen is
4/5 unfiltered by contrast) — this inflates every phosphorus-side rate.
Phosphorus-annotated target genes are significant in 40.7% of tests
(39.3% up), phosphorus positive controls in 49.2% (47.6% up). No
phosphorus background/H3 rate is reported: since the background pool
would be drawn from the same pre-filtered tables as the target genes,
it would not be a random sample of the genome and cannot serve as a
noise check.

**Cross-nutrient (H2).** Phosphorus-annotated genes under nitrogen
starvation: 15.7% significant across 395 tests (6.8% up, 8.9% down) —
between the nitrogen background rate (14.8%) and the nitrogen
target-gene rate (28.8%), closer to background. Nitrogen-annotated
genes under phosphorus starvation: only 6 total tests (2 of 93
gene-instances had any data at all) — too small to interpret, a direct
consequence of the phosphorus tables' narrow, pre-filtered coverage.

Full tables: `5_analyze/data/01_hit_rate_summary.csv`,
`01_hit_rate_per_gene.csv`. Figure: `5_analyze/figures/01_hit_rate_comparison.png`.

## Discussion

## References
