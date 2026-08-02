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
`4_methods/notebook.md`. The one significance test in this analysis
(Fisher's exact, nitrogen target genes vs. nitrogen background) uses
`scipy.stats.fisher_exact` (scipy 1.17.1). KG version: multiomics-kg
0.1.0-alpha.6 (production); explorer: multiomics_explorer 0.1.0a4;
template: 0.1.0-alpha.2.

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
starvation: 15.7% significant across 395 tests (6.8% up, 8.9% down),
tested directly against nitrogen background (14.8%) — fold enrichment
1.06x, p = 0.64 (Fisher's exact), statistically indistinguishable from
background. Nitrogen-annotated genes under phosphorus starvation: only
6 total tests (2 of 93 gene-instances had any data at all) — too small
to interpret, a direct consequence of the phosphorus tables' narrow,
pre-filtered coverage.

Full tables: `5_analyze/data/01_hit_rate_summary.csv`,
`01_hit_rate_per_gene.csv`. Figure: `5_analyze/figures/01_hit_rate_comparison.png`.

**Fold-change magnitude (not just the significance call).** All of the
above uses the categorical significant/not-significant/direction call;
the underlying log2 fold-change values show the same pattern more
sharply. Among significant hits, nitrogen target genes and positive
controls cluster tightly positive (median log2FC 2.3-2.4, almost no
negative points); the nitrogen background set's significant hits are
predominantly *negative* (median magnitude 1.7, downward-skewed) — the
mirror image of the target genes, not just a lower rate. This
direction pattern holds within each omics platform separately
(RNA-seq, proteomics, microarray), not only in aggregate — platform
magnitudes are not pooled since they are not directly comparable.
Phosphorus target genes and positive controls also skew positive
(median log2FC 2.6-3.1) with a wider spread, including 3 known outlier
genes (`PMM0707`, `PMM0708`=`phoA`, `PMM1416`; MED4 phosphate-
starvation microarray, log2FC 30-162) independently flagged as likely
artifacts by the prior triage analysis. Figures:
`5_analyze/figures/02_log2fc_distribution.png`,
`03_log2fc_by_platform.png`.

## Discussion

A Fisher's exact test on the one valid comparison in this analysis —
nitrogen-annotated target genes (116/403 significant) vs. a
same-experiment background gene set (236/1596 significant) — found a
1.95-fold enrichment (odds ratio 2.33, p = 3.2e-10). More informative
than the p-value itself is the direction split: target genes are
upregulated in 96% of their significant hits, while the background
set's significant hits lean the other way (33% up / 67% down). This
combination — higher hit rate *and* a specific, predicted direction —
is evidence against H3 (noise/method-limitation) and for H1
(nitrogen-annotated genes are upregulated under nitrogen starvation)
on the nitrogen side. H1 also holds for the 4 nitrogen positive
controls (`ntcA`, `glnA`, `amtB/amt1`, `ureA`), which score higher
still (57.9% significant, 100% up).

The phosphorus side shows the same qualitative pattern (39.3% up vs.
1.4% down for target genes; 47.6% vs. 1.6% for positive controls) but
cannot be statistically tested: every phosphorus experiment in this KG
build is pre-filtered by the source publication to already-significant
genes, so there is no unbiased phosphorus gene population to serve as
background. `[interpretation]`: the sharp up:down asymmetry survives
even though the underlying publications' inclusion filters are
symmetric in magnitude (e.g. fold-change >1.6 or <0.6, either
direction), which is suggestive that the phosphorus pattern reflects
real biology rather than being purely a curation artifact — but this
is not proof, and a properly powered test would need phosphorus data
this KG build does not have.

H2 (cross-nutrient) is not supported by the current data — and for the
nitrogen-experiment side, this is a tested negative result, not an
absence of evidence: phosphorus genes under nitrogen starvation (15.7%
significant) are statistically indistinguishable from nitrogen
background (14.8%; fold enrichment 1.06x, p = 0.64, Fisher's exact).
Nitrogen genes under phosphorus starvation cannot be assessed at all —
only 6 tests exist across the entire gene set, a direct consequence of
the phosphorus tables' narrow, pre-filtered coverage.

**Overall:** the nitrogen-side result is a clean, statistically
supported confirmation that this gene list's nitrogen-annotated
members respond to nitrogen starvation specifically, in the predicted
direction, more than a random gene would. The phosphorus-side result
is directionally consistent with the same claim but not independently
confirmable with this KG build's data. The cross-nutrient hypothesis
(H2) is a tested negative on the nitrogen-experiment side (phosphorus
genes do not respond to nitrogen starvation beyond background) and
untestable on the phosphorus-experiment side.

See `6_evaluate/notebook.md` for the full caveat list (table_scope
pre-filtering, excluded strains/experiments, uneven strain coverage,
pseudo-replication in the test's row count, and others).

## References

1. Weissberg O, Aharonovich D, Sher D (2025). Transcriptomic and
   Proteomic Analysis Reveals Nitrogen Recycling as a Core Mechanism
   for Prochlorococcus Prolonged Survival. *bioRxiv*.
   https://doi.org/10.1101/2025.11.24.690089
2. Read RW, Berube PM, Biller SJ, Neveux I, Cubillos-Ruiz A, Chisholm
   SW, Grzymski JJ (2017). Nitrogen cost minimization is promoted by
   structural changes in the transcriptome of N-deprived
   Prochlorococcus cells. *The ISME Journal*.
   https://doi.org/10.1038/ismej.2017.88
3. Lin X, Ding H, Zeng Q (2015). Transcriptomic response during phage
   infection of a marine cyanobacterium under phosphorus-limited
   conditions. *Environmental Microbiology*.
   https://doi.org/10.1111/1462-2920.13104
4. Fuszard MA, Wright PC, Biggs CA (2012). Comparative quantitative
   proteomics of Prochlorococcus ecotypes to a decrease in
   environmental phosphate concentrations. *Aquatic Biosystems*.
   https://doi.org/10.1186/2046-9063-8-7
5. Tolonen AC, Aach J, Lindell D, Johnson ZI, Rector T, Steen R,
   Church GM, Chisholm SW (2006). Global gene expression of
   Prochlorococcus ecotypes in response to changes in nitrogen
   availability. *Molecular Systems Biology*.
   https://doi.org/10.1038/msb4100087
6. Martiny AC, Coleman ML, Chisholm SW (2006). Phosphate acquisition
   genes in Prochlorococcus ecotypes: Evidence for genome-wide
   adaptation. *PNAS*. https://doi.org/10.1073/pnas.0601301103
