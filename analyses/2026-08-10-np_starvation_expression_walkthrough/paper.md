# N/P acquisition gene expression response under starvation — walkthrough analysis

## Question

For the 92 nitrogen- and phosphorus-acquisition genes in the researcher's
reference list, do genes annotated for a given nutrient show a
differential-expression response when axenic, uninfected Prochlorococcus —
restricted to the researcher's 15 study strains (MED4, MIT9312, MIT9313,
MIT1327, MIT0604, NATL2A, MIT9515, MIT9215, AS9601, PAC1, MIT9202, SB,
MIT9301, MIT1314, NATL1A) — is grown under literal starvation/limitation of
that nutrient (N or P withheld or limited in the medium) versus
nutrient-replete medium — and is that response nutrient-specific (matched:
N-genes respond to N starvation, P-genes to P starvation) or general
(cross: genes respond regardless of which nutrient is withheld)? Evidence
sources: RNA-seq, proteomics, and microarray differential expression in the
KG. A third question asks whether any observed matched/cross pattern is
statistically distinguishable from a background/noise rate, using
significance testing and bootstrap resampling.

## Background

The KG holds 12 Prochlorococcus nitrogen-treatment and 9 phosphorus-treatment
experiments. Applying the step-1 scope (axenic, uninfected, literal
starvation/limitation, restricted to the researcher's 15 study strains)
narrows this to 10 usable experiments across 4 strains (MED4, MIT9312,
MIT9313, NATL2A) and 6 publications: 5 nitrogen (Weissberg et al. 2025 —
proteomics and RNA-seq, axenic arm only; Read et al. 2017 — RNA-seq;
Tolonen et al. 2006 — microarray, MED4 and MIT9313) and 5 phosphorus (Lin,
Ding, Zeng 2015 — RNA-seq, uninfected arm only; Fuszard et al. 2012 —
proteomics, MIT9312 and NATL2A; Martiny, Coleman, Chisholm 2006 —
microarray, MED4 and MIT9313). Excluded: coculture- and
phage-infection-background arms, 2 strains not in the researcher's 15
(SS120), 4 alternate-N-source substrate-switch experiments, and 2
metabolomics-only experiments with no gene-level differential expression.

Of the researcher's 92-gene reference list, 91 carry a Cyanorak
ortholog-group ID; resolving each via `genes_by_homolog_group` against the
4 in-scope strains finds a locus tag for 61 genes in at least one strain
(MED4 53, NATL2A 48, MIT9312 44, MIT9313 42). A locus tag only means the
gene is present in that strain's genome, not that it was measured;
checking each resolved gene against the 10 in-scope experiments'
differential-expression rows narrows this to **60 genes (29
nitrogen-annotated, 31 phosphorus-annotated) with actual expression
evidence** in at least one in-scope experiment. The remaining genes have
no locus tag in any of the 4 strains — consistent with known
strain/ecotype-level gene content variation in Prochlorococcus rather than
a resolution gap (one exception, `phnW`, has no Cyanorak ID at all and
resolves only in *Pseudomonas putida* on direct name lookup — genuinely
absent from Prochlorococcus in this KG build). One strain, MIT9312,
resolves 44 genes but has evidence for only 2, because its sole in-scope
experiment's source table only reports genes that were already
significant in the original publication (see `gaps_and_friction.md`). See
`2_kg_selection/notebook.md` for the full experiment and gene resolution
tables.

## Methods

Three hypotheses: (H1, matched) nitrogen-annotated genes respond to
nitrogen starvation and phosphorus-annotated genes respond to phosphorus
starvation; (H2, cross) genes also respond to the non-matching nutrient's
starvation; (H3, noise check) any observed pattern is tested against a
background/bootstrap comparison — **nitrogen only**, since all 5
phosphorus experiments are pre-filtered by their source publication
(`table_scope` `significant_only`/`filtered_subset`) and have no unbiased
gene population to serve as background.

Each experiment's significance call was checked empirically (not assumed
from metadata): the DESeq2- and Rockhopper-based experiments (Weissberg
proteomics/RNA-seq, Read RNA-seq, Lin RNA-seq) use a dual criterion,
padj<0.05 **and** \|log2FC\|>~1; the Goldenspike microarray experiments
(Tolonen) use padj<0.01 alone; the Cyber-T microarray experiments
(Martiny) use padj<0.05 at the 48h timepoint; the two Fuszard iTRAQ
proteomics experiments (MIT9312, NATL2A phosphorus) report **no p-value
at all** — their significance flag is a pure fold-change cutoff
(>1.6-fold). Full table: `3_analysis_framing/data/01_significance_criteria.csv`.

Positive controls: `ntcA`, `glnA`, `amtB/amt1`, `ureA` (nitrogen);
`pstS`, `phoA`, `phoB`, `phoR` (phosphorus). Negative/background (nitrogen
only): the full population of genes tested in the 4 nitrogen
`all_detected_genes` experiments (1,387-2,196 genes each), excluding the
92 target genes and nitrogen-keyword product/name matches — no pre-set
size; the bootstrap draws random samples matching the nitrogen
target-gene-set size from this population. See
`3_analysis_framing/notebook.md` for full detail.

Two functions (`4_methods/np_response.py`) compute all three hypotheses:
`hit_rate()` (fraction of gene x experiment x timepoint tests
significant, split by direction — H1 and H2) and `bootstrap_pvalue()`
(nitrogen-only empirical significance test — H3), resampling
size-matched gene sets from the background pool 10,000 times per
comparison. Both verified against hand-computed toy data before use, and
`hit_rate()` additionally verified against an independent manual tally of
`ntcA`/`glnA`'s real DE rows in MED4's 4 nitrogen experiments (26 rows;
61.5% and 84.6% significant respectively, both consistent with their
known role as nitrogen-starvation markers). See `4_methods/notebook.md`.

## Results

Each of the 10 in-scope experiments was restricted to a single
representative starvation timepoint (rather than pooling across a time
course), chosen from the KG's `growth_phase` field where it distinguishes
a starved state, and from the source publication directly for 2 rapid-
onset nitrogen experiments where the KG's own annotation never reaches
"nutrient_limited" (Tolonen et al. 2006's MED4 and MIT9313 microarray
time courses mark every timepoint `acute_stress`; the source paper
reports cultures began declining at 12h). Full table and rationale:
`5_analyze/notebook.md`.

**Nitrogen (matched, H1):** 29 nitrogen-annotated target genes are
significant in 41.3% of (gene x experiment) tests (43/104), **entirely
upregulated (0 downregulated)**. The 4 canonical positive controls
(`ntcA`, `glnA`, `amtB/amt1`, `ureA`) score higher still: 80.0%
significant, 100% of that upregulated. Tested against a same-timepoint
background pool via 10,000-iteration bootstrap (nitrogen only — see
Methods): the background's null hit-rate distribution has mean 18.8%
(std 3.8%); the observed 41.3% rate has an empirical p-value < 0.0001 (0
of 10,000 random draws reached it). Fisher's exact test on the same
comparison agrees: odds ratio 3.01, p = 2.15e-07.

**Phosphorus (matched, H1, no noise check per step 3):** 18
phosphorus-annotated target genes are significant in 75.7% of tests
(26/37), mostly upregulated (70.3% up, 5.4% down). The 4 positive
controls (`pstS`, `phoA`, `phoB`, `phoR`) score 84.6% significant (76.9%
up, 7.7% down).

**Cross-nutrient (H2):** phosphorus-annotated genes tested under nitrogen
starvation are significant in 20.2% of tests (9 up, 13 down out of 109) —
no clean directional pattern like the matched-nitrogen result (roughly
balanced up/down). Nitrogen-annotated genes tested under phosphorus
starvation have only 1 test in the entire dataset — too sparse to
interpret, a direct consequence of the phosphorus tables' narrow,
pre-filtered gene coverage (step 3).

**Main figure:** `5_analyze/figures/01_gene_experiment_heatmap.png` — all
10 experiments (rows, shown separately per study, not averaged — a
concordance check on 9 genes measured by two independent MED4-nitrogen
studies found direction always agreed but 2 of 9 disagreed on the
significance call at the margin, which an average would have hidden) x
the 40 of 61 genes with a response in >=1 experiment (columns), colored
by upregulated / downregulated / tested-not-significant / no-data.

**Upregulation-rate figure:** `5_analyze/figures/02_pct_upregulated.png`
— for every gene, the percentage of its tested experiments (any nutrient)
that came back significantly upregulated. 6 nitrogen genes reach 100%
(`ntcA`, `cynA`, `cynD`, `focA`, `nirA`, `nirX`; the latter 3 from a
single test each). No phosphorus gene reaches 100%; the highest is
`PMM719` at 80% (4/5).

**Target-vs-background strip plot:** `5_analyze/figures/03_pct_scatter.png`
— the same per-gene percentages plotted as individual points against the
nitrogen background pool (n=4,019 background genes): mean 40.5% for
nitrogen target genes vs. 6.1% for the nitrogen background, a visual
complement to the H3 bootstrap p-value showing the whole distribution
shifted, not just a summary statistic. Phosphorus target genes shown
alongside (mean 19.4%) with no background comparison, since none is
statistically valid (step 3).

**Nutrient-collapsed heatmap:**
`5_analyze/figures/04_pct_heatmap_by_nutrient.png` — two stacked panels,
one per gene group (N-acquisition genes on top, P-acquisition genes
below), each with its own pair of rows ("N starvation experiments", "P
starvation experiments"). Cell = % of that gene's *tested* experiments
(no-data cells excluded from the denominator) that came back significant
in whichever direction had more hits, color intensity scaled continuously
by the percentage; the number is only printed inside the cell at >=55%
(a readability threshold, not a data cutoff — cells below it are still
colored by their exact percentage). Makes the H1-vs-H2 (matched-vs-cross)
pattern visible per gene: N-acquisition genes are red/upregulated in the
top panel's matched row and mostly gray/hatched in its cross row, and the
mirror holds for P-acquisition genes in the bottom panel.

## Discussion

## References

1. Weissberg O, Aharonovich D, Sher D (2025). Transcriptomic and Proteomic
   Analysis Reveals Nitrogen Recycling as a Core Mechanism for
   Prochlorococcus Prolonged Survival. *bioRxiv*.
   https://doi.org/10.1101/2025.11.24.690089
2. Read RW, Berube PM, Biller SJ, Neveux I, Cubillos-Ruiz A, Chisholm SW,
   Grzymski JJ (2017). Nitrogen cost minimization is promoted by
   structural changes in the transcriptome of N-deprived Prochlorococcus
   cells. *The ISME Journal*. https://doi.org/10.1038/ismej.2017.88
3. Tolonen AC, Aach J, Lindell D, Johnson ZI, Rector T, Steen R, Church
   GM, Chisholm SW (2006). Global gene expression of Prochlorococcus
   ecotypes in response to changes in nitrogen availability. *Molecular
   Systems Biology*. https://doi.org/10.1038/msb4100087
4. Lin X, Ding H, Zeng Q (2015). Transcriptomic response during phage
   infection of a marine cyanobacterium under phosphorus-limited
   conditions. *Environmental Microbiology*.
   https://doi.org/10.1111/1462-2920.13104
5. Fuszard MA, Wright PC, Biggs CA (2012). Comparative quantitative
   proteomics of Prochlorococcus ecotypes to a decrease in environmental
   phosphate concentrations. *Aquatic Biosystems*.
   https://doi.org/10.1186/2046-9063-8-7
6. Martiny AC, Coleman ML, Chisholm SW (2006). Phosphate acquisition
   genes in Prochlorococcus ecotypes: Evidence for genome-wide adaptation.
   *PNAS*. https://doi.org/10.1073/pnas.0601301103
