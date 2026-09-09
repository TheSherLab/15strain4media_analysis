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
resolves 44 genes but has evidence for only 3, because its sole in-scope
experiment's source table only reports genes that were already
significant in the original publication (see `gaps_and_friction.md`).

**Gene-identity corrections (2026-09-08).** The Cyanorak ID is not
reliably one gene: two independent step-2 checks — a multi-member audit
(`05_paralog_audit.py`: count `genes_by_homolog_group` members per strain;
>1 = a family, not a gene) and a name-vs-ID cross-check
(`06_name_vs_cyanorak_crosscheck.py`: resolve each gene by name too, diff)
— found 6 mis-resolved loci. `phoE`'s ID (`CK_00002330`) is a 3–6-member
"outer membrane porin" family per strain, and the resolver had been
silently keeping a non-pho-island paralog in MED4, MIT9312 and NATL2A;
`unkP2` was similar. `urtA` had `urtE`'s Cyanorak ID in the source
spreadsheet (so the analysis carried `urtE` twice, no `urtA`); `ptrA`'s
ID was a MED4-only singleton group, missing the NATL2A ortholog. Each
correction is confirmed by product annotation, ortholog-group membership,
genomic synteny to the pho/N operon, and expression behaviour (the
corrected `phoE` locus jumps ~100-fold at Martiny 48h). The 60-gene /
29-N / 31-P totals are unchanged — the fixes moved which locus a gene
points at, not the gene count. See `2_kg_selection/notebook.md`
("Reopened 2026-09-08") for the full method and the 6-correction table.

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

For each experiment one representative starvation timepoint is used (not a
pooled time course). Two choices required judgment beyond the KG's fields:
the Tolonen nitrogen microarrays' 12h point (KG marks all timepoints
`acute_stress`; source paper reports decline at 12h), and the Lin
phosphorus experiment's **59h** point — the latest `nutrient_limited`
timepoint in the arm that never had phosphate re-added; the same-hour
`59h (P added)` and `50h (P added)` recovery-arm timepoints are excluded.
(Lin was corrected from 46h to 59h on 2026-09-08: at 46h only 4 of Lin's
34 reported genes are significant, at 59h 18 are — the pho/pst regulon is
not fully engaged until 59h. See `gaps_and_friction.md`.)

**Table-absent reclassification.** Where a genome-wide assay is reported
only as a significance-filtered gene list, a target gene present in the
strain's genome but with no row in the table was measured and did not pass
— counted `not_significant`, not "no data". Applied to the 2 Martiny
phosphorus microarrays (2026-09-08; `table_scope = filtered_subset`,
q<0.05 at 48h) and, from 2026-09-09, to the Lin uninfected phosphorus
RNA-seq (`table_scope = significant_only`; RNA-seq measures the whole
transcriptome — same logic as Martiny, correcting an earlier over-cautious
"not applied to Lin" call). The Lin extension moved 35 cells from "no
data" to "not significant" (10 P-acquisition genes and 25 N-acquisition
genes have a NATL2A locus but sit outside Lin's 34-gene operon table); no
significant call changed, so it widened denominators only — H1 P-matched
60.6% → 52.6%, H2 N-in-P 4.2% → 2.7%, everything else unchanged.
Reclassification is *not* applied to Fuszard (iTRAQ holds a fixed
detected-protein list — absence there
means the peptide was not detected, verified against the live KG).

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

Full timepoint table and rationale: `5_analyze/notebook.md`.

All figures below are 2026-09-08 (gene-identity corrections, Lin at 59h,
Martiny reclassification) with a 2026-09-09 update that added the Lin
uninfected RNA-seq to the table-absent reclassification (see Methods) —
this widened the phosphorus tested denominators (no significant call
changed). Pre-2026-09-08 values are given in parentheses where a number
moved; 2026-09-08 → 2026-09-09 moves are called out inline.

**Nitrogen (matched, H1):** 29 nitrogen-annotated target genes are
significant in **44.2%** of (gene x experiment) tests (46/104), **entirely
upregulated (0 downregulated)** (was 41.3%, 43/104 — the real `urtA`
locus is more N-responsive than the `urtE` locus step 2 had been using).
The 4 canonical positive controls (`ntcA`, `glnA`, `amtB/amt1`, `ureA`)
score 80.0% significant, 100% of that upregulated. Tested against a
same-timepoint background pool via 10,000-iteration bootstrap (nitrogen
only — see Methods): the background's null hit-rate distribution has mean
18.8% (std 3.8%); the observed 44.2% rate has an empirical p-value <
0.0001 (0 of 10,000 random draws reached it). Fisher's exact test agrees:
**odds ratio 3.39, p = 4.98e-09** (was OR 3.01, p = 2.15e-07). The
nitrogen background pool is unchanged by every 2026-09-08 fix.

**Phosphorus (matched, H1, no noise check per step 3):** 31
phosphorus-annotated target genes are significant in **52.6%** of tests
(40/76), mostly upregulated (50.0% up, 2.6% down). The 40 significant
calls are the same as at 2026-09-08 (60.6%, 40/66); the 2026-09-09 Lin
reclassification added 10 tested-not-significant tests to the denominator.
Pre-2026-09-08 the rate was 75.7% (26/37 — the Martiny reclassification
removed an inflated all-significant denominator; the `phoE` fix added a
canonical responder). The 4 positive controls (`pstS`, `phoA`, `phoB`,
`phoR`) score **92.9%** significant (13/14; was 84.6%, 11/13 — `phoA` and
`phoR` are significant at Lin 59h but not 46h) — the controls are all in
Lin's table so the 2026-09-09 change did not touch them.

**Cross-nutrient (H2):** phosphorus-annotated genes tested under nitrogen
starvation are significant in 20.7% of tests (10 up, 13 down out of 111)
— no clean directional pattern like the matched result (roughly balanced
up/down). Nitrogen-annotated genes tested under phosphorus starvation are
significant in **2.7% of 73 tests (2 of 73, both down)** — was 4.2% of 48
(same 2 significant; the 2026-09-09 Lin reclassification added 25 flat
cross-tests), and 1 test total before the 2026-09-08 Martiny
reclassification. The genome-wide phosphorus assays' nitrogen-gene
cross-tests are almost all flat.

**Main figure:** `5_analyze/figures/01_gene_experiment_heatmap.png` (and
`.pptx` — an editable version, every cell a named rectangle, from
`scripts/07_heatmap_pptx.py`) — all 10 experiments (rows, shown
separately per study, not averaged; each labelled by its Oxford-style
citation `First author et al., YEAR - analysis type - strain`, grouped by
a left-side `Nitrogen starvation` / `Phosphorus starvation` bracket) x the
~40 of 61 genes with a response in >=1 experiment (columns), colored by
upregulated / downregulated / tested-not-significant / no-data (gene
present, not in that experiment's table) / gene-absent-from-strain (no
locus tag, drawn with three diagonal strokes on a near-white cell). Arial,
no bold; every cell has a light grey border so the three low-signal states
are separable. `[KG]`
The absent-gene cells fall on a consistent set of gene x strain combinations:
the cyanate operon (`cynA/B/D/S`) and several hypothetical P-region genes
(`PMM707/719/721`, `psiP1`, `phoA`, `ptrA`, `unkP5`) have no locus tag in
MIT9313; `nirA`/`nirX` and `focA` have none in MED4.
`[interpretation, researcher-provided]` These match known Prochlorococcus
strain/ecotype gene-content differences (e.g. `focA` as a
low-light-ecotype gene), not a resolution gap.

**Upregulation-rate figure:** `5_analyze/figures/02_pct_upregulated.png`
— for every gene, the percentage of its **matched-nutrient** tested
experiments (cross-nutrient tests excluded — 2026-09-08) that came back
significantly upregulated. 7 nitrogen genes reach 100% under N starvation
(`ntcA`, `glnA`, `cynA`, `cynD`, `focA`, `nirA`, `nirX` — the last 3 from
a single test each). Several phosphorus genes reach 100% under P
starvation, including `phoA`, `phoB`, `pstA`, `pstC`, `gap3` (all 3/3) and
`phoE` (4/4).

**Target-vs-background strip plot:** `5_analyze/figures/03_pct_scatter.png`
— the same per-gene matched-nutrient percentages plotted as individual
points against the nitrogen background pool (n=4,019 background genes):
**mean 43.1%** for nitrogen target genes vs. 6.1% for the nitrogen
background, a visual complement to the H3 bootstrap showing the whole
distribution shifted (the 43.1% per-gene mean tracks the 44.2% pooled hit
rate). Phosphorus target genes shown alongside (mean 50.4%) with no
background comparison, since none is statistically valid (step 3).
Restricting to matched-nutrient tests (2026-09-08) keeps this figure
aligned with the matched-nitrogen H3 comparison it illustrates.

**Nutrient-collapsed heatmap:**
`5_analyze/figures/04_pct_heatmap_by_nutrient.png` — two stacked panels,
one per gene group (N-acquisition genes on top, P-acquisition genes
below), each with its own pair of rows ("N starvation experiments", "P
starvation experiments"). Cell = % of that gene's *tested* experiments in
each nutrient that came back significant in whichever direction had more
hits, color intensity scaled by the percentage. This is where
cross-nutrient behaviour stays visible: N-acquisition genes are
upregulated in the top panel's matched row and almost entirely
gray/not-significant in its cross row, and the mirror holds for
P-acquisition genes below.

## Discussion

**Matched response (H1): supported for both nutrients.** Nitrogen-
annotated genes respond to nitrogen starvation in 44.2% of tests
(46/104), entirely upregulated (0 downregulated). Phosphorus-annotated
genes respond to phosphorus starvation in 52.6% of tests (40/76), mostly
upregulated (50.0% up, 2.6% down). [interpretation] The phosphorus rate
has stepped down twice, both times because a genome-wide phosphorus assay
reported only as a significance-filtered list had its table-absent genes
recounted as tested-not-significant rather than "no data": 75.7%
(pre-2026-09-08) → 60.6% (Martiny microarrays, 2026-09-08) → 52.6% (Lin
RNA-seq, 2026-09-09). The **same 40 genes** are significant at each step;
only the denominator grows. It is still above the nitrogen rate and still
not directly comparable: the 2 Fuszard tables remain detection-limited,
which is why the noise/bootstrap check (H3) stays nitrogen-only. Correcting `phoE`'s locus (it had been pointing
at a non-pho-island porin paralog) added a canonical pho-regulon
responder that is now significant in every phosphorus experiment testing
it.

Positive controls sharpen this per gene. Of the 4 canonical nitrogen
markers, `ntcA` and `glnA` are significant in 100% of their 5 tests each
(always upregulated), `ureA` in 80% (4/5), and `amtB/amt1` — the weakest
— in 40% (2/5); even `amtB/amt1` is never downregulated when significant.
Of the 4 phosphorus markers, `phoA` and `phoB` are significant in 100% of
3 tests each and `pstS` in 100% of 5 (4 up, 1 down — the only
downregulated hit among the 8 controls); `phoR` is at 66.7% (2/3).
`phoA` and `phoR` moved up from the pre-2026-09-08 figures because Lin now
sits at 59h, where the pho regulon is fully engaged.

**Noise/bootstrap check (H3, nitrogen only): supported, and not an
artifact of the 4 positive-control genes.** The observed 44.2% matched-
nitrogen rate exceeds every one of 10,000 random size-matched draws from
the background pool (null distribution: mean 18.8%, std 3.8%, maximum
32.6%) — empirical p<0.0001. Fisher's exact test agrees (odds ratio 3.39,
p=4.98e-09). A stability check reran both tests with the 4 positive
controls excluded: the remaining 25 genes (84 tests) still show 35.7%
significant against the same background's null mean of 18.7% — bootstrap
p<0.0001, Fisher's odds ratio 2.37, p=3.77e-04. The signal weakens
without the well-established markers but does not disappear — the matched-
nitrogen response is a property of the target gene set as a whole, not 4
famous genes carrying the rest. (All H3 numbers are marginally firmer
than the pre-2026-09-08 values — OR 3.01 → 3.39 — because the corrected
`urtA` locus is more N-responsive; the nitrogen background pool itself is
untouched by every 2026-09-08 fix.)

**Cross-nutrient response (H2): not supported as a general stress
response, and now demonstrated in both directions.** Phosphorus-annotated
genes tested under nitrogen starvation are significant in 20.7% of tests
(10 up, 13 down of 111) — well below the matched-nitrogen rate and, unlike
the matched result, split roughly evenly up/down. Nitrogen-annotated genes
tested under phosphorus starvation are significant in **2.7% of 73 tests**
(2 down of 73) — before the 2026-09-08 Martiny reclassification this
direction had only 1 test and could not be interpreted; Martiny then Lin
(2026-09-09) made the genome-wide phosphorus assays' nitrogen-gene
cross-tests visible, and they are almost entirely flat. [interpretation] A
low rate with no consistent direction in *both* cross directions is what
a nutrient-specific response predicts and a general-stress response would
not.

**What the negative background actually looks like.** The nitrogen
background pool's "noise" rate is not one number — it ranges from 3.0%
(Tolonen MIT9313 microarray) to 47.2% (Weissberg RNA-seq) across the 4
experiments it draws from, a roughly 16-fold spread. The Weissberg
RNA-seq background rate alone is higher than the observed matched-
nitrogen target rate (44.2%). [interpretation] This likely reflects a
broad transcriptional response to nitrogen starvation affecting many
genes genome-wide in that dataset, or a comparatively liberal DESeq2 call
specific to it, rather than "noise" in the everyday sense. The bootstrap
is not misled by this — it draws each iteration's sample from the same
per-experiment pools the real target genes were tested in, not from one
pooled rate — but the single pooled mean (18.8%) understates how much the
baseline "significant" rate varies by platform and experiment.

**Gene identity is only as good as the supplied Cyanorak IDs.**
[interpretation] The single biggest correction in this analysis was not
statistical — it was finding that 6 of the 61 resolved gene loci pointed
at the wrong gene, because a Cyanorak ortholog-group ID is not reliably
one gene (`phoE`'s ID is a porin family; `urtA`'s ID in the source list
was actually `urtE`'s). Two mechanical checks now catch this — a
member-count audit and a resolve-by-name cross-check — and both should be
standard for any analysis built on a supplied gene list. The failure was
silent: the original step-2 QC only verified that each ID matched *a*
group, never that it matched *one gene*.

**Caveats.**
- The phosphorus matched rate (52.6%) still cannot be tested against a
  background/noise rate — no unfiltered phosphorus table exists in scope
  (see Methods) — so it should be read as descriptive, not as a
  statistically confirmed enrichment the way the nitrogen result is. The
  rate also depends on the table-absent reclassification choices; the 40
  significant calls behind it do not.
  `[KG, verified 2026-09-06]` Each of the 5 phosphorus experiments' KG
  table covers only 1.5-6.0% of its strain's genome (Lin RNA-seq 1.5%,
  Fuszard proteomics 1.9%/2.8%, Martiny microarray 1.7%/6.0%), regardless
  of method — versus 76-86% for the nitrogen microarray tables on the same
  platform type. The table-absent reclassification (Martiny 2026-09-08,
  Lin 2026-09-09) counts genome-present table-absent genes as
  tested-not-significant for those 3 genome-wide experiments (the Lin step
  added 35 such calls, widening the phosphorus denominators without
  changing a significant hit); the 2 Fuszard tables remain
  detection-limited. See `gaps_and_friction.md` (2026-09-06, 2026-09-08,
  2026-09-09 entries).
- Significance criteria differ by platform (dual padj+fold-change vs.
  padj-only vs. fold-change-only with no p-value at all — see Methods),
  so raw hit-rate percentages are not fully comparable test-for-test
  across experiments; the H3 bootstrap and Fisher's test correct for this
  for nitrogen specifically, but no equivalent correction exists for
  phosphorus.
- Effective strain coverage is narrow relative to the researcher's
  original 15-strain list: the nitrogen side draws on only 2 strains
  (MED4, MIT9313), and while the phosphorus side draws on 4, MIT9312
  contributes almost no data (3 genes) because its sole in-scope
  experiment's source table only reports genes already significant in the
  original publication.
- Individual (gene x experiment) tests are not independent — the same
  gene is tested across multiple experiments and strains. The bootstrap
  and Fisher's test address this at the aggregate level for the nitrogen
  H3 comparison; the raw per-gene and per-experiment percentages
  elsewhere should be read as descriptive summaries, not independent
  trials.
- Two timepoint choices required judgment beyond the KG's own fields: the
  Tolonen nitrogen microarrays' 12h point came from the source paper (KG
  marks all timepoints `acute_stress`); and the Lin phosphorus
  experiment's 59h point (corrected from 46h on 2026-09-08) — the KG's
  `growth_phase` distinguishes the two experimental arms but the bare
  timepoint label does not, and 46h under-represented the response Lin
  reported. The Martiny MIT9313 experiment's stated "significant at 48h"
  criterion also has no matching 48h KG row for that strain (24h used,
  flagged `*` in the figure). All documented in `5_analyze/notebook.md`.
- Gene identity rests on the researcher's Cyanorak IDs plus two
  cross-checks (member-count audit, resolve-by-name). The checks caught
  6 loci; a residual failure mode remains — a wrong ID that also has no
  `gene_name` in the KG and resolves to one plausible locus would slip
  through both. `PMM719` in NATL2A is left unresolved (ambiguous;
  immaterial — no in-scope data).

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
