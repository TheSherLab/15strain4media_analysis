# COG-based starvation-sensitivity candidates — in vivo expression validation

## Question

The researcher's comparative-genomics work (`Dataset 3.xlsx`, "Significant
NorMixed COGs") identified **41 COGs** whose cross-strain **copy number**
associates with a strain-level "starvation sensitivity" phenotype, tagged
`N` (13 COGs, nitrogen-linked) or `mixed` (28 COGs, more than one nutrient
/ unclear). This analysis asks whether that genomic signal has in-vivo
expression support: do the genes making up these COGs show a
differential-expression response under literal nitrogen or phosphorus
starvation vs. replete medium, in axenic, uninfected Prochlorococcus,
restricted to the researcher's 15 study strains? Reported for all 41 COGs
pooled and for the `N` / `mixed` subsets, per nutrient. A second question
asks whether the nitrogen response is distinguishable from a random
background rate (bootstrap + Fisher's exact), with this 41-COG gene set
and the prior analysis's 92-gene N/P-acquisition list both excluded from
that background.

Evidence sources: the 10 experiments (5 nitrogen, 5 phosphorus) validated
in the prior `2026-08-10-np_starvation_expression_walkthrough/` analysis.

*Scope history: originally scoped to both nutrients; narrowed to nitrogen
only during step 2 (2026-09-06) when raw phosphorus evidence for the 41
COGs was near-zero; phosphorus restored (2026-09-09) once the walkthrough
analysis's later Martiny reclassification made a flat phosphorus result
reportable. See `1_question/notebook.md`.*

## Background

This analysis reuses the 10-experiment scope validated in the walkthrough
analysis (5 nitrogen — MED4, MIT9313; 5 phosphorus — NATL2A, MIT9312,
MED4, MIT9313), re-confirmed present in the live KG.

**Gene resolution (rebuilt 2026-09-09).** An earlier pass resolved each
COG through a single researcher-supplied Cyanorak ID. Two identity checks
(`2_kg_selection/scripts/04_paralog_audit.py`,
`05_name_vs_cyanorak_crosscheck.py`) showed each such ID is a real member
of its COG but only one of several — and sometimes a poor representative
(`COG0443` "DnaK" resolved to a DUF3181 "conserved hypothetical protein"
while the actual `dnaK1/2/3` chaperones sit in three other Cyanorak
groups). Because the researcher's significance test was run on COG **copy
number** — every gene carrying the COG number in a strain — resolution was
rebuilt from her own whole-genome annotation
(`WhloeGenome_AllStrains_Concated.xlsx`): for each COG, in each strain,
every CDS tagged with that COG number, its Cyanorak cluster mapped to a
multiomics-KG locus. This grows the gene universe from ~95 loci to **427
(COG × strain × locus) rows across the 4 evidence strains** (389 distinct
loci). `COG0443` now correctly includes `dnaK1/2/3`; `COG0188` includes
the real `gyrA` alongside the topoisomerase-IV paralog; and so on.

Of the 41 COGs, **37 have differential-expression evidence** in ≥1
in-scope experiment (10 `N`, 27 `mixed`). The 4 with none (`COG0367` AsnB,
`COG1454` EutG, `COG3206` GumC, `COG3727` Vsr) resolve only in strains
whose in-scope experiment does not cover them. The gene set is uneven:
15 COGs resolve to 1 gene per strain, 3 (`COG0477`, `COG0697`, `COG0845`)
to 10+ — `COG0477` (every MFS permease) alone is 25–34 genes per strain.
11 Cyanorak clusters carry more than one of the 41 COG numbers
(multi-domain genes such as `zipN`); these are kept under each COG,
matching the copy-number test.

## Methods

**Hypothesis:** the 41-COG list — using the full per-strain COG
membership — shows a nitrogen-starvation response elevated above a random
background, tested for all 41 pooled and for the `N` (13) / `mixed` (28)
Direction subsets. Phosphorus reported the same way but descriptively. No
positive-control subset: this list is candidate genes without established
literature confirmation.

Significance criteria for the 10 experiments are reused unchanged from the
walkthrough analysis (same platforms, same empirically-derived
padj/log2FC thresholds). One representative starvation timepoint per
experiment, reused from the walkthrough's final (2026-09-08) table —
including **Lin at 59h** and Martiny at 48h (MED4) / 24h (MIT9313, flagged
— no 48h row). The two Martiny genome-wide phosphorus microarrays get the
walkthrough's reclassification: a COG-gene present in the strain but
absent from Martiny's q<0.05 table is counted `not_significant`, not "no
data". Not applied to Lin (curated operon table) or Fuszard (iTRAQ fixed
detected-protein list).

The negative/background pool (`3_analysis_framing`) reuses the walkthrough's
4 unfiltered nitrogen experiments and single-timepoint rule, excluding
this analysis's full 41-COG gene set, the walkthrough's 92-gene list (read
live, so it carries that analysis's 2026-09-08 corrections), and the same
13 nitrogen-keyword text matches — 6,664 background genes, overall
significant rate 19.1%. **Nitrogen only:** the phosphorus tables remain
pre-filtered / detection-limited even after the Martiny reclassification,
so there is no unbiased phosphorus population to resample.

`hit_rate()` and `bootstrap_pvalue()` (Fisher's exact alongside) are
reused unchanged from the walkthrough's `4_methods/np_response.py`,
verified in step 4 against a real driving example (`COG1403`/McrA).

## Results

**Nitrogen** (`5_analyze/data/02_hit_rate_results.csv`; background 19.1%):

| Group | Tests | Up | Down | % significant |
|---|---|---|---|---|
| all 41 pooled | 438 | 37 | 31 | **15.5%** |
| Direction N (13) | 99 | 12 | 11 | 23.2% |
| Direction mixed (28) | 339 | 25 | 20 | 13.3% |

**Nitrogen background comparison** (10,000-iteration bootstrap + Fisher's
exact):

| Group | Observed | Bootstrap null mean (max) | Bootstrap p | Fisher OR (p) |
|---|---|---|---|---|
| pooled (41) | 15.5% | 18.3% (25.2%) | 0.93 | 0.78 (0.067) |
| Direction N (13) | 23.2% | 22.9% (39.8%) | 0.44 | 1.28 (0.30) |
| Direction mixed (28) | 13.3% | 16.8% (27.4%) | 0.97 | 0.65 (**0.0067**) |

**Neither the pooled list nor either subset clears background.** The `N`
subset trends above (23.2%) but is squarely inside the bootstrap null; the
`mixed` subset is significantly *below* background.

**Phosphorus** (descriptive, no background test): **5.3%** significant
(13/246), 8 up / 5 down — flat.

**Robustness** (`6_evaluate/data/01_sensitivity_summary.txt`): the
nitrogen headline is stable — 13.9% with the 3 biggest COGs
(`COG0477`/`COG0697`/`COG0845`) removed (bootstrap p=0.97), 14.9% with
each COG collapsed to one value per experiment. Every cut sits at or below
background.

**Main figure:** `5_analyze/figures/01_gene_experiment_heatmap.png` — 10
experiments (rows, N block then P block) × 41 COGs (columns), the COGs
bracketed and labelled by the researcher's "General annotation" functional
category (Quality control – DNA / Protein level, Amino acid & mixotrophy
biosynthesis, LPS, Membrane linker, Exopolysaccharide, Biofilm/Attachment,
Energy production, Translation/transcription/signalling, RTX toxins,
mixed); COG label colour = Direction. Coloured cell text is `k / n / m` —
`k` significant / `n` genes measured in that experiment / `m` genome
copies in that strain (`m` shown only when `n < m`); colour = dominant
direction, intensity scales with `k/n`; walkthrough palette; COGs absent
from a strain's genome hatched. `[KG]` The field is mostly grey
(tested, no gene significant); the Weissberg RNA-seq row is the hottest
(that experiment has a high genome-wide baseline, established in the
walkthrough); the entire phosphorus block is grey/pale.
`5_analyze/figures/02_pct_significant_vs_background.png` — the hit-rate
bars against the 19.1% background line.

Top nitrogen "responsive" COGs (`6_evaluate/data/02_per_cog_nitrogen_fractions.csv`):
`COG1596` (Wza, 2/5), `COG1787` (Mrr, 1/3), `COG2192` (carbamoyltransferase,
2/7), `COG0188` (gyrase/topo, 7/26) — all `N`-tagged, all small, all with
balanced up/down.

## Discussion

**The 41-COG candidate list shows no in-vivo expression response to
nutrient starvation above background.** The pooled nitrogen hit rate
(15.5%) is below the random background (19.1%); bootstrap p = 0.93,
Fisher's exact p = 0.067. The `N`-tagged subset trends above background
(23.2%) but is not statistically distinguishable from a random
size-matched draw (bootstrap p = 0.44). The `mixed` subset is
significantly *below* background (Fisher p = 0.0067). Phosphorus is flat
(5.3%). The result is robust to dropping the largest COGs and to counting
COGs rather than genes.

`[interpretation]` These COGs are general cellular-maintenance machinery —
chaperones (DnaK, DnaJ, CbpA), DNA repair (gyrase, Vsr, McrA, RM systems),
membrane transport (MFS/DMT/ABC), LPS/exopolysaccharide biosynthesis. As a
class they are not transcriptionally or translationally mobilised by
nitrogen or phosphorus starvation any more than a random gene, and the
`mixed` half slightly less so — consistent with housekeeping functions
that are held roughly constant while the cell reallocates resources.

**What this means for the comparative-genomics finding.** `[interpretation]`
The 41 COGs were selected because their cross-strain gene *copy number*
tracks the starvation-sensitivity phenotype. This analysis tested a
different, independent prediction — that starvation changes these genes'
*expression* — and it does not hold. This does not refute the copy-number
association: gene-dosage differences between strains can matter for a
phenotype without the genes being starvation-*regulated* within a strain.
But the two lines of evidence are not mutually reinforcing, and should be
reported as separate observations rather than one supporting the other.

**Rebuilding the gene set mattered, and did not rescue the signal.** The
2026-09-09 reopen replaced one representative Cyanorak ID per COG with the
full per-strain COG membership — the same gene set the copy-number test
used, and the correct one for this comparison. It added the real DnaK/DnaJ
chaperones and other family members that the single-ID route had missed.
The pooled result got marginally *more* negative on the larger set (16.2%
→ 15.5%), because those added chaperones are not strongly
starvation-regulated in these experiments either.

**Caveats.**
- The pooled rate mixes very unequal COGs — `COG0477` alone is ~112 of
  438 nitrogen tests. The "use every gene" choice matches the copy-number
  test's gene set; the sensitivity checks (drop-big-COGs, per-COG unit)
  are why the headline is trustworthy despite this.
- Nitrogen draws on 2 strains (MED4, MIT9313); phosphorus effectively 2
  (MED4, MIT9313 via Martiny), with Lin/Fuszard contributing a handful of
  table genes. Narrow relative to the 13–15 strains the copy-number test
  spanned.
- A few genes carry 2–3 COG numbers (multi-domain proteins); their
  expression is counted once per COG, matching the copy-number test — a
  small non-independence.
- Individual (gene × experiment) tests are not independent; the bootstrap
  and Fisher's test address this at the aggregate for nitrogen only.
- The phosphorus 5.3% is descriptive — no valid phosphorus background
  exists, so it is not a confirmed "below background".
- Significance criteria differ by platform (padj + fold-change vs
  padj-only vs fold-change-only), carried unchanged from the walkthrough;
  corrected for by the bootstrap/Fisher for nitrogen only.

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
