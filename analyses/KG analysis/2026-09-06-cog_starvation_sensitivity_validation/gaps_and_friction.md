# Gaps and friction log

Append-only. Each entry: date, short name, what happened, workaround/impact if any.

**2026-09-06 — Phosphorus scope dropped after step 2's evidence check
showed it was effectively unusable for this gene list.** Step 1 locked
"test both nitrogen and phosphorus" per the researcher's explicit choice.
Step 2 found 32 of 41 COGs have DE evidence, but 32 of those 32 come from
the nitrogen side and only 1 from phosphorus; NATL2A resolves 26 of 41
COGs to a locus tag but contributes zero evidence, because none of those
26 loci happen to fall inside its two in-scope experiments' own narrow,
pre-filtered gene tables (1.5-2.8% of genome — the same phosphorus-table
narrowness quantified in the prior walkthrough analysis's post-closure
verification). This is a second, independent case of the same
methodology gap: any gene list checked against these particular
phosphorus source tables will see this same near-total blackout unless it
happens to overlap the source papers' own pre-selected significant-gene
lists. Impact: reopened step 1 (`1_question/notebook.md`, "Reopened
2026-09-06") to drop phosphorus from scope; the 5 phosphorus experiments
stay documented in step 2's data as the record of what was checked. Future
analyses reusing this KG's phosphorus Prochlorococcus experiments should
expect the same problem and budget for it up front rather than discovering
it after gene resolution.

**2026-09-09 — A single Cyanorak ID per COG is a poor proxy for the COG's
gene set; the copy-number test used the whole membership.** The 2026-09-06
resolution used one researcher-supplied Cyanorak ID per COG. Two identity
checks (`2_kg_selection/scripts/04_paralog_audit.py`,
`05_name_vs_cyanorak_crosscheck.py`) showed every supplied ID is a *real*
member of its COG but only one of several — and sometimes a bad pick:
`COG0443` ("DnaK") pointed at `CK_00000457`, a DUF3181 "conserved
hypothetical protein" (COG category S), while `dnaK1/2/3` sit in three
other Cyanorak groups the single-ID route never saw. Root cause: the
researcher's comparative-genomics significance test was run on COG **copy
number** (every gene carrying the COG number in a strain), so the
expression validation must use the same gene set. Fix: resolution rebuilt
from the researcher's own whole-genome annotation
(`WhloeGenome_AllStrains_Concated.xlsx`) — `00_wholegenome_cog_members.csv`,
`02_resolve_genes.py`. Gene universe grew ~95 → 427 (COG × strain × locus)
rows in the 4 evidence strains. Impact: any analysis validating a
COG-level or ortholog-group-level comparative-genomics result against
per-gene data must resolve to the *full* group/COG membership, not a
representative — and should get that membership from the same annotation
the original test used, not re-derive it. A representative-ID shortcut
silently changes which gene is tested.

**2026-09-09 — COG-level heatmap: "most significant gene in the COG"
inflates big COGs.** First step-5 heatmap pass coloured each (COG,
experiment) cell by its single most-significant gene. With COGs of 25–34
genes (COG0477 = every MFS permease), almost every big-COG cell rendered
as significant, because at least one of 30 genes clears threshold by
chance. Fixed: cell shows the *fraction* of the COG's tested genes that
are significant (`k/n`), colour intensity scaled by it — same fix pattern
as the walkthrough's Figure 4. Impact: any per-group summary cell over a
variable-size group needs a rate, not an extremum; flag whenever a
"collapse to one cell" figure sits on top of groups that range from 1 to
30+ members.

**2026-09-09 — Phosphorus restored after the walkthrough's later Martiny
decision.** This analysis dropped phosphorus on 2026-09-06 because raw DE
evidence for the 41 COGs was near-zero on that side. The walkthrough
analysis then (2026-09-08) adopted the researcher's argument that a
genome-wide microarray filtered purely on significance (Martiny et al.
2006) yields real "not significant" calls for table-absent genes. Applying
that here restored a usable — though flat — phosphorus result. Impact:
scope decisions driven by "no evidence" should be revisited when an
upstream analysis changes how "tested-absent" is classified; the
`table_scope` semantics are still evolving across analyses.

**2026-09-09 — "table-absent = not significant" extended from Martiny to
Lin; the walkthrough's own "not applied to Lin" note was over-cautious.**
The 2026-09-08 walkthrough reopen applied the table-absent
reclassification to the two Martiny microarrays but explicitly *not* to
Lin et al. 2015, reasoning that Lin's 34-gene table "is not a
threshold-defined set". A live `list_experiments` check
(`publication_doi=10.1111/1462-2920.13104`) shows the in-scope Lin
uninfected experiment is tagged `table_scope = significant_only` and is
genome-wide RNA-seq — so a target gene with a NATL2A locus but no Lin row
was measured and not significant, exactly the Martiny logic. The earlier
note conflated "few of the 34 listed genes are significant at every
timepoint" (true, and irrelevant) with "the table isn't significance-
filtered" (false — the KG tag says it is). Fix: Lin added to
`RECLASSIFY_ABSENT_AS_NOT_SIG` in both analyses' `01_extract_target_de.py`
(2026-09-09, researcher-directed). COG phosphorus hit rate 5.3% → 3.8%;
walkthrough P-matched rate recomputed. Impact: when a prior analysis
hedges a `table_scope`-driven call, re-check the KG's own `table_scope`
field before carrying the hedge forward — it is authoritative and the
hedge may not be.

**2026-09-06 — `differential_expression_by_gene(experiment_ids=[])` is a
no-op filter, not "match nothing".** While building step 4's driving
example, iterating per-organism (not per-experiment) produced an empty
`experiment_ids` list for every strain with no in-scope nitrogen
experiment (13 of the 15 study strains). The call didn't return zero rows
as expected — it silently returned whatever unrelated experiment data
existed for that organism/locus elsewhere in the KG (e.g. a 1-row hit for
AS9601 from a completely different, out-of-scope treatment experiment),
corrupting the driving example's row count before it was caught by the
hand-tally cross-check. Confirmed directly: `experiment_ids=[]` behaves
identically to `experiment_ids=None` (no filter applied at all), not an
empty match set. Workaround: never pass an empty `experiment_ids` list —
filter to only the organisms that actually have an in-scope experiment
*before* calling the function, and assert the list is non-empty at the
call site (done in `4_methods/scripts/01_driving_example.py`). Impact:
any script that derives `experiment_ids` from a per-organism lookup
(rather than iterating per-experiment the way the prior analysis's
extraction scripts do) is at risk of this same silent scope leak — worth
a standing caution for future analyses' extraction scripts.

**2026-09-09 — A per-group summary cell that reports one `k` and one
colour hides within-group direction disagreement.** Figure 1 collapses
each COG's genes in a cell to `k` significant / `n` measured / `m` copies,
coloured by the majority direction. The researcher asked whether any cell
mixes up- and down-regulated genes: 11 do (all nitrogen, all MED4;
`COG0477` alone in 3 experiments, up to 3 up / 6 down). The single colour
showed only the majority and `k` only the total — a reader could not see
the split. Fix: those 11 cells are drawn **orange** (a third fill state,
not the dominant-direction colour) and `k` is written `Nup Mdn`
(`5_analyze/scripts/05_within_cell_direction_splits.py` enumerates them;
`03_figures.py` / `04_heatmap_pptx.py` render them). Same shape as the
earlier "extremum vs rate" fix on this figure — a collapse-to-one-cell
summary over a heterogeneous group needs to carry enough detail that the
collapse is visible, not just a point estimate. Numbers unchanged; the
per-COG sensitivity check already bounds the effect on the headline.
