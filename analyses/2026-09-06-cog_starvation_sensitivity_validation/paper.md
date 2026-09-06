# COG-based starvation-sensitivity candidates — in vivo expression validation

## Question

Do the 41 COGs in the researcher's "Significant NorMixed COGs" list
(`Dataset 3.xlsx`) — genes whose cross-strain copy-number/presence pattern
her comparative-genomics work associates with starvation sensitivity,
tagged `N` (12 genes, nitrogen-linked), `mixed` (28 genes, associated with
more than one nutrient/unclear), or `unlabeled` (1 gene) — show an in vivo
differential-expression response under literal **nitrogen** starvation vs.
replete medium, in axenic, uninfected Prochlorococcus restricted to the
researcher's 15 study strains? Reported for all 41 genes pooled and for
the `N`/`mixed`/`unlabeled` subsets separately. Evidence sources: the 5
nitrogen experiments (RNA-seq, proteomics, microarray) validated in the
prior `2026-08-10-np_starvation_expression_walkthrough/` analysis. A
second question asks whether any observed response is statistically
distinguishable from a background/noise rate (bootstrap + Fisher's exact),
with both this 41-gene list and the prior analysis's 92-gene list excluded
from that background population.

*(Originally scoped to test both nitrogen and phosphorus; reopened during
step 2 and narrowed to nitrogen-only after the phosphorus side proved to
carry essentially no usable evidence for this gene list — see
`1_question/notebook.md`, "Reopened 2026-09-06".)*

## Background

This analysis reuses the 10-experiment scope validated in the prior
`2026-08-10-np_starvation_expression_walkthrough` analysis (5 nitrogen
experiments — MED4, MIT9313; 5 phosphorus experiments — NATL2A, MIT9312,
MED4, MIT9313), re-confirmed present in the live KG. The 41 COGs resolve
via the researcher's own COG-to-Cyanorak-ID mapping (verified complete:
41/41 COGs, 364 gene rows across the 15 study strains, zero unmatched
groups/strains) — all 41 resolve to at least one locus tag in at least one
strain. Checking each resolved gene against the 10 in-scope experiments'
differential-expression rows narrows this to **32 of 41 COGs with actual
expression evidence**, and that evidence is heavily one-sided: all 32 come
from the nitrogen side (MED4/MIT9313); only 1 comes from phosphorus
(MIT9312). NATL2A resolves 26 of the 41 COGs to a locus tag but
contributes zero evidence — none of those 26 loci fall inside either of
NATL2A's two in-scope experiments' own narrow, pre-filtered gene tables
(consistent with the 1.5-2.8%-of-genome phosphorus coverage already
established in the prior analysis). See `2_kg_selection/notebook.md` for
the full resolution and evidence tables.

## Methods

## Results

## Discussion

## References
