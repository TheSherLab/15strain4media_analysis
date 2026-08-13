# N/P starvation sensitivity — genome-wide gene association

## Question

Across the 15 *Prochlorococcus* strains split by an independently measured survival phenotype
under N/P limitation (9 "N-sensitive": MED4, MIT9515, MIT9202, MIT9215, MIT0604, AS9601,
MIT9312, MIT1314, NATL1A; 6 "mixed": MIT9301, SB, NATL2A, PAC1, MIT9313, MIT1327), scan every
gene in each strain's genome — not just the researcher's 92-gene N/P-acquisition reference list,
which showed no correlation with this phenotype — for a presence/absence and/or copy-number
pattern that associates with the N-sensitive vs. mixed split. Genes that show such a pattern,
including ones with no established function, become candidate hypotheses for what drives the
phenotype. Candidates are then validated using the differential-expression-under-starvation
method built in `2026-08-10-np_starvation_expression_walkthrough/4_methods/np_response.py`,
checking whether they actually respond to real nutrient starvation.

## Background

All 15 study strains have genome data in the KG regardless of experiment coverage: 32,053 genes
total (1,883–2,948 per strain). Cross-strain gene identity is resolved via the KG's
`Gene_in_ortholog_group` relationship, which carries two independent grouping systems: Cyanorak
(curated specifically for cyanobacteria, one group per gene) and Eggnog (automated, broader,
up to 3 nested groupings per gene at increasing taxonomic breadth). 29,356 of the 32,053 genes
(91.6%) have a Cyanorak group; the remaining 2,697 are resolved via Eggnog's tightest grouping
level where available (2,179 genes), leaving 518 "true orphan" genes with no cross-strain group
in either system (of which 389 are functionally uncharacterized). Cyanorak's 5,003 distinct
groups average only 5.73 of the 15 strains each, confirming substantial real presence/absence
variation across strains rather than a largely-universal gene set. See
`2_kg_selection/notebook.md` for the full extraction, assignment-rule rationale, and the
per-strain true-orphan breakdown (uneven across strains — MIT1314, MIT1327, and MIT9313 notably
higher than the rest, not yet interpreted).

## Methods

**Clade confound.** The sensitivity phenotype is not independent of strain ancestry: HLI/HLII-
clade strains lean heavily N-sensitive (8 of 10 HL strains) and LLIV-clade strains are entirely
mixed (2 of 2), so a naive genome-wide scan risks surfacing ecotype-driven genes rather than
sensitivity-specific ones. Every strain's clade (HLI, HLII, LLI, LLIV) is collapsed to a
high-light/low-light (HL/LL) ecotype split — HL: MED4, MIT9515, AS9601, MIT0604, MIT1314,
MIT9202, MIT9215, MIT9301, MIT9312, SB (10 strains); LL: NATL1A, NATL2A, PAC1, MIT9313, MIT1327
(5 strains) — and every gene is tested against **both** the sensitivity grouping (9 N-sensitive
vs. 6 mixed) and this ecotype grouping, using the same statistical test each time.

**Statistical tests.** Presence/absence per ortholog group: Fisher's exact test on the 2x2 count
of strains carrying vs. not carrying the group (appropriate for the small strain counts here).
Copy number per ortholog group (where variable): Mann-Whitney U test, with copy number normalized
per strain by that strain's total gene count first, to prevent MIT9313's larger genome from
mechanically inflating its counts. Benjamini-Hochberg FDR correction is applied to both the
sensitivity and ecotype p-values across all groups tested (several thousand simultaneous tests);
significant means adjusted p < 0.05. This FDR correction is the noise control for this analysis —
unlike the prior walkthrough, there is no separate background gene pool, since every gene is
tested at once rather than a curated subset compared against a background. Both the raw and
FDR-adjusted p-value are retained in the output for every test, not just the adjusted value.

**Candidate definition.** Every gene keeps both test results (sensitivity, ecotype) side by side
in the output table — nothing is filtered out. A gene is tagged `candidate` when its sensitivity
result is significant and its ecotype result is not; genes significant for both (or for ecotype
only) remain fully visible in the table, just labeled differently. Alongside each test's p-value,
**direction** is recorded independently for both groupings — which side (N-sensitive vs. mixed;
HL vs. LL) has the higher presence rate or higher median normalized copy number — so a candidate's
readout states not just "significant" but which way it points.

**Controls.** No sensitivity-specific positive control exists (that is what this analysis
searches for, and naming one from memory would violate the KG-as-sole-source rule this project
runs on). The ecotype test itself serves as the method sanity check: the clade cross-tabulation
already shows a real structural ecotype signal in this strain set, so reliably recovering strong,
FDR-significant ecotype hits validates the pipeline before it is trusted on the unvalidated
sensitivity question.

**Methods module.** The two tests, FDR correction, and classification logic are implemented in
`4_methods/np_sensitivity_scan.py` (`presence_absence_test`, `copy_number_test`,
`benjamini_hochberg`, `classify_significance`; scipy's `fisher_exact` / `mannwhitneyu`,
statsmodels' `multipletests` with `method="fdr_bh"`). Every test function takes the two strain
groups to compare as arguments, so the same code path runs both the sensitivity and ecotype
comparisons. Copy number is judged "variable" (and so eligible for the Mann-Whitney test) when
its per-strain count — absent strains counted as 0 — takes more than one distinct value across
the strains being compared; FDR correction is applied separately per (test type × grouping)
combination (4 corrections total), not pooled, since presence/absence and copy-number test
different, only partly-overlapping sets of ortholog groups. Both co-defined with the researcher
(`4_methods/notebook.md`).

Verified against 8 hand-built toy cases (clear/null presence signal, a degenerate
universal-presence edge case, clear/absent copy-number variation, variability-via-absence, FDR
wiring, all 4 classification outcomes) before running on real data — all 8 passed. Worked example
on real data: ortholog group `cyanorak:CK_00000001` (`rpoD8`, RNA polymerase sigma factor type
II), carried by 14 of 15 strains (absent only in MIT1314, 2 copies in MIT1327 and MIT9313, 1
copy elsewhere) — presence/absence p=1.0 on both groupings, copy-number p=0.33 (sensitivity) and
p=0.86 (ecotype), no signal on either axis. This group was picked by a mechanical rule (first
ortholog group, sorted by ID, with both presence/absence and copy-number variation), not for
biological relevance, so the sanity check here is methodological rather than biological: the
near-universal-presence case (14/15) correctly resolves to p=1.0 rather than erroring, confirming
the degenerate-table handling verified on toy data also holds on real data.

## Results

## Discussion

## References
