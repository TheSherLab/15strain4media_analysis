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

## Results

## Discussion

## References
