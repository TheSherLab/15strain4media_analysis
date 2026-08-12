# Step 1 — Research question

## Context

Follow-up to `analyses/2026-08-10-np_starvation_expression_walkthrough/`, which tested whether
the researcher's 92 N/P-acquisition genes show a differential-expression response under N/P
starvation. Separately (outside this KG), the researcher tested whether the same 92-gene list's
copy number per genome predicts a survival phenotype she measured directly: the 15 study strains
split into those that die faster under nitrogen limitation ("N-sensitive") versus those that die
faster under phosphorus limitation or at similar rates under both ("mixed"). That test found no
correlation. This analysis flips the search: instead of testing the curated 92-gene list against
the phenotype, it scans the full genome of each of the 15 strains for genes whose presence or
copy-number pattern associates with the phenotype grouping, then validates any candidates using
the prior walkthrough's differential-expression-under-starvation method.

## KG context

Grounding queries run before locking scope:

- `kg_schema` — confirmed the `Gene_in_ortholog_group` relationship (Gene → OrthologGroup /
  `GroupingClass`, source `cyanorak` or `eggnog`) as the mechanism for comparing "the same gene"
  across strains regardless of locus tag — the same mechanism the prior analysis used
  (`genes_by_homolog_group`) to resolve the 92-gene list across strains.
- `list_organisms(organism_names=[the 15 strains], verbose=true)` then re-run with `limit=15` —
  confirmed genome data exists for **all 15 strains** (`organism_type=genome_strain`,
  `gene_count` populated for each): AS9601 1,951; MED4 1,976; MIT0604 2,137; MIT1314 1,883;
  MIT1327 2,452; MIT9202 2,027; MIT9215 2,030; MIT9301 1,935; MIT9312 1,978; MIT9313 2,948;
  MIT9515 1,949; NATL1A 2,226; NATL2A 2,214; PAC1 2,370; SB 1,977 — **~32,053 gene rows total**
  across all 15. Only 7 of 15 have any expression experiments at all (matches the prior
  analysis's step-1 finding) — this constrains phase B (DE validation) but not phase A (the
  genome-wide scan), which needs only genome/gene data, present for all 15 regardless of
  experiment coverage.
- Structural surprise worth flagging now: **MIT9313's genome (2,948 genes) is substantially
  larger than the other 14 (1,883–2,452)**. A bigger genome mechanically has higher gene/copy
  counts across many families — a naive copy-number comparison could mistake "big genome" for
  "mixed-sensitivity signal." Flagged for step 3 framing, not resolved here.

## What I did

Worked through the phenotype-grouping definition and the two prior tests (this repo's 92-gene DE
walkthrough; the researcher's own out-of-KG acquisition-gene-copy-number test against the
survival phenotype) with the researcher via three clarifying questions, then ran the KG grounding
queries above before proposing and locking the question below.

Clarifying questions and answers:
1. **Would testing all genes (not just 92) take much longer?** No, not computationally — a bulk
   pull (all genes for all 15 strains, ~32k rows) replaces 92 individual lookups, and the stats
   run in seconds regardless of gene count. What it adds is a real methodological requirement
   (multiple-testing correction across however many genes get tested) and a real constraint
   (only 15 strains split 9-vs-6 limits statistical power no matter how many genes are tested) —
   both flagged for step 3, not resolved here.
2. **What should the genome-wide scan check — presence/absence, copy number, or both?** Both are
   in scope. The researcher's original acquisition-gene test used copy number; either data type
   could carry signal the other misses. Which is more informative for whatever candidates surface
   is a step-3/4 methods decision.
3. **Should this be a new analysis or an extension of the walkthrough?** New analysis — different
   research question (genotype-vs-survival-phenotype across all genes, not expression-vs-
   starvation for a curated list) and different data scope (genome content across all 15 strains,
   vs. DE data limited to the 4 strains with usable starvation experiments). Phase B explicitly
   reuses `2026-08-10-np_starvation_expression_walkthrough/4_methods/np_response.py`
   (`hit_rate`, `bootstrap_pvalue`) rather than reinventing it.

## Decisions

**2026-08-12 — New analysis, not an extension of the 92-gene walkthrough.** See clarifying
question 3 above.

**2026-08-12 — Two data types in scope for the genome-wide scan: presence/absence and copy
number.** See clarifying question 2 above. No decision yet on which is primary — deferred to
step 3/4.

**2026-08-12 — Genome-size confound (MIT9313) flagged, not resolved.** MIT9313 carries
~500–1,000 more genes than any other strain in scope. Any copy-number comparison must account
for this or risk mistaking genome size for phenotype signal. Deferred to step 3 framing.

## Locked research question

Across the 15 *Prochlorococcus* strains split by an independently measured survival phenotype
under N/P limitation (9 "N-sensitive": MED4, MIT9515, MIT9202, MIT9215, MIT0604, AS9601,
MIT9312, MIT1314, NATL1A; 6 "mixed": MIT9301, SB, NATL2A, PAC1, MIT9313, MIT1327), scan every
gene in each strain's genome — not just the researcher's 92-gene N/P-acquisition reference list,
which showed no correlation with this phenotype — for a presence/absence and/or copy-number
pattern that associates with the N-sensitive vs. mixed split. Genes that show such a pattern,
including ones with no established function, become candidate hypotheses for what drives the
phenotype. Candidates are then validated using the differential-expression-under-starvation
method built in `2026-08-10-np_starvation_expression_walkthrough/4_methods/np_response.py`
(hit-rate + bootstrap significance test), checking whether they actually respond to real
nutrient starvation.

## Decide-gate checklist

- **Outputs produced:** this notebook; analysis scaffold (`paper.md`, `gaps_and_friction.md`,
  `.gitignore`).
- **Results presented:** KG grounding (genome data for all 15 strains, ortholog-group mechanism,
  genome-size range); locked question (all above).
- **QC gate:** all 15 strains checked against `list_organisms` → all have genome data
  (`gene_count` populated, `organism_type=genome_strain`); phenotype-grouping membership
  cross-checked against the researcher's 9+6 lists — all 15 strains accounted for, no overlap
  between groups.
- **Decisions made this step:** new-analysis framing; both presence/absence and copy number kept
  in scope; genome-size confound flagged (all above, dated 2026-08-12).
- **Advance rationale:** question and scope are locked with the researcher's explicit input on
  the two-phase structure (genome-wide scan, then DE validation) and the genome-size caveat;
  ready for step 2 to enumerate what the KG actually offers for a genome-wide, cross-strain
  gene-content comparison.
