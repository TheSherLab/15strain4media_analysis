# Step 1 — Research question

## Context

Follow-on from the lightweight triage analysis
`analyses/2026-07-13-n_p_genes_in_vivo_experiments/`, which checked all 91
genes in the researcher's "N and P acquisition genes" spreadsheet
(`Dataset 2.xlsx`) for the mere presence of in vivo expression data in the
KG, and found 54 with evidence (29 nitrogen-annotated, 25
phosphorus-annotated — recounted precisely in step 2 below; see that analysis's
`data/n_p_genes_experimental_evidence_final.csv`,
`has_expression_evidence == True`). That analysis stopped at "data exists
for this gene" — it did not test whether the gene actually responds to
nutrient starvation. This analysis picks up that question as a full
6-step study.

## What I did

Clarifying dialogue with the researcher to lock the question, plus KG
grounding queries to check what nitrogen/phosphorus starvation
experiments the KG holds for Prochlorococcus before locking scope.

**Dialogue summary:**
- Starting point: the 54-gene evidence list from the prior analysis, split
  29 nitrogen-annotated / 25 phosphorus-annotated.
- Researcher's question: do these genes actually change expression in
  vivo under N and P starvation, and — expanding on the first pass —
  how do phosphorus-annotated genes behave under nitrogen starvation and
  vice versa (cross-nutrient response), not just the matched-nutrient
  case.
- Culture background: KG has both axenic (pure-culture) and
  coculture-with-*Alteromonas macleodii* nitrogen-starvation experiments.
  Researcher chose **axenic only**, consistent with the prior analysis's
  handling (coculture rows were dropped before computing log2FC there).

## KG context

Grounding queries (`list_experiments`, run 2026-08-02, no organism filter
beyond treatment_type — Prochlorococcus strains identified via the
`by_organism` breakdown):

- `list_experiments(treatment_type=['nitrogen'], summary=True)` →
  16 experiments total. By organism: Prochlorococcus MED4 (8), MIT9313
  (3), SS120/CCMP1375 (1) — 12 Prochlorococcus experiments, all tagged
  `axenic`; the remaining 4 are tagged organism "Alteromonas macleodii
  HOT1A3" with background factor `coculture` (excluded per the scope
  decision above). By omics type: MICROARRAY 6, RNASEQ 5, PROTEOMICS 5.
- `list_experiments(treatment_type=['phosphorus'], summary=True)` →
  9 experiments total, all Prochlorococcus, all `axenic` (no coculture
  phosphorus experiments exist in this KG build). By organism: NATL2A
  (3), MIT9301 (2), MIT9313 (1), MED4 (1), SS120/CCMP1375 (1), MIT9312
  (1). By omics type: PROTEOMICS 3, RNASEQ 2, MICROARRAY 2, METABOLOMICS 2.

This matches the prior analysis's finding that only a handful of the
15-strain paper set (MED4, MIT9301, MIT9312, MIT9313, NATL2A, plus
AS9601 for other treatments) actually carry nitrogen/phosphorus
significant-DE data in the current KG build — step 2 will enumerate the
exact experiment list and publications.

## Research question (locked)

Do the 54 N/P-acquisition genes with confirmed in vivo evidence (from
the prior analysis) show a significant expression change under nutrient
starvation in axenic Prochlorococcus — both **matched** (nitrogen-
annotated genes under nitrogen starvation, phosphorus-annotated genes
under phosphorus starvation) and **cross** (nitrogen-annotated genes
under phosphorus starvation, phosphorus-annotated genes under nitrogen
starvation)?

## Decisions

- **2026-08-02** — Scope restricted to axenic experiments only (excludes
  4 coculture-with-*Alteromonas* nitrogen experiments), for consistency
  with the prior triage analysis and to avoid a coculture-partner
  confound.
- **2026-08-02** — Question expanded mid-dialogue from matched-nutrient
  only to include the cross-nutrient case (N genes under P starvation,
  P genes under N starvation), at researcher's request.

## Decide-gate checklist

- **Outputs produced** — `paper.md` (skeleton), `gaps_and_friction.md`
  (header only), `.gitignore`, `1_question/notebook.md` (this file). No
  scripts this step (step 1 is conversation + grounding, not computation).
- **Results presented** — locked research question above; KG grounding
  counts (16 nitrogen / 9 phosphorus experiments, organism and omics
  breakdowns) shown inline above, matching what was shown to the
  researcher in chat.
- **QC gate** — confirmed via `list_experiments` that both nitrogen and
  phosphorus axenic Prochlorococcus experiments exist in the KG in
  sufficient number to support the matched and cross comparisons →
  question is answerable with current KG data.
- **Decisions made this step** — see Decisions above.
- **Advance rationale** — question is locked, scope (axenic-only,
  matched + cross) agreed with researcher, and grounding confirms the KG
  holds enough nitrogen/phosphorus experiments to proceed to step 2 (KG
  entries).
