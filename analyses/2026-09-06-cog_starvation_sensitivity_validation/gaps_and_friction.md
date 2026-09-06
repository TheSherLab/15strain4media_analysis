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
