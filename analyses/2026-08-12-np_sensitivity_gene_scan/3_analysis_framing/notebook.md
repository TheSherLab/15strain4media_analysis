# Step 3 — Analysis framing

## Context

Step 2 built the full gene x ortholog-group assignment table (32,053 genes; 29,356 Cyanorak,
2,179 Eggnog-rank1 fallback, 518 true orphans). This step turns that raw material into an actual
statistical test, and confronts a real confound the raw data surfaced: the phenotype grouping is
not independent of strain ancestry.

## What I did

**Clade-confound discussion** (co-defined with the researcher, no script — a real judgment call,
not a computation): cross-referencing each strain's KG-recorded `clade` against the sensitivity
grouping showed HLI/HLII-clade strains lean heavily N-sensitive and LLIV-clade strains are
entirely mixed (see Results). The researcher confirmed this was already known and had been
handled before by testing association against *both* the sensitivity phenotype and the ecotype
(HL/LL) grouping, then distinguishing sensitivity-specific genes from ecotype-driven ones —
without discarding the ecotype-associated genes, just tagging them separately so they stay
accessible.

**Built the strain-group membership table** (`scripts/01_build_strain_groups.py`): pulled each
strain's `clade` from the KG (`list_organisms`, not assumed from memory) and collapsed it to the
standard two-ecotype split (HLI+HLII → HL, LLI+LLIV → LL), alongside the existing
sensitivity-phenotype grouping.

## Results

**Clade vs. sensitivity cross-tabulation** — the confound that motivated this step:

| Clade | Strains | N-sensitive | Mixed |
|---|---|---|---|
| HLI | 2 | 2 (100%) | 0 |
| HLII | 8 | 6 (75%) | 2 |
| LLI | 3 | 1 | 2 |
| LLIV | 2 | 0 | 2 (100%) |

**Strain group membership** (full table: `data/01_strain_group_membership.csv`):

| Ecotype | Strains (10 + 5 = 15) |
|---|---|
| HL (HLI+HLII) | MED4, MIT9515, AS9601, MIT0604, MIT1314, MIT9202, MIT9215, MIT9301, MIT9312, SB |
| LL (LLI+LLIV) | NATL1A, NATL2A, PAC1, MIT9313, MIT1327 |

**Sensitivity x ecotype cross-tab** (confirms the two groupings are correlated but not
identical — real room for a gene to differ on one axis without differing on the other):

| | HL | LL |
|---|---|---|
| N-sensitive | 8 | 1 |
| Mixed | 2 | 4 |

## Framing

**Hypothesis (prose):** some gene(s) beyond the researcher's 92-gene N/P-acquisition list — which
already showed no correlation with the sensitivity phenotype — show a presence/absence or
copy-number pattern across the 15 strains' genomes that associates with the N-sensitive vs. mixed
split specifically, distinguishable from both chance and from strain ancestry (ecotype).

**Test design (co-defined):**
- **Presence/absence**, per ortholog group (Cyanorak or Eggnog-rank1 assigned; true orphans
  excluded from this test — see step 2): Fisher's exact test on the 2x2 count of strains
  carrying vs. not carrying the group, run once against the sensitivity split (9 vs. 6) and once
  against the ecotype split (10 vs. 5).
- **Copy number**, per ortholog group with variable copy count across strains: Mann-Whitney U
  test on copy-number distributions, same two groupings. Copy number is normalized per strain by
  that strain's total gene count (copies per 1,000 genes) before comparison, to keep MIT9313's
  larger genome from mechanically inflating its copy counts; any candidate whose signal only
  survives on raw (non-normalized) counts is flagged as a genome-size artifact, not kept as a
  clean candidate.
- **Direction, alongside every test result:** a p-value says whether a group differs between the
  two strains-groups, not which way. For presence/absence, direction compares the presence *rate*
  (fraction of strains carrying the group) between the two sides — e.g. "enriched in N-sensitive"
  if the rate is higher there than in mixed. For copy number, direction compares the median
  normalized copy number between the two sides — e.g. "higher in HL". Recorded for both the
  sensitivity test and the ecotype test independently, so a candidate's full readout is: does it
  differ (p-value), how strongly (FDR-adjusted p-value), and which side it leans toward, for each
  of the two groupings.
- **Multiple-testing correction:** Benjamini-Hochberg FDR applied to both the sensitivity and
  ecotype p-values across all groups tested; "significant" means adjusted p < 0.05. This plays
  the same noise-control role the walkthrough's background/bootstrap comparison played — there is
  no separate background pool here because every gene is tested at once, not a curated subset
  against a background. **Both the raw (pre-correction) and FDR-adjusted p-value are kept in the
  output table for every test**, not just the adjusted one — so the correction's effect on any
  given gene stays checkable rather than hidden behind a single pass/fail number.
- **Candidate tag, not a filter:** every gene keeps both test results (sensitivity, ecotype) in
  the same output table. A gene is tagged `candidate` when sensitivity is significant *and*
  ecotype is not. Genes significant for both, or for ecotype only, remain in the full table,
  labeled accordingly — nothing is discarded.

**Positive control:** no KG-verified gene is already known to drive N/P sensitivity specifically
— that is what this analysis is searching for, so naming one from memory would break the
KG-is-sole-source rule this whole project runs on. Instead, the **ecotype test itself is the
sanity check**: the clade cross-tabulation above already shows a real, structural ecotype signal
exists in this strain set, so if the pipeline reliably recovers strong, FDR-significant ecotype
hits, that validates the method before it's trusted on the harder, unvalidated sensitivity
question.

**Negative control / background:** none separate from the FDR correction (see above) — explicitly
different from the prior walkthrough, where a curated target list needed an external background
pool. Here every gene is tested, so the multiple-testing correction across the full set of tests
is the background control.

**Expected outcome (operational):** a (possibly empty) list of ortholog groups tagged
`candidate` — sensitivity-significant, ecotype-not-significant, in either the presence/absence or
copy-number test. Each candidate carries its Fisher's/Mann-Whitney statistics *and direction* for
both groupings (e.g. "enriched in N-sensitive, no ecotype direction"), its `gene_category`
(flagging unknown-function ones), and enough detail to hand to step 4/5's validation phase
(checking whether it responds to real nutrient starvation via the prior walkthrough's method).

## Surprises

- None new this step — the clade confound was surfaced from step 2's raw data and is documented
  in the Results table above; no additional surprises arose while framing the test around it.

## Decisions

**2026-08-13 — Test both sensitivity and ecotype associations for every gene; tag candidates
rather than filtering ecotype hits out.** Co-defined with the researcher, matching her prior
approach to this same confound. Full rationale above.

**2026-08-13 — Fisher's exact (presence/absence) and Mann-Whitney U (copy number, normalized by
strain gene count), both with Benjamini-Hochberg FDR correction.** Co-defined; chosen for
appropriateness to small strain counts (exact test, not large-sample approximations) and to
control false discoveries across several thousand simultaneous gene tests.

**2026-08-13 — No sensitivity positive control; ecotype test used as the method sanity check
instead.** Co-defined; avoids inventing a control gene from memory, which the repo's KG-as-sole-
source rule prohibits.

**2026-08-13 — Record direction alongside every test result, not just significance.** Co-defined
with the researcher: a p-value alone says whether a group differs between strain-groups, not
which way. Direction is computed independently for the sensitivity test and the ecotype test
(presence rate or median normalized copy number, whichever side is higher), so a candidate's full
readout states not just "significant" but "enriched in N-sensitive" (or mixed / HL / LL).

**2026-08-13 — Keep both raw and FDR-adjusted p-values in the output, not just the adjusted
value.** Co-defined with the researcher, so the correction's effect on any specific gene stays
inspectable rather than collapsed into a single pass/fail number.

## Decide-gate checklist

- **Outputs produced:** `scripts/01_build_strain_groups.py`; `data/01_strain_group_membership.csv`.
- **Results presented:** clade-vs-sensitivity cross-tab, strain group membership, sensitivity-vs-
  ecotype cross-tab (all above).
- **QC gate:** all 15 strains' clades pulled fresh from the KG (not assumed from memory) and
  every clade successfully mapped to HL or LL (assertion passed, no unrecognized clade values);
  sensitivity x ecotype cross-tab sums to 15 with no strain double-counted or missing.
- **Decisions made this step:** test-both-groupings-and-tag design; Fisher's exact + Mann-Whitney
  U with FDR correction; no fabricated sensitivity positive control (all above, dated 2026-08-13).
- **Advance rationale:** the hypothesis is operational (exactly what "candidate" means, in KG-
  testable terms), the confound that step 2's data surfaced has a concrete, co-defined handling
  plan, and the control strategy is honest about what is and isn't KG-verifiable. Ready for step
  4 to build the actual test functions.
