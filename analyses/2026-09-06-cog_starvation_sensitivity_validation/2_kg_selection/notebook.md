# Step 2 — KG entries

## Context

Step 1 locked the question and decided to reuse the prior walkthrough
analysis's 10-experiment scope (axenic/uninfected, literal starvation,
15-study-strain restriction) rather than re-derive it, and to resolve the
41 COGs via the researcher's own COG-to-Cyanorak-ID mapping rather than
name-based search. This step mechanically applies both: re-verifies the
10 experiments still exist in the live KG, resolves all 41 COGs to locus
tags in the 15 strains, and checks which resolved genes actually have DE
evidence in the 10 in-scope experiments.

## What I did

**Source data** (`scripts/00_prepare_source_data.py`): copied and cleaned
the "Significant NorMixed COGs" sheet (41 rows) and the researcher's
COG-to-Cyanorak-ID mapping (41 rows) into this analysis's `data/`. Sanity
check: the two files cover the identical 41 COGs (no gap either
direction). One `Direction` value is blank (`COG0640`, `smtB/ArsR`) —
relabeled `unlabeled` rather than folded into `mixed`, since it isn't
actually tagged either way; reported as its own small category below.

**Experiment scope re-verification** (`scripts/01_verify_experiment_scope.py`):
pulled the prior analysis's 10 in-scope experiment IDs and re-queried
`list_experiments` live for nitrogen- and phosphorus-treatment
Prochlorococcus experiments — all 10 (5 nitrogen, 5 phosphorus) still
present, same KG build. No re-derivation of the selection filter itself
was needed.

**Gene resolution** (`scripts/02_resolve_genes.py`): ran
`genes_by_homolog_group` with the 40 unique Cyanorak group IDs (41 COGs,
1 shared pair — `COG4252`/`COG2114`, see step 1) against the 15 study
strains. First pass silently dropped a second gene whenever one Cyanorak
group matched 2 genes in the same strain (a bug, not a data property);
fixed to keep every match as a separate row instead of collapsing to one
— same fix pattern as the prior analysis's `urtA`/`urtE` collision.
Confirmed 7 such cases are genuine in-genome gene duplications (e.g. MED4
has two separate carbamoyltransferase-family loci, `PMM1246` and
`PMM1251`, both correctly assigned to the same Cyanorak group), not an
artifact — both kept.

**Evidence check** (`scripts/03_check_evidence.py`): cross-referenced
every resolved (COG, strain, locus) row against the 10 in-scope
experiments' `differential_expression_by_gene` rows.

## Results

**Direction tag breakdown** (41 COGs): `mixed` 28, `N` 12, `unlabeled` 1
(`smtB/ArsR`, COG0640).

**Experiment scope:** all 10 experiments from the prior analysis
(2026-08-10-np_starvation_expression_walkthrough) confirmed still present
in the live KG — 5 nitrogen (MED4 x3, MIT9313 x2), 5 phosphorus (NATL2A
x2, MIT9312 x1, MED4 x1, MIT9313 x1).

**Gene resolution funnel:**

| Stage | Count |
|---|---|
| COGs in the researcher's source list | 41 |
| ...with a researcher-supplied Cyanorak ID | 41 |
| ...resolved to >=1 locus tag in >=1 of the 15 study strains | 41 |

**Per-strain resolved COG counts** (of 41):

| Strain | Resolved COGs |
|---|---|
| NATL1A | 29 |
| MIT1327 | 28 |
| SB | 28 |
| MIT9312 | 27 |
| MIT9313 | 27 |
| PAC1 | 27 |
| MIT0604 | 26 |
| NATL2A | 26 |
| MIT9215 | 25 |
| AS9601 | 24 |
| MIT9301 | 24 |
| MED4 | 24 |
| MIT9202 | 23 |
| MIT9515 | 21 |

Full table: `data/02_gene_locus_resolution.csv` (622 rows — includes the 7
genuine same-strain gene-duplication cases as separate rows, and
`not_matched_in_strain` rows for strain-absence).

**Evidence check** (restricted to the 4 strains the 10 in-scope
experiments actually cover — MED4, MIT9313 nitrogen; NATL2A, MIT9312
phosphorus):

| Strain | Resolved COGs queried | ...with >=1 DE row |
|---|---|---|
| MIT9313 | 27 | 27 |
| MED4 | 26 | 24 |
| MIT9312 | 27 | 1 |
| NATL2A | 26 | 0 |

**32 of 41 COGs have DE evidence in >=1 in-scope experiment — and 32 of
those 32 come from the nitrogen side (MED4/MIT9313); only 1 comes from
phosphorus (MIT9312).** NATL2A resolves 26 of 41 COGs to a locus tag but
has **zero** overlap between those 26 loci and the genes actually reported
in either of its 2 in-scope experiments (Lin RNA-seq: 34 genes total,
`significant_only`; Fuszard proteomics: covers the rest of NATL2A's 92
combined distinct genes across both tables, also `significant_only`) —
consistent with the narrow (1.5-2.8% of genome) phosphorus table coverage
already established in the prior analysis; with tables this narrow, a
26-gene panel unrelated to either paper's own significance list has a low
chance of overlapping by chance.

**By Direction tag, COGs with evidence:** `mixed` 22 of 28, `N` 9 of 12,
`unlabeled` 1 of 1.

**9 COGs have zero DE evidence anywhere in the 10 in-scope experiments:**
`COG0367` (AsnB), `COG1454` (EutG), `COG2114`/`COG4252` (AcyC/CHASE2 — the
shared-gene pair, so both zero together), `COG2931` (unnamed, "RTX
toxin-related domain"), `COG3206` (GumC), `COG3727` (Vsr), `COG3754`
(RgpF), `COG5651` (PPE).

## Surprises

- **The phosphorus side contributes almost nothing evidence-wise**, even
  though NATL2A resolves the second-most COGs to a locus tag (26 of 41,
  more than MED4's 24). This is a direct consequence of how narrow the
  in-scope phosphorus tables are (established in the prior analysis:
  1.5-2.8% of genome), not a resolution problem — 0 of NATL2A's 26
  resolved loci happen to fall inside either phosphorus table's own
  filtered gene list. `[interpretation]` This means step 3/5 will likely
  have very little phosphorus-side signal to report for this gene list,
  regardless of framing choices — worth flagging to the researcher before
  step 3 locks the framing, not discovering it after step 5 runs.
- **`COG2114` and `COG4252` are the same gene** (confirmed in step 1 —
  adenylate cyclase + CHASE2 sensor domain fusion), so their identical
  zero-evidence result is expected, not two independent data points.

## Decide-gate checklist

- **Outputs produced:** `scripts/00_prepare_source_data.py`,
  `scripts/01_verify_experiment_scope.py`, `scripts/02_resolve_genes.py`,
  `scripts/03_check_evidence.py`; `data/00_source_normixed_cogs.csv`,
  `data/00_cog_to_ck_id_mapping.csv`, `data/01_np_experiments.csv`,
  `data/02_gene_locus_resolution.csv`, `data/03_gene_evidence.csv`.
- **Results presented:** Direction breakdown; experiment re-verification;
  gene resolution funnel and per-strain table; evidence table by strain
  and by Direction tag; zero-evidence COG list (all above).
- **QC gate:** source sheet and CK_ID mapping confirmed to cover the
  identical 41 COGs (script assertion). All 10 prior-analysis experiment
  IDs confirmed present in the live KG. Same-strain multi-gene-per-group
  cases (7 found) investigated and confirmed genuine gene duplications,
  not a query bug — fixed to keep both rather than silently dropping one.
  NATL2A's zero-evidence result cross-checked directly (confirmed true
  zero-overlap, not a locus-tag mismatch bug).
- **Decisions made this step:** none new (mechanical application of
  step-1 decisions); the same-strain-duplicate handling was a bug fix, not
  a judgment call, and `unlabeled` (vs. folding into `mixed`) is a
  reporting choice noted above rather than a scope decision.
- **Advance rationale:** all 41 COGs resolve to >=1 locus tag; 32 have
  actual DE evidence, heavily concentrated on the nitrogen side (32 of 32
  vs. 1 of 32 from phosphorus). This evidence imbalance triggered
  reopening step 1 (see `1_question/notebook.md`, "Reopened 2026-09-06")
  — the researcher chose to drop phosphorus from scope entirely rather
  than report a near-empty result. Ready for step 3 to set the framing
  (background-pool construction excluding both this 41-COG list and the
  prior 92-gene list, per step 1) against the nitrogen-only evidence base.
