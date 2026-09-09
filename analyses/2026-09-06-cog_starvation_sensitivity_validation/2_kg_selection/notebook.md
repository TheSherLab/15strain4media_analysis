# Step 2 — KG entries

## Context

Step 1 locked the question (nitrogen-only after the 2026-09-06 reopen) and
decided to reuse the prior walkthrough analysis's 10-experiment scope. This
step resolves the 41 COGs to expression-queryable genes and checks which
have differential-expression evidence.

## Reopened 2026-09-09 — gene resolution rebuilt around COG membership

**Original approach (2026-09-06):** one researcher-supplied Cyanorak ID per
COG (`00_cog_to_ck_id_mapping.csv`), resolved via `genes_by_homolog_group`.

**Why it was reopened:** two identity checks on that mapping
(`04_paralog_audit.py`, `05_name_vs_cyanorak_crosscheck.py`) showed the
single ID is only a *representative* member of each COG, and sometimes a
poor one — `COG0443` ("DnaK") pointed at `CK_00000457`, a DUF3181
"conserved hypothetical protein" (COG category: function unknown), while
the actual `dnaK1/2/3` chaperones sit in three other Cyanorak groups that
the single-ID route never saw.

**The fix (researcher-directed):** the researcher's comparative-genomics
significance test was run on COG **copy number** — every gene carrying a
COG number in a strain. The in-vivo expression validation must use the
same gene set. Resolution now comes straight from the researcher's own
whole-genome annotation (`WhloeGenome_AllStrains_Concated.xlsx`): for each
COG, in each strain, take every CDS tagged with that COG number, and map
its Cyanorak cluster to a multiomics-KG locus tag.

- `COG0443` now resolves to `dnaK1`, `dnaK2`, `dnaK3` **and** the
  hypothetical, in every strain.
- The gene universe grows from ~95 resolved loci to **427 (COG × strain ×
  locus) rows across the 4 evidence strains** (389 distinct loci).
- A few physical genes carry more than one COG number in the annotation
  (e.g. `zipN` is tagged `COG0188`, `COG0484`, `COG2214` — a multi-domain
  protein). These are kept under each COG, matching how the copy-number
  test counted them. 11 Cyanorak clusters span >1 COG this way.

`00_cog_to_ck_id_mapping.csv` is kept in `data/` marked `SUPERSEDED` so the
step-1 → step-2 trail stays in the repo; nothing downstream reads it.

## What I did

**Source data** (`scripts/00_prepare_source_data.py`):
- `data/00_source_normixed_cogs.csv` — the 41 COGs with Direction tag and
  the researcher's functional ordering (sheet order). `COG0640` (smtB/ArsR)
  had a blank Direction cell; the researcher confirmed (2026-09-09, on the
  ordered screenshot) it is **N**. Final tags: **N = 13, mixed = 28**
  (no "unlabeled"). Cross-checked against `cogs_with_sig_p_values.xlsx`
  (which uses N / "Other"): the only difference is `COG0640` (blank there).
- `data/00_wholegenome_cog_members.csv` — every CDS carrying one of the 41
  COGs, across all 13 annotated study strains: **1,254 gene rows, 165
  Cyanorak clusters, 41/41 COGs**.

**Gene resolution** (`scripts/02_resolve_genes.py`): one batched
`genes_by_homolog_group` call over the 165 clusters × 13 strains, joined
back to the whole-genome gene rows by (cluster, strain). Every whole-genome
gene row mapped to a KG locus (0 `unmapped_in_kg`).

**Evidence check** (`scripts/03_check_evidence.py`): each resolved (COG,
strain, locus) against the 10 in-scope experiments' `differential_
expression_by_gene` rows.

**Identity checks that motivated the reopen** (kept for the record):
- `scripts/04_paralog_audit.py` → `data/05_paralog_audit.csv` — member
  count per supplied Cyanorak ID per strain. Every group was internally
  consistent across strains; 7 genuine in-genome duplications.
- `scripts/05_name_vs_cyanorak_crosscheck.py` →
  `data/06_name_vs_cyanorak_crosscheck.csv` — name-route vs ID-route.
  0 cases where a resolved locus didn't carry its supplied ID; surfaced
  the `COG0443`/DnaK representative-gene problem.

## Results

**Direction tags (41 COGs):** N 13, mixed 28.

**Experiment scope:** all 10 experiments from the walkthrough analysis
confirmed present in the live KG — 5 nitrogen (MED4 ×3, MIT9313 ×2),
5 phosphorus (NATL2A ×2, MIT9312 ×1, MED4 ×1, MIT9313 ×1).

**Gene resolution funnel:**

| Stage | Count |
|---|---|
| COGs in the researcher's list | 41 |
| Whole-genome gene rows carrying one of the 41 COGs (13 strains) | 1,254 |
| ...mapped to a KG locus tag | 1,254 (100%) |
| (COG × strain × locus) rows in the 4 evidence strains | 427 |
| distinct loci in the 4 evidence strains | 389 |

**Genes per COG, 4 evidence strains (largest):**

| COG | dir | MED4 | MIT9313 | NATL2A | MIT9312 |
|---|---|---|---|---|---|
| COG0477 | mixed | 25 | 34 | 28 | 25 |
| COG0697 | mixed | 9 | 11 | 9 | 8 |
| COG0845 | mixed | 3 | 11 | 5 | 4 |
| COG0454 | mixed | 5 | 8 | 5 | 5 |
| COG2214 | mixed | 5 | 7 | 6 | 5 |
| COG0188 | N | 6 | 5 | 5 | 6 |
| COG0443 | mixed | 4 | 5 | 4 | 4 |

15 COGs resolve to 1 gene/strain, 16 to 2–4, 7 to 5–9, 3 to ≥10. The 10
largest COGs hold ~305 of the 427 rows.

**Evidence check:**

| Strain | Distinct loci queried | ...with ≥1 DE row |
|---|---|---|
| MED4 (5 N + P experiments) | 84 | 84 |
| MIT9313 (2 N + P experiments) | 131 | 123 |
| MIT9312 (1 P experiment) | 85 | 4 |
| NATL2A (2 P experiments) | 89 | 5 |

**37 of 41 COGs have ≥1 DE row in ≥1 in-scope experiment** (N 10, mixed 27)
— up from 32 under the single-ID approach. The nitrogen side (MED4,
MIT9313) is well covered; the phosphorus side is sparse at the raw level
(Martiny reclassification in step 5 addresses this — see step 5).

**4 COGs have zero evidence anywhere:** `COG0367` (AsnB), `COG1454` (EutG),
`COG3206` (GumC), `COG3727` (Vsr). Each resolves only in strains whose
in-scope experiment doesn't cover it (e.g. AsnB/GumC/Vsr resolve in MIT9312
only, whose sole P experiment is the 38-protein Fuszard table).

## Surprises

- **The single-ID mapping was systematically narrow, not wrong per se.**
  Every supplied ID is a real member of its COG; it's just one of several.
  The DnaK case is the clearest — the representative was a hypothetical
  while the named chaperones were excluded.
- **11 Cyanorak clusters carry more than one of the 41 COG numbers.** These
  are multi-domain genes (e.g. `zipN` → COG0188/COG0484/COG2214; the
  adenylate-cyclase/CHASE2 fusion → COG2114/COG4252). Kept under each COG,
  consistent with the copy-number test; a small non-independence the
  pooled hit-rate carries.
- **COG0477 is 25–34 genes per strain** (every MFS permease). Included per
  the researcher's decision that the expression set must match the
  copy-number test's set; flagged for step 6 evaluation.

## Decisions

**2026-09-09 — Resolve COGs to their full per-strain gene membership from
the researcher's whole-genome annotation, not to one Cyanorak ID each.**
Co-defined with the researcher: "when I check the copy number and did the
statistical test I used all, this needs to be the same." Supersedes the
2026-09-06 single-ID mapping.

**2026-09-09 — COG0640 Direction = N.** Researcher-confirmed on the ordered
screenshot; the sheet cell was blank.

## Decide-gate checklist

- **Outputs produced:** `scripts/00_prepare_source_data.py` (rewritten),
  `scripts/02_resolve_genes.py` (rewritten), `scripts/03_check_evidence.py`,
  `scripts/04_paralog_audit.py`, `scripts/05_name_vs_cyanorak_crosscheck.py`;
  `data/00_source_normixed_cogs.csv`, `data/00_wholegenome_cog_members.csv`,
  `data/00_cog_to_ck_id_mapping.csv` (superseded),
  `data/02_gene_locus_resolution.csv`, `data/03_gene_evidence.csv`,
  `data/05_paralog_audit.csv`, `data/06_name_vs_cyanorak_crosscheck.csv`.
- **Results presented:** resolution funnel; genes-per-COG table; evidence
  table; zero-evidence COG list (all above).
- **QC gate:** 1,254/1,254 whole-genome gene rows mapped to a KG locus
  (0 unmapped). Direction cross-check vs `cogs_with_sig_p_values.xlsx`:
  1 expected difference (COG0640). Multi-COG clusters (11) inspected and
  confirmed genuine multi-domain genes, kept under each COG. DnaK
  resolution manually confirmed (dnaK1/2/3 + hypothetical, all strains).
- **Decisions made this step:** full-membership resolution (2026-09-09);
  COG0640 = N (2026-09-09).
- **Advance rationale:** the gene universe now matches the copy-number
  test's gene set (every gene per COG per strain), 37 of 41 COGs have
  expression evidence, and identity is triangulated. Ready for step 3 to
  set the framing (background pool, phosphorus handling).
