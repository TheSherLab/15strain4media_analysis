# Step 2 — KG entries

## Context

Step 1 locked the question: do the 54 N/P-acquisition genes with in vivo
evidence show a significant expression change under nitrogen/phosphorus
starvation in axenic Prochlorococcus, matched and cross? This step
identifies exactly which experiments, publications, and gene loci that
question will be evaluated against.

## Co-define

Proposed to the researcher: enumerate axenic Prochlorococcus N/P
starvation experiments, then re-resolve the 54 genes to locus tags per
relevant strain (fresh from the KG, not reused from the prior triage
analysis's cached strain-name-only results). Agreed before starting.

A scoping question came up while enumerating nitrogen experiments (see
Surprises) and was resolved with the researcher before finalizing the
experiment list.

**Redo 1 (2026-08-02):** during step 3 co-define, the researcher
clarified that SS120 (CCMP1375) is not part of their study — removed
from the experiment list and gene-locus resolution. Prior version (14
experiments / 6 strains, including SS120) preserved in git history at
commit `2da4c3b`.

**Redo 2 (2026-08-02):** while building step 3's background gene set,
discovered MIT9301's 2 phosphorus experiments are METABOLOMICS
(metabolite abundance) and carry zero gene-level DE data —
`differential_expression_by_gene` returns nothing for them. Restricted
the experiment list to gene-expression/proteomics omics types only
(RNASEQ, PROTEOMICS, MICROARRAY), which drops MIT9301 entirely (it had
no nitrogen experiments either). Prior version (12 experiments / 5
strains, including MIT9301's metabolomics rows) preserved in git history
at commit `2dc2e86`.

## What I did

- `01_export_evidence_gene_list.py` — froze the 54-gene evidence list
  (from the prior triage analysis's
  `n_p_genes_experimental_evidence_final.csv`, `has_expression_evidence
  == True`) as this analysis's own input.
- `02_list_np_experiments.py` — `list_experiments(treatment_type=[...],
  organism='Prochlorococcus', background_factors=['axenic'],
  verbose=True)` for nitrogen and phosphorus, excluding: the 4
  alternate-N-source experiments, all SS120/CCMP1375 experiments, and
  all METABOLOMICS-omics-type experiments (see Surprises and Decisions).
- `03_resolve_gene_loci.py` — resolved each of the 54 genes to locus tags
  in the 4 strains that appear in the final experiment list (MED4,
  MIT9312, MIT9313, NATL2A), via `resolve_gene` by name for 53 genes, the
  KEGG-orthology route (K01077, phoB-subtracted, phosphatase-product-
  filtered) for `phoA`, and a direct zero-padded `gene_overview` lookup
  for the 4 genes whose "name" is actually an unpadded MED4 locus tag
  (`PMM707`→`PMM0707`, etc. — a known issue from the prior analysis).

## Results

**10 axenic Prochlorococcus N/P starvation experiments** (5 nitrogen, 5
phosphorus; gene-expression/proteomics data only), across 6
publications, spanning 4 strains:

| n_or_p | organism_name | omics_type | treatment | control | publication_doi | is_time_course |
|---|---|---|---|---|---|---|
| N | MED4 | PROTEOMICS | PRO99-lowN nutrient starvation | PRO99-lowN exponential growth | 10.1101/2025.11.24.690089 | True |
| N | MED4 | RNASEQ | PRO99-lowN nutrient starvation | PRO99-lowN exponential growth | 10.1101/2025.11.24.690089 | False |
| N | MED4 | RNASEQ | N-depleted Pro99 medium | N-replete Pro99 medium | 10.1038/ismej.2017.88 | True |
| N | MED4 | MICROARRAY | Nitrogen deprivation (MED4) | N-replete Pro99 medium (MED4) | 10.1038/msb4100087 | True |
| N | MIT9313 | MICROARRAY | Nitrogen deprivation (MIT9313) | N-replete Pro99 medium (MIT9313) | 10.1038/msb4100087 | True |
| P | NATL2A | RNASEQ | P-limited | P-replete control | 10.1111/1462-2920.13104 | True |
| P | MIT9312 | PROTEOMICS | Phosphate deplete (10 uM NaH2PO4) | Phosphate replete (50 uM NaH2PO4) | 10.1186/2046-9063-8-7 | False |
| P | NATL2A | PROTEOMICS | Phosphate deplete (10 uM NaH2PO4) | Phosphate replete (50 uM NaH2PO4) | 10.1186/2046-9063-8-7 | False |
| P | MED4 | MICROARRAY | Phosphate starvation | P-replete Pro99 medium | 10.1073/pnas.0601301103 | True |
| P | MIT9313 | MICROARRAY | Phosphate starvation | P-replete Pro99 medium | 10.1073/pnas.0601301103 | True |

Full detail in `data/02_np_experiments.csv`.

**53/54 genes resolved to at least one locus tag** among the 4 relevant
strains (`ptxB/phnD2` has none — its only prior evidence was in MIT9301,
now out of scope; see step 3's notebook). Loci per (nutrient side,
strain):

| n_or_p | strain | distinct loci |
|---|---|---|
| N | MED4 | 22 |
| N | MIT9312 | 22 |
| N | MIT9313 | 24 |
| N | NATL2A | 25 |
| P | MED4 | 25 |
| P | MIT9312 | 19 |
| P | MIT9313 | 20 |
| P | NATL2A | 21 |

(A gene can resolve to more than one locus per strain when paralogs
exist, e.g. `ureA`; a strain's count can also be below the 29 N / 25 P
totals when a gene has no ortholog annotated in that strain.) Full
detail in `data/03_gene_loci.csv`.

## Surprises

- **Nitrogen "treatment_type" is not homogeneous.** Of the Prochlorococcus
  nitrogen experiments, some are coculture (out of scope per step 1) and,
  of the axenic ones, 4 are "growth on alternate N source
  (cyanate/urea/nitrite) vs N-replete" rather than starvation/deprivation
  — a different biological question. This mattered because several
  nitrogen-annotated genes in the 54-gene list are themselves
  alternate-N-source-utilization genes (`cynA/B/D/S`, `ureA–G`, `nirA`).
  Raised with the researcher; decided to exclude them and keep nitrogen
  strictly to starvation/deprivation contrasts.
- Phosphorus has one non-axenic experiment (NATL2A, phage-infected,
  `10.1111/1462-2920.13104`) alongside its axenic sibling in the same
  publication — excluded per the axenic-only scope from step 1.
- `PMM707/PMM719/PMM721/PMM722` are unpadded MED4 locus tags, not gene
  symbols — same issue the prior triage analysis found; fixed the same
  way (zero-pad, direct `gene_overview` lookup).
- **METABOLOMICS experiments carry no gene-level DE data.** MIT9301's
  two phosphorus experiments are metabolomics assays (metabolite
  abundance) — `differential_expression_by_gene` returns zero rows for
  them, since there is no Gene-level DE edge for a metabolomics
  experiment in this KG. Discovered in step 3 while building the
  background gene set; not visible from `list_experiments` alone (it
  reports `omics_type` but the practical consequence — "this experiment
  cannot answer a gene-expression question" — only became clear when the
  DE query itself returned empty). Restricting to RNASEQ/PROTEOMICS/
  MICROARRAY drops MIT9301 entirely.

## Decisions

- **2026-08-02** — Nitrogen scope narrowed to strict starvation/
  deprivation contrasts only; "growth on alternate N source" experiments
  excluded.
- **2026-08-02 (redo 1)** — SS120 (CCMP1375) excluded entirely at the
  researcher's request (not part of their study).
- **2026-08-02 (redo 2)** — METABOLOMICS-omics-type experiments excluded
  (no gene-level DE data); MIT9301 drops out of the analysis entirely as
  a consequence, since both its phosphorus experiments were metabolomics
  and it had no nitrogen experiments.

## Decide-gate checklist

- **Outputs produced** —
  `scripts/01_export_evidence_gene_list.py` → `data/01_evidence_gene_list.csv`
  (54 rows);
  `scripts/02_list_np_experiments.py` → `data/02_np_experiments.csv`
  (10 rows);
  `scripts/03_resolve_gene_loci.py` → `data/03_gene_loci.csv`
  (178 rows). All run via `uv run python analyses/2026-08-02-np_gene_expression_full_analysis/2_kg_selection/scripts/<script>.py`
  from the repo root (must run from repo root, not the script's own
  directory, for `.env` KG-credential discovery to work).
- **Results presented** — experiment table and gene-resolution counts
  shown inline above, matching what was shown to the researcher in chat.
- **QC gate** — confirmed 53/54 genes resolve to >=1 locus tag in the 4
  relevant strains (`ptxB/phnD2` correctly has none, given MIT9301's
  exclusion — not a resolution failure); confirmed all 10 remaining
  experiments are RNASEQ/PROTEOMICS/MICROARRAY (0 metabolomics, 0 SS120).
- **Decisions made this step** — see Decisions above.
- **Advance rationale** — the KG entries this analysis will run against
  are enumerated, restricted to omics types that actually carry
  gene-level DE data, and gene→locus mapping is complete for the
  53 usable genes; step 3 (analysis framing) is already underway against
  this corrected KG-entries set.
