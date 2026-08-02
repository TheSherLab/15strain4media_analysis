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

## What I did

- `01_export_evidence_gene_list.py` — froze the 54-gene evidence list
  (from the prior triage analysis's
  `n_p_genes_experimental_evidence_final.csv`, `has_expression_evidence
  == True`) as this analysis's own input.
- `02_list_np_experiments.py` — `list_experiments(treatment_type=[...],
  organism='Prochlorococcus', background_factors=['axenic'],
  verbose=True)` for nitrogen and phosphorus, with the 4 alternate-N-source
  experiments excluded per the researcher's decision (see Surprises).
- `03_resolve_gene_loci.py` — resolved each of the 54 genes to locus tags
  in the 6 strains that appear in the experiment list (MED4, MIT9301,
  MIT9312, MIT9313, NATL2A, SS120/CCMP1375), via `resolve_gene` by name
  for 53 genes, the KEGG-orthology route (K01077, phoB-subtracted,
  phosphatase-product-filtered) for `phoA`, and a direct zero-padded
  `gene_overview` lookup for the 4 genes whose "name" is actually an
  unpadded MED4 locus tag (`PMM707`→`PMM0707`, etc. — a known issue from
  the prior analysis).

## Results

**14 axenic Prochlorococcus N/P starvation experiments** (6 nitrogen, 8
phosphorus), across 8 publications:

| n_or_p | organism_name | omics_type | treatment | control | publication_doi | is_time_course |
|---|---|---|---|---|---|---|
| N | MED4 | PROTEOMICS | PRO99-lowN nutrient starvation | PRO99-lowN exponential growth | 10.1101/2025.11.24.690089 | True |
| N | MED4 | RNASEQ | PRO99-lowN nutrient starvation | PRO99-lowN exponential growth | 10.1101/2025.11.24.690089 | False |
| N | MED4 | RNASEQ | N-depleted Pro99 medium | N-replete Pro99 medium | 10.1038/ismej.2017.88 | True |
| N | SS120 (CCMP1375) | PROTEOMICS | Azaserine 10 uM (ferredoxin-GOGAT inhibitor, proxy for N-limitation) | Untreated control | 10.1128/mSystems.00008-17 | False |
| N | MED4 | MICROARRAY | Nitrogen deprivation (MED4) | N-replete Pro99 medium (MED4) | 10.1038/msb4100087 | True |
| N | MIT9313 | MICROARRAY | Nitrogen deprivation (MIT9313) | N-replete Pro99 medium (MIT9313) | 10.1038/msb4100087 | True |
| P | MIT9301 | METABOLOMICS | replete + P-limited | replete (Pro99) | 10.1128/msystems.01261-22 | False |
| P | MIT9301 | METABOLOMICS | replete + P-limited | replete (Pro99) | 10.1128/msystems.01261-22 | False |
| P | NATL2A | RNASEQ | P-limited | P-replete control | 10.1111/1462-2920.13104 | True |
| P | MIT9312 | PROTEOMICS | Phosphate deplete (10 uM NaH2PO4) | Phosphate replete (50 uM NaH2PO4) | 10.1186/2046-9063-8-7 | False |
| P | NATL2A | PROTEOMICS | Phosphate deplete (10 uM NaH2PO4) | Phosphate replete (50 uM NaH2PO4) | 10.1186/2046-9063-8-7 | False |
| P | SS120 (CCMP1375) | PROTEOMICS | Phosphate deplete (10 uM NaH2PO4) | Phosphate replete (50 uM NaH2PO4) | 10.1186/2046-9063-8-7 | False |
| P | MED4 | MICROARRAY | Phosphate starvation | P-replete Pro99 medium | 10.1073/pnas.0601301103 | True |
| P | MIT9313 | MICROARRAY | Phosphate starvation | P-replete Pro99 medium | 10.1073/pnas.0601301103 | True |

Full detail in `data/02_np_experiments.csv`.

**54/54 genes resolved to at least one locus tag** among the 6 relevant
strains (29 nitrogen-annotated, 25 phosphorus-annotated — corrects the
"28/26" approximation used in step 1's dialogue). Loci per (nutrient
side, strain):

| n_or_p | strain | distinct loci |
|---|---|---|
| N | MED4 | 22 |
| N | MIT9301 | 19 |
| N | MIT9312 | 22 |
| N | MIT9313 | 24 |
| N | NATL2A | 25 |
| N | SS120 | 10 |
| P | MED4 | 25 |
| P | MIT9301 | 23 |
| P | MIT9312 | 19 |
| P | MIT9313 | 20 |
| P | NATL2A | 21 |
| P | SS120 | 17 |

(A gene can resolve to more than one locus per strain when paralogs
exist, e.g. `ureA`; a strain's count can also be below 29/25 when a gene
has no ortholog annotated in that strain.) Full detail in
`data/03_gene_loci.csv`.

## Surprises

- **Nitrogen "treatment_type" is not homogeneous.** Of the 12
  Prochlorococcus nitrogen experiments, 2 are coculture (already out of
  scope per step 1) and, of the remaining 10 axenic ones, 4 are "growth
  on alternate N source (cyanate/urea/nitrite) vs N-replete" rather than
  starvation/deprivation — a different biological question (can the
  strain use compound X as its N source, not "what happens when N runs
  out"). This mattered here because several nitrogen-annotated genes in
  the 54-gene list are themselves alternate-N-source-utilization genes
  (`cynA/B/D/S`, `ureA–G`, `nirA`), so those 4 experiments looked like
  the most on-target evidence for exactly those genes. Raised with the
  researcher; decided to exclude them and keep nitrogen strictly to the
  6 starvation/deprivation contrasts, for symmetry with the phosphorus
  side (which has no equivalent complication — all 8 axenic phosphorus
  experiments are limitation/deplete-vs-replete).
- Phosphorus has one non-axenic experiment (NATL2A, phage-infected,
  `10.1111/1462-2920.13104`) alongside its axenic sibling in the same
  publication — excluded per the axenic-only scope from step 1.
- `PMM707/PMM719/PMM721/PMM722` are unpadded MED4 locus tags, not gene
  symbols — same issue the prior triage analysis found; fixed the same
  way (zero-pad, direct `gene_overview` lookup).

## Decisions

- **2026-08-02** — Nitrogen scope narrowed to strict starvation/
  deprivation contrasts only; the 4 "growth on alternate N source"
  experiments are excluded from this analysis (see Surprises). This
  means genes whose only relevant nitrogen evidence was an
  alternate-N-source experiment (e.g. `cynA/B/D`, `ureA–G`, `nirA`) may
  show fewer or no significant nitrogen-side hits here than in the prior
  triage analysis, which did not distinguish starvation from
  alternate-source experiments.

## Decide-gate checklist

- **Outputs produced** —
  `scripts/01_export_evidence_gene_list.py` → `data/01_evidence_gene_list.csv`
  (54 rows);
  `scripts/02_list_np_experiments.py` → `data/02_np_experiments.csv`
  (14 rows);
  `scripts/03_resolve_gene_loci.py` → `data/03_gene_loci.csv`
  (250 rows). All run via `uv run python analyses/2026-08-02-np_gene_expression_full_analysis/2_kg_selection/scripts/<script>.py`
  from the repo root (must run from repo root, not the script's own
  directory, for `.env` KG-credential discovery to work).
- **Results presented** — experiment table and gene-resolution counts
  shown inline above, matching what was shown to the researcher in chat.
- **QC gate** — confirmed 54/54 genes resolve to >=1 locus tag in the 6
  relevant strains (0 unresolved, after fixing the PMM7xx zero-padding
  case); confirmed experiment counts (6 N / 8 P) match the researcher-
  approved scope (axenic-only, strict-starvation nitrogen).
- **Decisions made this step** — see Decisions above.
- **Advance rationale** — the KG entries this analysis will run against
  are enumerated and gene→locus mapping is complete for all 54 genes;
  ready to move to step 3 (analysis framing: hypothesis, controls,
  expected outcome).
