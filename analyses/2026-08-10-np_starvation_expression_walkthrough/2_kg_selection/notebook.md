# Step 2 — KG entries

## Context

Step 1 locked the question and four scope decisions: axenic/uninfected
only, literal starvation/limitation only (no alternate-N-source
substrate-switch experiments), restricted to the researcher's 15 study
strains. This step applies those filters mechanically to the KG's
Prochlorococcus nitrogen/phosphorus experiments, and separately resolves
the researcher's 92-gene list to locus tags in the strains that survive.

## What I did

**Experiments.** `scripts/01_select_experiments.py` pulled all
Prochlorococcus experiments with `treatment_type=["nitrogen"]` (12) and
`treatment_type=["phosphorus"]` (9) via `list_experiments(verbose=True)`,
then applied the step-1 scope as a per-experiment classifier: strain must
be one of the 15 study strains; omics type must not be metabolomics-only
(no gene-level DE); background factors must not include `coculture` or
`viral`; treatment/experiment-name text must not indicate an
alternate-N-source substrate switch or a chemical-inhibitor proxy
(`"as sole n source"`, `"azaserine"`).

**Genes.** `scripts/02_resolve_genes.py` read the researcher's source gene
list (`data/00_source_gene_list.csv`, copied from `Dataset 2.xlsx`, 92
rows, 91 carrying a Cyanorak ortholog-group ID). For each gene with a
Cyanorak ID, `genes_by_homolog_group(group_ids=["cyanorak:<CK_id>"],
organisms=[the 4 strains below])` returns the locus tag in each strain
where a member exists — this is the ortholog-group route, not name-based
lookup, so it doesn't depend on gene-name spelling and reports each strain
separately rather than picking one match (Rule 2: locus tags, not gene
names — paralogs and strain-specific presence/absence are both real
biology here, not noise to collapse away).

## Results

**Experiment filter funnel:**

| Nutrient | KG total (Prochlorococcus, treatment-typed) | In scope after step-1 filters |
|---|---|---|
| Nitrogen | 12 | 5 |
| Phosphorus | 9 | 5 |

**Included experiments (10 total, 4 strains, 6 publications):**

| Nutrient | Strain | Omics | Publication |
|---|---|---|---|
| N | MED4 | PROTEOMICS | Weissberg et al. 2025 (bioRxiv, 10.1101/2025.11.24.690089) |
| N | MED4 | RNASEQ | Weissberg et al. 2025 (bioRxiv, 10.1101/2025.11.24.690089) |
| N | MED4 | RNASEQ | Read et al. 2017 (ISME J, 10.1038/ismej.2017.88) |
| N | MED4 | MICROARRAY | Tolonen et al. 2006 (Mol Syst Biol, 10.1038/msb4100087) |
| N | MIT9313 | MICROARRAY | Tolonen et al. 2006 (Mol Syst Biol, 10.1038/msb4100087) |
| P | NATL2A | RNASEQ | Lin, Ding, Zeng 2015 (Environ Microbiol, 10.1111/1462-2920.13104) |
| P | MIT9312 | PROTEOMICS | Fuszard et al. 2012 (Aquat Biosyst, 10.1186/2046-9063-8-7) |
| P | NATL2A | PROTEOMICS | Fuszard et al. 2012 (Aquat Biosyst, 10.1186/2046-9063-8-7) |
| P | MED4 | MICROARRAY | Martiny, Coleman, Chisholm 2006 (PNAS, 10.1073/pnas.0601301103) |
| P | MIT9313 | MICROARRAY | Martiny, Coleman, Chisholm 2006 (PNAS, 10.1073/pnas.0601301103) |

Full table with exclusion reasons: `data/01_np_experiments.csv`.

**Publication funnel** (`scripts/04_publication_table_and_funnel.py`, using
`list_publications` counts gathered interactively): 43 publications total in
this KG build -> 33 involving Prochlorococcus -> 8 with a nitrogen- or
phosphorus-treatment experiment -> 6 with >=1 experiment surviving the
step-1 scope. Figure: `figures/01_publication_funnel.png`.

**The 6 publications used, plus the 2 that had N/P-treatment experiments
but were fully excluded** (full table: `data/04_publications_table.csv`):

| Nutrient | Publication | What it tested | Data type | Strains | Status |
|---|---|---|---|---|---|
| N | Weissberg et al. 2025, bioRxiv | MED4 nutrient starvation vs. exponential growth, axenic + coculture (90-day time course) | RNA-seq, Proteomics | MED4 | Included (axenic arm) |
| N | Read et al. 2017, ISME J | MED4 N-depleted vs. N-replete Pro99, TSS mapping | RNA-seq | MED4 | Included |
| N | Tolonen et al. 2006, Mol Syst Biol | MED4+MIT9313 N-deprivation time course vs. replete (+ alt-N-source arms, excluded) | Microarray | MED4, MIT9313 | Included (deprivation arms) |
| P | Lin, Ding, Zeng 2015, Environ Microbiol | NATL2A P-limited vs. P-replete, +/- phage infection | RNA-seq | NATL2A | Included (uninfected arm) |
| P | Fuszard et al. 2012, Aquat Biosyst | MIT9312, NATL2A, SS120 Pi-replete vs. Pi-deplete (iTRAQ) | Proteomics | MIT9312, NATL2A | Included (SS120 arm excluded) |
| P | Martiny, Coleman, Chisholm 2006, PNAS | MED4+MIT9313 phosphate starvation vs. replete, time course | Microarray | MED4, MIT9313 | Included |
| N | Dominguez-Martin et al. 2017, mSystems | SS120 azaserine (chemical N-limitation proxy) vs. untreated | Proteomics | SS120 | **Excluded** (strain + proxy design) |
| P | Kujawinski et al. 2023, mSystems | MIT9301 metabolite profiling, replete vs. P-limited | Metabolomics | MIT9301 | **Excluded** (no gene-level DE) |

**Excluded experiments (11 total):** 2 coculture-background (Weissberg
MED4, nitrogen); 1 strain not in the 15-strain list (SS120 azaserine,
nitrogen); 4 alternate-N-source substrate-switch (Tolonen cyanate/urea x2/
nitrite); 2 metabolomics-only (Kujawinski, MIT9301, phosphorus); 1
phage-infection background (Lin et al., NATL2A infected arm, phosphorus);
1 strain not in the 15-strain list (SS120, phosphorus).

**Gene resolution funnel:**

| Stage | Count |
|---|---|
| Genes in researcher's source list | 92 |
| ...with a Cyanorak ortholog-group ID | 91 |
| ...resolved to >=1 locus tag in >=1 of the 4 in-scope strains | 61 |

**Per-strain resolved gene counts** (of the 92-gene list, how many have a
locus tag in that strain's genome):

| Strain | Resolved genes |
|---|---|
| MED4 | 53 |
| NATL2A | 48 |
| MIT9312 | 44 |
| MIT9313 | 42 |

Full gene x strain resolution table: `data/02_gene_locus_resolution.csv`.

**Evidence check** (`scripts/03_check_evidence.py`): a locus tag only means
the gene exists in a strain's genome. Checking each resolved (gene,
strain) pair against the 10 in-scope experiments'
`differential_expression_by_gene` rows narrows this to genes that were
actually measured:

| Strain | Resolved genes | ...with >=1 DE row in an in-scope experiment |
|---|---|---|
| MED4 | 53 | 52 |
| MIT9313 | 42 | 41 |
| NATL2A | 48 | 14 |
| MIT9312 | 44 | 2 |

**60 of the 92 genes (29 N-annotated, 31 P-annotated) have actual
differential-expression evidence in at least one in-scope experiment.**
Full table: `data/03_gene_evidence.csv`.

**Decision-tree artifact:** a running decision log (updated after every
step's approval) is published at
https://claude.ai/code/artifact/031b1561-15a1-4655-a05c-2ee0db6fe7d6 —
covers steps 1-2 so far, including the experiment-inclusion rule as a
flowchart and the publication table above.

## Surprises

- **31 of 92 genes have no locus tag in any of the 4 in-scope strains.**
  For 90 of these genes this reflects genuine strain-level gene
  presence/absence (well documented in Prochlorococcus ecotypes — e.g.
  phosphorus-acquisition gene content varies substantially between HL and
  LL ecotypes), not a resolution failure: their Cyanorak ortholog group
  exists in the KG but has no member gene in these 4 particular strains.
- **`phnW` has no Cyanorak ID in the source spreadsheet at all**, so it
  fell through ortholog-based resolution. Direct name lookup
  (`resolve_gene("phnW")`) finds it only in *Pseudomonas putida* KT2440 —
  confirms it is not present in the KG for any Prochlorococcus strain, not
  a lookup gap on our end.
- Strain resolution counts track roughly with each strain's total gene
  count in the KG (MED4 1,976 genes / 53 resolved; MIT9313 2,948 genes /
  42 resolved) — MIT9313 is genomically the largest of the 4 strains but
  resolves fewer target genes, consistent with it being a low-light
  ecotype with different nutrient-acquisition gene content than the
  high-light strains (MED4, MIT9312) the gene list may be more weighted
  toward. `[interpretation]` — not verified against ecotype/clade
  annotation in this step; a candidate check for step 3 if it matters for
  interpreting cross-strain patterns.
- **MIT9312's evidence rate collapses from 44 resolved genes to only 2
  with actual DE evidence.** Cause: MIT9312's only in-scope experiment
  (Fuszard et al. 2012 proteomics) has `table_scope = significant_only`
  ("only proteins reaching the paper's fold-change cutoff (>1.6 or <0.6)
  are reported") — the source table itself only contains 38 genes total,
  8 of them significant. This is a data-scarcity property of the
  publication, not a resolution problem; it means MIT9312 will contribute
  almost nothing to this analysis regardless of what step 3 decides. Flag
  for `gaps_and_friction.md`.

## Decide-gate checklist

- **Outputs produced:** `scripts/01_select_experiments.py`,
  `scripts/02_resolve_genes.py`, `scripts/03_check_evidence.py`,
  `scripts/04_publication_table_and_funnel.py`;
  `data/00_source_gene_list.csv` (copied input),
  `data/01_np_experiments.csv`, `data/02_gene_locus_resolution.csv`,
  `data/03_gene_evidence.csv`, `data/04_publications_table.csv`;
  `figures/01_publication_funnel.png`.
- **Results presented:** experiment filter funnel and included-experiment
  table; gene resolution funnel and per-strain table; gene evidence
  funnel and per-strain table (all above).
- **QC gate:** all 91 Cyanorak IDs in the source list matched a known
  ortholog group in the KG (zero `group_not_found_in_kg`) → resolution
  failures are strain-absence, not ID mismatches. `phnW`'s absence
  cross-checked by direct name lookup → confirmed absent, not a tooling
  gap. MIT9312's low evidence rate traced to its one in-scope experiment's
  `table_scope = significant_only` → data-scarcity property of the source
  publication, not a bug.
- **Decisions made this step:** none (mechanical application of step-1
  scope; no new judgment calls forced by the data).
- **Advance rationale:** the usable experiment set (10 experiments, 4
  strains, 6 publications) and the usable gene universe (60 of 92 genes
  with actual DE evidence in >=1 in-scope experiment: 29 N-annotated, 31
  P-annotated) are both established and QC'd; ready for step 3 to set the
  hypotheses' operational controls (positive/negative controls, what
  "responds" means numerically) against this 60-gene universe.
