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
organisms=[the 4 strains below])` returns the group's member gene(s) in
each strain — the ortholog-group route, not name-based lookup, so it
doesn't depend on gene-name spelling and reports each strain separately
rather than picking one match (Rule 2: locus tags, not gene names —
paralogs and strain-specific presence/absence are both real biology here,
not noise to collapse away). **A Cyanorak ID with exactly one member per
strain resolves to that locus; multi-member groups and wrong/incomplete
IDs in the source spreadsheet are handled by the 2026-09-08 reopen — see
that section below.** All tables in this notebook show the post-reopen
resolution.

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
locus tag in that strain's genome; post-reopen):

| Strain | Resolved genes |
|---|---|
| MED4 | 53 |
| NATL2A | 48 |
| MIT9312 | 44 |
| MIT9313 | 42 |

Full gene x strain resolution table: `data/02_gene_locus_resolution.csv`
(now carries `n_group_members_in_strain` / `all_members_in_strain`
columns). One (gene x strain) pair — `PMM719` in NATL2A — is
`ambiguous_multi_member` and left unresolved (2 candidate loci, neither
with any in-scope DE data, so immaterial).

**Evidence check** (`scripts/03_check_evidence.py`): a locus tag only means
the gene exists in a strain's genome. Checking each resolved (gene,
strain) pair against the 10 in-scope experiments'
`differential_expression_by_gene` rows narrows this to genes that were
actually measured:

| Strain | Resolved genes | ...with >=1 DE row in an in-scope experiment |
|---|---|---|
| MED4 | 53 | 52 |
| MIT9313 | 42 | 42 |
| NATL2A | 48 | 15 |
| MIT9312 | 44 | 3 |

**60 of the 92 genes (29 N-annotated, 31 P-annotated) have actual
differential-expression evidence in at least one in-scope experiment.**
Full table: `data/03_gene_evidence.csv`. (Post-reopen: `phoE` and `unkP2`
gained evidence via their corrected loci; `ptrA` gained a NATL2A locus.
The 60-gene / 29-N / 31-P totals are unchanged — the fixes moved which
locus each gene points at, not the gene count.)

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
- **MIT9312's evidence rate collapses from 44 resolved genes to only 3
  with actual DE evidence** (2 before the `phoE` fix — `phoE` = PMT9312_0721
  is one of the 38 proteins Fuszard detected). Cause: MIT9312's only
  in-scope experiment (Fuszard et al. 2012 proteomics) has
  `table_scope = significant_only`
  ("only proteins reaching the paper's fold-change cutoff (>1.6 or <0.6)
  are reported") — the source table itself only contains 38 genes total,
  8 of them significant. This is a data-scarcity property of the
  publication, not a resolution problem; it means MIT9312 will contribute
  almost nothing to this analysis regardless of what step 3 decides. Flag
  for `gaps_and_friction.md`.

## Reopened 2026-09-08 — gene identity (the Cyanorak ID is not always 1:1)

**Original lock (2026-08-10):** every gene resolved through its source-list
Cyanorak ID via `genes_by_homolog_group`; the QC gate checked only that
each ID matched *a* known ortholog group. Decisions: "none — mechanical."

**The data reveal:** while tracing a Figure 1 cell during step 6, the
researcher found `phoE` reading as "not significant" in Lin et al. when
the KG clearly has a strongly P-induced porin (log2FC 6.46, padj 2e-4) in
the NATL2A pho island. `phoE`'s Cyanorak ID (`CK_00002330`) turned out to
be a **3–6-member "outer membrane porin" family per strain**, not one
gene, and `02_resolve_genes.py` was silently keeping the *last* member
(built a dict keyed by `organism_name`). The pick was a non-pho-island
porin in MED4, MIT9312 and NATL2A.

**Two checks added to find every similar case:**

1. **`scripts/05_paralog_audit.py`** — for each of the 91 Cyanorak IDs,
   count `genes_by_homolog_group` members per in-scope strain. >1 member =
   the ID is a family and the "pick one" logic was guessing. Found **4**:
   `phoE` (`CK_00002330`), `unkP2` (`CK_00003432`), `pstS` in MIT9313
   (`CK_00043821`), `PMM719` in NATL2A (`CK_00044628`). Output:
   `data/05_paralog_audit.csv`.
2. **`scripts/06_name_vs_cyanorak_crosscheck.py`** — an independent second
   identifier. Resolve each gene by *name* (and normalised ORF forms,
   `PMM707` → `PMM0707`) with `resolve_gene`, compare to the locus the
   Cyanorak route produced. `DISAGREE` (both resolve, to different loci)
   and `name_only` (name finds a locus the ID missed) are the red flags.
   Found **2**: `urtA` (DISAGREE, all 4 strains) and `ptrA` (`name_only`,
   NATL2A). Output: `data/06_name_vs_cyanorak_crosscheck.csv`. After the
   fixes below, re-run is clean (0 DISAGREE, 0 name_only).

**Each hit confirmed by hand** — `gene_overview`/`gene_details` (what the
locus *is*), `gene_homologs` (which group it *really* belongs to),
`gene_neighbors` (genomic synteny — the pho regulon is a physical
cluster: `phoB`-`phoR`-`phoA`-porin-`pstS`), and
`differential_expression_by_gene` (does the chosen locus behave like the
gene should — the corrected `phoE` locus jumps ~100-fold at Martiny 48h).

**The 6 corrections** (`MANUAL_LOCUS_OVERRIDE` / `CYANORAK_ID_OVERRIDE` /
`FORCE_LOCUS` in `02_resolve_genes.py`, each with its one-line basis):

| gene | was | now | why |
|---|---|---|---|
| `phoE` MED4 | PMM1121 | **PMM0709** | multi-member family; PMM0709 sits between `phoA` (PMM0708) & `pstS` (PMM0710) — log2FC 117 at Martiny 48h |
| `phoE` MIT9312 | PMT9312_1515 | **PMT9312_0721** | in the `phoA`(0720)/`pstS`(0722) cluster; rank-1 sig_up in Fuszard |
| `phoE` MIT9313 | PMT_2631 | **PMT0998** | next to `phoB` (PMT0994)/`phoR`/`pstS` |
| `phoE` NATL2A | PMN2A_1757 | **PMN2A_0440** | between `phoA` (PMN2A_0439) & `pstS` (PMN2A_0441) |
| `unkP2` MED4 | PMM2011 | **PMM0715** | multi-member; PMM0715 is in the pho island (`arsR`/`arsB` neighbours), sig_up Martiny 48h |
| `urtA` (4 strains) | PMM0974 &c. (= `urtE`) | **PMM0970 &c.** | source spreadsheet gave `urtA` the ID `CK_00008074`, which is `urtE`'s group; both rows carried it, so the analysis had `urtE` twice and no `urtA`. Real `urtA` group is `CK_00000076` |
| `ptrA` NATL2A | (not matched) | **PMN2A_0435** | spreadsheet ID `CK_00056804` is a MED4-only singleton group; NATL2A's `ptrA` (in the pho island, in Lin's gene set) is in `CK_00001606`. MIT9312/MIT9313 genuinely lack the pho-island Crp regulator |

`pstS`-MIT9313 was multi-member but the resolver already had the right one
(`PMT0993`, next to `phoB`); `PMM719`-NATL2A is left `ambiguous` (no
in-scope data).

**Downstream impact** (steps consumed step 2's resolution — cascade
approved by the researcher): step 5 re-extracted and re-ran; step 6
re-evaluated. Net: nitrogen matched hit rate 41.3% → 44.2% (real `urtA` is
more N-responsive than `urtE`); H3 bootstrap p<0.0001 throughout (Fisher
OR 3.01 → 3.39); phosphorus matched rate moved (see step 5). `gaps_and_friction.md`
(2026-09-08) carries the methodology takeaways.

**Re-locked 2026-09-08.**

## Decide-gate checklist

- **Outputs produced (original + reopen):**
  `scripts/01_select_experiments.py`, `scripts/02_resolve_genes.py`
  (rewritten: multi-member handling + override maps),
  `scripts/03_check_evidence.py`,
  `scripts/04_publication_table_and_funnel.py`,
  `scripts/05_paralog_audit.py`, `scripts/06_name_vs_cyanorak_crosscheck.py`;
  `data/00_source_gene_list.csv` (copied input),
  `data/01_np_experiments.csv`, `data/02_gene_locus_resolution.csv`,
  `data/03_gene_evidence.csv`, `data/04_publications_table.csv`,
  `data/05_paralog_audit.csv`, `data/06_name_vs_cyanorak_crosscheck.csv`;
  `figures/01_publication_funnel.png`.
- **Results presented:** experiment filter funnel and included-experiment
  table; gene resolution funnel and per-strain table; gene evidence
  funnel and per-strain table; the 6-correction table (all above).
- **QC gate:** all 91 Cyanorak IDs matched a known ortholog group (zero
  `group_not_found_in_kg`). `phnW`'s absence cross-checked by direct name
  lookup → confirmed absent. MIT9312's low evidence rate traced to its one
  in-scope experiment's `table_scope = significant_only`. **Paralog audit:
  4 multi-member IDs found, all resolved or flagged. Name-vs-Cyanorak
  cross-check: 2 mis-IDs found and fixed; re-run clean.**
- **Decisions made this step:** (reopen, 2026-09-08) — the 6 gene-identity
  corrections above, each with a synteny/expression-confirmed basis;
  `PMM719`-NATL2A left ambiguous. Co-defined with the researcher, who also
  approved cascading the fix through steps 5–6.
- **Advance rationale:** the usable experiment set (10 experiments, 4
  strains, 6 publications) and the usable gene universe (60 of 92 genes
  with DE evidence: 29 N, 31 P) are established, QC'd, and now
  identity-triangulated by two independent checks; ready for step 3's
  operational controls.
