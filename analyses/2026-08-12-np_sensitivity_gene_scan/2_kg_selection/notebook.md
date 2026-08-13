# Step 2 — KG entries

## Context

Step 1 locked the question: scan every gene across the 15 study strains' genomes for a
presence/copy-number pattern that lines up with the N-sensitive (9 strains) vs. mixed (6 strains)
survival phenotype. Unlike the prior walkthrough, this step isn't about enumerating publications
or experiments — it's about pulling the raw genome-content data (genes + cross-strain ortholog-
group identity) the genome-wide scan in step 5 will run on, and characterizing what the KG
actually offers for that comparison.

## What I did

**Bulk pull** (`scripts/01_pull_gene_ortholog_data.py`): no single higher-level tool covers
"every gene x its ortholog group(s), scoped to a set of organisms" — `genes_by_homolog_group`
takes group IDs as input (wrong direction), `gene_homologs` takes locus tags one at a time. Used
`run_cypher` as the documented escape hatch (`docs://guide/start_here`) for exactly this shape of
question, after checking `kg_schema` for the `Gene_in_ortholog_group` relationship and verifying
both queries on a small scope first. Pulled (a) every gene for the 15 strains and (b) every
`Gene_in_ortholog_group` edge for those genes, across both `OrthologGroup` sources (`cyanorak`,
`eggnog`).

**Consolidation** (`scripts/02_build_gene_group_assignment.py`): built one row per gene with its
final cross-strain group assignment, co-defined assignment rule:
1. Cyanorak group (curated, `specificity_rank=0`) if the gene has one — primary, for consistency
   with the prior 92-gene analysis.
2. Else, Eggnog `specificity_rank=1` (tightest/family-level) group if it has one — fallback, used
   only for genes Cyanorak didn't curate. Rank 1 chosen over the looser rank 2/3 because broader
   ranks average more strains per group (avg. 7.5–8.1 of 15) and so carry less power to
   distinguish the two phenotype groups.
3. Else, "true orphan" — no group in either system; excluded from the per-ortholog-group scan
   (step 5), tracked separately as a per-strain count.

## Results

**Gene coverage funnel** (full data: `data/03_gene_group_assignment.csv`):

| Stage | Count | % of 32,053 |
|---|---|---|
| Total genes (all 15 strains) | 32,053 | 100% |
| ...with a Cyanorak group | 29,356 | 91.6% |
| ...without a Cyanorak group | 2,697 | 8.4% |
| ...of those, covered by Eggnog rank 1 | 2,179 | 6.8% |
| ...of those, no group in either system ("true orphans") | 518 | 1.6% |
| ...of true orphans, unknown function (`gene_category=Unknown`) | 389 | 1.2% |

**Ortholog-group landscape** (full data: `data/04_group_landscape_summary.csv`):

| Source | Rank | Meaning | Distinct groups | Avg. strains/group (of 15) |
|---|---|---|---|---|
| Cyanorak | 0 | curated, tightest | 5,003 | 5.73 |
| Eggnog | 1 | tight (family-level) | 3,374 | 8.00 |
| Eggnog | 2 | medium (order-level) | 3,664 | 7.54 |
| Eggnog | 3 | broadest (domain-level) | 2,986 | 8.05 |

Cyanorak groups averaging only 5.73 of 15 strains means most groups are **not** universal across
the 15 strains — real presence/absence variation exists for the step-5 scan to search through,
not a KG where every gene trivially appears everywhere.

**True orphans, per strain** (full data: `data/05_orphan_summary_by_strain.csv`):

| Strain | Group | True orphans | ...of which unknown function |
|---|---|---|---|
| MIT1314 | N-sensitive | 146 | 73 |
| MIT1327 | mixed | 95 | 76 |
| MIT9313 | mixed | 83 | 72 |
| PAC1 | mixed | 38 | 32 |
| NATL2A | mixed | 34 | 32 |
| MIT9515 | N-sensitive | 23 | 21 |
| MIT9202 | N-sensitive | 21 | 16 |
| NATL1A | N-sensitive | 18 | 18 |
| MIT0604 | N-sensitive | 14 | 10 |
| SB | mixed | 9 | 5 |
| MIT9301 | mixed | 8 | 7 |
| MED4 | N-sensitive | 5 | 5 |
| AS9601 | N-sensitive | 5 | 5 |
| MIT9312 | N-sensitive | 3 | 2 |

## Surprises

- Each gene can belong to **multiple** ortholog-group memberships simultaneously (1 Cyanorak +
  up to 3 Eggnog, at 3 different taxonomic breadths) — not a single group per gene as the step-1
  framing implicitly assumed. This is why the consolidation script needed an explicit priority
  rule rather than a simple join.
- True-orphan counts are very uneven by strain: MIT1314 (146), MIT1327 (95), and MIT9313 (83)
  sit well above the rest (single digits to ~38). Not interpreted yet — this is raw-data
  structure, not a tested result — but worth carrying into step 3 framing since it could reflect
  genome-assembly/annotation completeness differences between strains as easily as biology, and
  that distinction matters before reading anything into it.

## Decisions

**2026-08-13 — Eggnog rank 1 (tightest) used as the Cyanorak-orphan fallback**, not rank 2 or 3.
Co-defined with the researcher after seeing the real distinct-group and avg-strains-per-group
numbers above; broader ranks would have covered more of the 2,697 Cyanorak-orphans but with
groups too broad (avg. ~8 of 15 strains) to carry much discriminating power.

## Decide-gate checklist

- **Outputs produced:** `scripts/01_pull_gene_ortholog_data.py`,
  `scripts/02_build_gene_group_assignment.py`; `data/01_genes_all_strains.csv`,
  `data/02_ortholog_edges_all_strains.csv`, `data/03_gene_group_assignment.csv`,
  `data/04_group_landscape_summary.csv`, `data/05_orphan_summary_by_strain.csv`.
- **Results presented:** gene coverage funnel, ortholog-group landscape, true-orphan-by-strain
  table (all above).
- **QC gate:** per-strain gene counts cross-checked against step 1's `list_organisms` grounding
  → exact match (32,053 total, same per-strain breakdown). Cyanorak group assignment checked for
  uniqueness (0 genes with >1 Cyanorak group, as expected for a curated primary system).
- **Decisions made this step:** Eggnog rank-1 fallback choice (above, dated 2026-08-13).
- **Advance rationale:** the full gene x ortholog-group assignment table is built and QC'd, the
  group landscape confirms real cross-strain variation exists to search through, and the
  true-orphan set is characterized (including the uneven per-strain distribution flagged for
  step 3). Ready for step 3 to turn this raw material into an actual statistical test.
