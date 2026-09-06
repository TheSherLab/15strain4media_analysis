# Step 1 — Research question

## Context

The researcher's own comparative-genomics work (`Dataset 3.xlsx`, sheet
"Significant NorMixed COGs", from the "15 strains 4 media" paper revision)
identified 41 COGs whose cross-strain copy-number/presence pattern
associates significantly with a strain-level trait the researcher labels
"starvation sensitivity" — tagged `N` (12 COGs, nitrogen-linked) or `mixed`
(29 COGs, associated with more than one nutrient/unclear). The researcher
wants to know whether any of these 41 genes have prior in vivo evidence
(RNA-seq/proteomics/microarray differential expression under literal
nutrient starvation) supporting the hypothesis that they matter for
starvation sensitivity — i.e. does existing published expression data give
this comparative-genomics finding any biological support. This reuses the
method (`hit_rate()`, `bootstrap_pvalue()`) built in the prior
`2026-08-10-np_starvation_expression_walkthrough/` analysis, applied to
this new gene list, with the explicit requirement that these 41 genes be
excluded from whatever background/randomness-check population is used
(mirroring how that analysis excluded its own 92-gene list from its
background pool).

## KG context

Grounding queries run before locking scope:

- Source file check (`Dataset 3.xlsx`, sheet "Significant NorMixed COGs"):
  47 rows, 41 with data (6 trailing rows empty). Columns: General
  annotation, COGs (NCBI COG accession, e.g. `COG0188`), p-value,
  Direction (`N`/`mixed`/blank — 12/28/1... see below), p-value (HL LL),
  Protien [sic] (gene/protein name), Annotation, Wider annotation.
  `Direction` value counts: `mixed` 28, `N` 12, blank 1 (`smtB/ArsR`, COG0640).
- `kg_schema()` — checked whether NCBI COG accessions are indexed anywhere
  in this KG build. Result: **no.** `OrthologGroup` nodes come from two
  sources only — `cyanorak` (7,640 groups, `id` format `cyanorak:CK_...`,
  same system the 92-gene list used) and `eggnog` (40,748 groups, `id`
  format `eggnog:<cluster>@<taxon>`, e.g. `eggnog:1H1HZ@1129`). Neither
  uses NCBI COG accessions as an identifier. `CogFunctionalCategory` nodes
  exist but only encode single-letter functional categories (e.g. "J" =
  translation), not individual COG numbers. No gene's `all_identifiers`
  list contains an NCBI COG accession (checked `COG0188` directly — zero
  hits).
- Name-based resolution feasibility, tested on 12 of the 41 protein names
  two ways:
  - `resolve_gene(name)` (whole-KG, all 47 organisms): unreliable for this
    list. `TonB` → 1 hit total, in *Pseudomonas putida* (not
    Prochlorococcus at all). `CDC9` (yeast naming for DNA ligase I) → 0
    hits (bacterial name is `ligA`). Several other names return hits
    concentrated in *Alteromonas*/other genera rather than
    Prochlorococcus.
  - `search_homolog_groups(name, source="cyanorak")`: better-scoped
    (Cyanorak groups cover ~21 cyanobacterial/related organisms).
    `ligA` → 2 candidate groups; `DnaK` → 3 paralogs (`dnaK1/2/3`, all
    "chaperone protein DnaK"); `GyrA` → 2 candidate groups; `TonB` → 3
    candidates, none a clean single match (`exbB`/`exbD`/an
    uncharacterized "TonB box-containing" group — TonB-system components,
    not necessarily TonB itself); `Wza` → 0 matches.

Structural implication (superseded — see Decisions below): initial testing
suggested this list would need per-gene ortholog-group resolution built
from scratch in step 2 via Cyanorak name/product search, since it doesn't
ship with Cyanorak IDs the way the 92-gene list did. Systematic testing
across all 41 protein names found name search alone insufficient: only 10
of 41 names return exactly one Cyanorak group, 12 return multiple (one,
`DnaJ`, returns an implausible 12 candidates — a search-breadth artifact,
not real paralogy), and 17 of 41 return zero matches even after trying a
Prochlorococcus-scoped gene-name fallback and a function-text search (some
zero-match cases, e.g. `PPE`-repeat proteins, are plausibly genuinely
absent from a marine cyanobacterium; others may just use different KG
vocabulary). This was resolved by the discovery below — the researcher
located her own COG-to-Cyanorak-ID mapping, making the tiered fuzzy-search
plan unnecessary.

## COG-to-Cyanorak-ID mapping (supersedes name-based resolution)

The researcher provided a pre-existing 41-row mapping (`cog_numbers` ->
`CK_ID`) built during her own comparative-genomics work. Verified before
adopting it:

- **Complete coverage:** all 41 COGs from `Dataset 3.xlsx` appear in the
  mapping; the two lists match exactly (no COG missing either direction).
- **One shared CK_ID, confirmed legitimate, not a data-entry error:**
  `COG4252` (CHASE2, "Extracytoplasmic sensor domain") and `COG2114`
  (AcyC, "Adenylate cyclase") both map to `cyanorak:CK_00002384`. Checked
  this group directly in the KG: its members' `gene_summary` field reads
  "adenylate cyclase :: CHASE2" — a single fused multi-domain protein
  (adenylate cyclase with a CHASE2 sensor domain), correctly carrying two
  separate NCBI COG numbers for its two domains. Both COGs resolve to the
  same one gene per strain, as expected for a domain fusion, not an error
  to fix.
- **Resolution against the 15 study strains:** `genes_by_homolog_group`
  with all 40 unique CK_IDs (41 COGs, 1 shared pair), restricted to the 15
  strains — **364 gene rows, zero unmatched groups, zero unresolved
  strains.** Every one of the 40 groups has >=1 member gene in >=1 of the
  15 strains (per-strain counts: 22-29 genes). This is a dramatic
  improvement over name-based search's 41%-zero-match rate.

Given this, step 2 uses the supplied CK_ID mapping directly (the same
mechanical join the 92-gene list used), not the tiered
name-search-with-disambiguation approach originally planned.

## What I did

Reviewed the source Excel file and confirmed which sheet/columns the
researcher meant. Ran the KG grounding queries above before proposing
scope, then worked through two clarifying questions with the researcher:

1. **What "Direction" (`N`/`mixed`) means** — researcher confirmed: `N` =
   this COG's cross-strain pattern links specifically to nitrogen: `mixed`
   = links to both nutrients or is unclear. This mirrors the matched/cross
   structure of the prior analysis's 92-gene list, so results will be
   reported for all 41 pooled and for the two subsets separately.
2. **Which nutrient(s) "starvation sensitivity" refers to** — researcher
   chose to test both: run the method against both the nitrogen- and
   phosphorus-starvation experiment sets already validated in the prior
   walkthrough (same axenic/uninfected, literal-starvation, 15-study-strain
   scope — no reason to redo that scoping work, since these are the same
   organism/treatment-type constraints).
3. **Paralog handling** — when a protein name matches more than one
   Cyanorak ortholog group (DnaK, GyrA, TonB, etc.), researcher chose to
   keep all matching paralogs as separate rows/entities, consistent with
   Rule 2 (locus tags, not gene names) and how the prior analysis handled
   `urtA`/`urtE`'s shared-ID collision.

## Decisions

**2026-09-06 — Reuse the prior walkthrough's experiment scope (10
experiments, axenic/uninfected, literal starvation, 15-study-strain
restriction) rather than re-deriving it.** The organism/treatment-type
constraints are identical (same 15 strains, same "literal starvation vs.
replete" definition, same axenic/uninfected requirement) — re-running step
1/2's scoping logic from scratch would reproduce the same 10 experiments
already validated. Step 2 of this analysis will re-verify this rather than
assume it silently.

**2026-09-06 — No matched/cross split by nutrient annotation for this
list.** The 92-gene list's genes were each annotated as N-specific or
P-specific acquisition genes, enabling a matched-vs-cross test. This list's
genes (DNA repair, chaperones, membrane transport, LPS biosynthesis) are
general cellular-maintenance/stress machinery, not nutrient-acquisition
genes — the `N`/`mixed` tags describe the *statistical* association
pattern from the comparative-genomics work, not a biological
nitrogen-specific vs. phosphorus-specific function. So both nutrients are
tested as separate direct hit-rate checks (does this gene set respond to N
starvation? to P starvation?), not as a matched/cross pair.

**2026-09-06 — Background/randomness check excludes both gene lists, not
just the 41 target genes.** The researcher pointed out that the prior
analysis's 92-gene N/P-acquisition list is also, by definition, a set of
genes already hypothesized to respond to starvation — leaving them in the
"random" background pool would let known responders inflate the null
baseline rate, working against detecting a real signal in the 41-COG set.
So the nitrogen background pool (built in step 3) will exclude: the 41
target COGs' resolved loci, the 92-gene list's resolved loci (reusing
`2026-08-10-.../2_kg_selection/data/02_gene_locus_resolution.csv`), and
the nitrogen-keyword product/name matches used in the prior analysis. As
in the prior analysis, this check is only statistically valid on the
nitrogen side (the only nutrient with an unfiltered "every gene tested"
background population in scope) — to be confirmed, not assumed, in step 3.

**2026-09-06 — Gene resolution uses the researcher's own COG-to-Cyanorak-ID
mapping, not name-based search.** See "COG-to-Cyanorak-ID mapping" above —
verified complete (41/41 COGs covered) and clean (364 gene rows, zero
unresolved groups/strains across the 15 study strains) before adopting it.

## Locked research question

Do the 41 COGs in "Significant NorMixed COGs" (`Dataset 3.xlsx`) — genes
whose cross-strain copy-number/presence pattern the researcher's
comparative-genomics work associates with starvation sensitivity, tagged
`N` (12 genes) or `mixed` (29 genes) — show an in vivo
differential-expression response under literal nitrogen or phosphorus
starvation vs. replete medium, in axenic, uninfected Prochlorococcus
restricted to the researcher's 15 study strains? Reported for all 41
pooled and for the `N`/`mixed` subsets separately. Evidence sources:
RNA-seq, proteomics, and microarray differential expression in the KG,
using the same 10 in-scope experiments validated in
`2026-08-10-np_starvation_expression_walkthrough/`. A second question asks
whether any observed nitrogen-side response is statistically
distinguishable from a background/noise rate (bootstrap + Fisher's exact),
with these 41 genes excluded from that background population.

## Decide-gate checklist

- **Outputs produced:** this notebook; analysis scaffold (`paper.md`,
  `gaps_and_friction.md`, `.gitignore`).
- **Results presented:** source-file structure (47 rows/41 populated,
  Direction value counts); KG schema check for COG-accession indexing
  (none found); name-resolution feasibility comparison (whole-KG
  `resolve_gene` vs. Cyanorak-scoped `search_homolog_groups`, both found
  insufficient); the researcher-supplied COG-to-CK_ID mapping, verified
  complete (41/41) and clean (364 gene rows, 0 unmatched groups/strains
  across the 15 study strains) against the live KG; locked question.
- **QC gate:** confirmed no NCBI-COG-accession node/property exists in this
  KG build (`kg_schema()`, direct `all_identifiers` check on `COG0188`).
  Confirmed the researcher's supplied CK_ID mapping covers all 41 COGs
  exactly (no gaps either direction) and resolves fully against the 15
  study strains. Investigated the one shared CK_ID (`COG4252`/`COG2114` ->
  `CK_00002384`) and confirmed it's a genuine domain-fusion gene
  (adenylate cyclase + CHASE2 sensor), not a mapping error.
- **Decisions made this step:** reuse prior experiment scope rather than
  re-derive; no matched/cross split (test both nutrients directly); keep
  all paralogs as separate rows; background exclusion covers both the
  41-COG list and the prior analysis's 92-gene list; gene resolution uses
  the supplied CK_ID mapping, not name-based search (all dated 2026-09-06,
  co-defined with researcher).
- **Advance rationale:** question and scope are locked, including every
  genuine judgment call this list forced (Direction-column meaning, which
  nutrient(s) to test, paralog handling, background-pool scope) and the
  gene-resolution path is now a clean mechanical join like the original
  analysis's — ready to move to step 2.
