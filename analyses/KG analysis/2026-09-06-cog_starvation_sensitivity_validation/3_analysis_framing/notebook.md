# Step 3 — Analysis framing

## Context

Step 2 (reopened 2026-09-09) rebuilt gene resolution: every gene carrying
each of the 41 COG numbers in each strain, from the researcher's
whole-genome annotation — 427 (COG × strain × locus) rows in the 4
evidence strains, 37 of 41 COGs with expression evidence. Step 1 (reopened
2026-09-09) restored phosphorus, reported descriptively. This step
operationalizes the hypothesis and rebuilds the nitrogen background pool
against the new (larger) gene set.

## What I did

**Significance criteria** (`scripts/01_copy_significance_criteria.py`):
reused the walkthrough analysis's empirically-derived
significance-criteria table for these same 10 experiments — no new
derivation (same platforms, same thresholds).

**Background pool** (`scripts/02_build_background_pool.py`): every gene
tested in the 4 unfiltered ("all_detected_genes") **nitrogen** experiments
at the walkthrough's chosen starvation timepoint per experiment,
excluding:
1. this analysis's full 41-COG gene set, by (organism, locus) — now
   ~1,141 (organism, locus) pairs across 13 strains (was ~360 under the
   single-ID approach);
2. the walkthrough's 92-gene N/P-acquisition list, by (organism, locus),
   read live from its `2_kg_selection/data/02_gene_locus_resolution.csv`
   (so it picks up that analysis's 2026-09-08 gene-identity corrections);
3. the same 13 nitrogen-keyword product/name text matches.

## Results

**Background pool** (`data/02_negative_background_pool.csv`, 6,664 rows):

| Experiment | Distinct background genes |
|---|---|
| Weissberg proteomics (MED4) | 1,320 |
| Weissberg RNA-seq (MED4) | 1,711 |
| Tolonen microarray (MED4) | 1,560 |
| Tolonen microarray (MIT9313) | 2,073 |

Overall background significant rate: **19.1%** (1,274 / 6,664).

**Framing (co-defined with the researcher):**

- **Hypothesis:** the 41-COG candidate list — using the full per-strain COG
  membership, matching the copy-number test's gene set — shows a
  nitrogen-starvation response elevated above the background/random rate.
  Reported for all 41 pooled and for the `N` (13) / `mixed` (28) Direction
  subsets. Phosphorus reported the same way but descriptively.
- **No positive controls** (unchanged): this list is candidate genes, not
  textbook markers — nothing to sanity-check against.
- **Nitrogen background only** (unchanged): the phosphorus tables are
  pre-filtered / detection-limited even after the Martiny reclassification
  (Lin 34 genes, Fuszard 38/63 proteins), so there is no unbiased
  phosphorus population to resample. The Martiny reclassification makes the
  phosphorus *hit rate* computable but not a *background* — a genome-wide
  microarray filtered to q<0.05 has no "tested-and-listed non-significant"
  genes to draw from.
- **Expected outcome, operationally:** `hit_rate()` on the 41-COG set (and
  subsets) vs. the 5 nitrogen experiments, compared via `bootstrap_pvalue()`
  and Fisher's exact to the background pool; phosphorus hit rate reported
  without a test. Same two functions from the walkthrough's
  `4_methods/np_response.py`, reused (step 4 verified, not reimplemented).

## Surprises

- The exclusion set grew ~3× (the full COG membership is ~1,141 loci vs.
  ~360 for one ID per COG). Background pool shrank only modestly
  (6,914 → 6,664) — most COG genes fall outside the 4 background
  experiments' tested-gene lists anyway.

## Decisions

**2026-09-09 — Background pool rebuilt against the full COG gene set;
phosphorus stays background-free.** Follows mechanically from the step-1
and step-2 reopens; the no-phosphorus-background rationale is unchanged
from 2026-09-06, only reconfirmed now that phosphorus is back in scope
descriptively.

## Decide-gate checklist

- **Outputs produced:** `scripts/01_copy_significance_criteria.py`,
  `scripts/02_build_background_pool.py`; `data/01_significance_criteria.csv`,
  `data/02_negative_background_pool.csv` (6,664 rows).
- **Results presented:** background pool per experiment; overall background
  rate; framing (above).
- **QC gate:** significance-criteria table has the expected rows.
  Exclusion counts per experiment nonzero and proportionate. 92-gene list
  read from the walkthrough's live (2026-09-08-corrected) resolution file.
- **Decisions made this step:** background rebuilt against the full COG
  gene set; phosphorus background-free (2026-09-09).
- **Advance rationale:** hypothesis is operational for both nutrients, the
  nitrogen background is rebuilt against the correct gene set, and the
  phosphorus-descriptive-only choice is reconfirmed. Ready for step 5.
