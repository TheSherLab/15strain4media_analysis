# Step 5 — Analyze

## Context

Redo (2026-09-09) after the step-1 (phosphorus restored), step-2 (full COG
membership) and step-3 (background rebuilt) reopens. This step runs the
hit-rate + background comparison on the new 41-COG gene set across all 10
experiments and produces the figures.

## What I did

**Extraction** (`scripts/01_extract_target_de.py`): for every gene
carrying one of the 41 COGs (step-2 resolution), pulled DE data from each
of the 10 experiments at its chosen starvation timepoint. One row per
(COG, locus, experiment). Timepoints reused verbatim from the walkthrough
analysis's final (2026-09-08) table, including **Lin at 59h** and the two
**Martiny** microarrays at 48h / 24h. The two Martiny genome-wide P
microarrays get `RECLASSIFY_ABSENT_AS_NOT_SIG` — a COG-gene present in the
strain but absent from Martiny's q<0.05 table is counted `not_significant`,
not "no data". Not applied to Lin (curated operon table) or Fuszard (iTRAQ
fixed detected-protein list).

**Hit rates** (`scripts/02_compute_hit_rates.py`): `hit_rate()` for
nitrogen and phosphorus, each pooled and split by Direction (N 13 / mixed
28). `bootstrap_pvalue()` (10,000 iterations, seed 42) and Fisher's exact
against the step-3 background — **nitrogen only** (pooled + both subsets).

**Figures** (`scripts/03_figures.py`):
- `figures/01_gene_experiment_heatmap.png` — experiments as rows (5 N
  block, 5 P block), the 41 COGs as columns in the researcher's functional
  order, label `[N]`/`[m]` + colour = Direction. A COG has several genes
  per strain, so each cell shows the **share of the COG's tested genes
  that are significant** (`k/n` printed when >1 gene), colour = dominant
  direction, intensity scaled by the fraction. Walkthrough palette; COGs
  absent from a strain's genome drawn hatched.
- `figures/02_pct_significant_vs_background.png` — background line vs. the
  6 hit-rate bars.

## Results

**Hit rates** (`data/02_hit_rate_results.csv`):

| Group | Distinct COGs | Distinct loci | Tests | Up | Down | % significant |
|---|---|---|---|---|---|---|
| Nitrogen — all 41 pooled | 37 | ~230 | 438 | 37 | 31 | **15.5%** |
| Nitrogen — Direction N (13) | 10 | 55 | 99 | 12 | 11 | 23.2% |
| Nitrogen — Direction mixed (28) | 27 | ~175 | 339 | 25 | 20 | 13.3% |
| Phosphorus — all 41 pooled | ~20 | ~150 | 246 | 8 | 5 | **5.3%** |
| Phosphorus — Direction N (13) | — | — | 40 | 2 | 0 | 5.0% |
| Phosphorus — Direction mixed (28) | — | — | 206 | 6 | 5 | 5.3% |

**Nitrogen background comparison** (`data/02_h_summary.txt`; background
overall rate **19.1%**, 1,274/6,664):

| Group | Observed | Bootstrap null mean (max) | Bootstrap p | Fisher OR | Fisher p |
|---|---|---|---|---|---|
| Nitrogen pooled (41) | 15.5% | 18.3% (25.2%) | 0.930 | 0.78 | 0.067 |
| Nitrogen Direction N (13) | 23.2% | 22.9% (39.8%) | 0.443 | 1.28 | 0.30 |
| Nitrogen Direction mixed (28) | 13.3% | 16.8% (27.4%) | 0.969 | 0.65 | **0.0067** |

**Phosphorus** is descriptive only (no valid background): **5.3%**
significant (13/246), 8 up / 5 down — essentially flat.

## Surprises

- **The pooled nitrogen hit rate (15.5%) is *below* background (19.1%)** —
  the same direction as the 2026-09-06 result (16.2% vs 19.0%), now on the
  full COG gene set. Neither the pooled test nor either subset clears
  background: the `N` subset trends above (23.2%) but the bootstrap says
  that is well within the null (p=0.44), and the `mixed` subset is
  significantly *below* background (Fisher p=0.0067) — these
  general-maintenance genes (chaperones, transporters, DNA repair) respond
  to nitrogen starvation *less* often than random genes.
- **The COG-level heatmap needed a fraction, not a "most significant"
  collapse.** A first pass coloured each cell by the single
  most-significant gene in the COG; with COGs of 25–34 genes (COG0477
  especially) that made almost every big-COG cell look significant. The
  figure now shows `k/n` significant and scales colour by the fraction —
  COG0477's Weissberg-RNA-seq cell is `9/25`, not "red".
- **Phosphorus really is flat.** Even with the Martiny reclassification
  adding ~230 tested data points, only 13 are significant, split up/down.
  The candidate list shows no phosphorus-starvation response.

## Decisions

**2026-09-09 — Redo the whole step on the reopened inputs.** Not a new
judgment call — step 5 consumes steps 2 and 3. Phosphorus now included
(step-1 reopen); gene set is the full COG membership (step-2 reopen);
Lin 59h + Martiny reclassification carried from the walkthrough.

**2026-09-09 — Heatmap cell = fraction of the COG's genes significant, not
the single most-significant gene.** Co-defined implicitly by the
"use all genes" decision — a per-COG summary that lets a 30-gene COG and a
1-gene COG sit in the same figure honestly. Colours and legend still match
the walkthrough (up red / down blue / not-sig grey / no-data pale /
absent hatched).

**2026-09-09 — `N` subset tested against background directly.** The
2026-09-06 notebook flagged this as a follow-up; it is now in
`02_compute_hit_rates.py` (bootstrap p=0.44 — not distinguishable). The
separate `04_stability_check_n_subset.py` is removed as redundant.

## Decide-gate checklist

- **Outputs produced:** `scripts/01_extract_target_de.py` (phosphorus +
  Lin 59h + Martiny reclassification), `scripts/02_compute_hit_rates.py`
  (nutrient split, N/mixed subsets, N-subset background test),
  `scripts/03_figures.py` (flipped, fraction-based, walkthrough palette);
  `data/01_target_gene_experiment_matrix.csv` (1,179 rows),
  `data/02_hit_rate_results.csv`, `data/02_bootstrap_null_distribution.csv`,
  `data/02_h_summary.txt`; `figures/01_gene_experiment_heatmap.png`,
  `figures/02_pct_significant_vs_background.png`.
- **Results presented:** hit-rate table, background comparison table,
  phosphorus descriptive result, both figures (all above).
- **QC gate:** matrix row count reconciled (1,179 = per-experiment
  (COG,locus) pairs + no-locus rows). COG1212 / COG3839 P-block cells
  spot-checked against the raw matrix after the first heatmap pass showed
  a big-COG rendering artifact (fixed — see Surprises). Bootstrap reused
  seed 42 / 10,000 iterations.
- **Decisions made this step:** redo on reopened inputs; heatmap
  fraction-based; N-subset background test folded in (all 2026-09-09).
- **Advance rationale:** both nutrients computed on the full COG gene set,
  every hypothesis has concrete numbers and a figure, the negative pooled
  result is stable across the resolution change. Ready for step 6.
