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
(COG, locus, experiment). Timepoints from the walkthrough analysis's final
table, including **Lin at 59h** and the two **Martiny** microarrays at
48h / 24h. `RECLASSIFY_ABSENT_AS_NOT_SIG` — a COG-gene present in the
strain but absent from the source table is counted `not_significant`, not
"no data" — is applied to the two Martiny genome-wide P microarrays
(q<0.05 filter) **and, from 2026-09-09, to the Lin uninfected P RNA-seq**:
the KG marks it `table_scope = significant_only` and RNA-seq measures the
whole transcriptome, so a COG gene not in Lin's reported set was tested
and not significant — same logic as Martiny (see Decisions). Still not
applied to Fuszard (iTRAQ fixed detected-protein list).

**Hit rates** (`scripts/02_compute_hit_rates.py`): `hit_rate()` for
nitrogen and phosphorus, each pooled and split by Direction (N 13 / mixed
28). `bootstrap_pvalue()` (10,000 iterations, seed 42) and Fisher's exact
against the step-3 background — **nitrogen only** (pooled + both subsets).

**Figures** (`scripts/03_figures.py`):
- `figures/01_gene_experiment_heatmap.png` — experiments as rows (5 N
  block, 5 P block), the 41 COGs as columns in the researcher's functional
  order, **bracketed and labelled by her "General annotation" category**
  (Quality control – DNA level / Protein level, Amino acid & mixotrophy
  biosynthesis, LPS, Membrane linker, Exopolysaccharide,
  Biofilm/Attachment, Energy production, Translation/transcription/signal
  transduction, RTX toxins, mixed), with vertical dividers between groups.
  All COG/protein labels are black; a per-COG **N / mixed Direction strip**
  sits directly under each category bracket. Coloured cell text is
  **`k / n / m`** (always all three) — `k` significant genes / `n` genes
  measured in that experiment / `m` copies of the COG in that strain's
  genome. `n == m` means the experiment covered every copy (the genome-wide
  Martiny arrays); `n << m` means it didn't (the pre-filtered phosphorus
  tables often list only 1–2 of a COG's genes, so `prop` under Lin reads
  `1/2/28`). Colour = dominant direction, intensity scaled by `k/n`.
  Walkthrough palette; COGs absent from a strain's genome drawn hatched.
- `figures/01_gene_experiment_heatmap.pptx` (`scripts/04_heatmap_pptx.py`,
  `uv run --with python-pptx`) — the same heatmap as an editable
  PowerPoint: 410 named cell rectangles + per-cell `k/n/m` text boxes +
  group brackets, so every element can be recoloured / moved / relabelled.
- `figures/02_pct_significant_vs_background.png` — background line vs. the
  6 hit-rate bars.

**Figure refinement (2026-09-09, later — researcher-requested).** Three
rounds. Rounds 1–2 are presentation only; round 3 changes the phosphorus
numbers (Lin reclassification).

*Round 1 — within-cell direction splits.* The `k / n / m` cell text
collapses a COG's significant genes into one number `k`, and the cell is
coloured by whichever direction has more genes. That hides cells where the
COG's own paralogs move opposite ways. Added
`scripts/05_within_cell_direction_splits.py` →
`data/03_within_cell_direction_splits.csv`, which enumerates every such
cell. **11 cells** qualify — all nitrogen, all MED4, 7 of them in the
Weissberg RNA-seq row.

*Round 2 — split cells drawn orange + a general style pass.* In those 11
cells the fill is now **orange** (not the dominant-direction colour) and
`k` is written `Nup Mdn` on its own line (`COG0477` under Weissberg
RNA-seq reads `3up 6dn / 25/25`); an "Both up- and downregulated genes"
legend entry was added. Same style pass applied to this figure and to
Figure 1 in the sibling `2026-08-10-np_starvation_expression_walkthrough/`
analysis: Arial throughout, no bold, larger fonts; the "tested, no gene
significant" grey darkened (`#c3c2b7` → `#a8a69b`); experiment rows
relabelled to the Oxford-style citation `First author et al., YEAR -
analysis type - strain` in black (the `N:` / `P:` prefixes replaced by a
left-side `Nitrogen starvation` / `Phosphorus starvation` bracket).

*Round 3 — legibility + Lin reclassification.* Fonts enlarged again;
orange split-cell `k` written with `↑n↓n` arrows instead of `Nup Mdn`;
"absent from strain genome" cells drawn with **three** diagonal strokes on
a near-white fill and every cell given a visible grey border, so the three
low-signal states (tested-not-significant grey / no-data beige / absent
white+strokes) are separable — the "no data" colour was also darkened to
`#e7e4d8`. **And** the Lin uninfected P RNA-seq joined
`RECLASSIFY_ABSENT_AS_NOT_SIG` (see Decisions) — this is not cosmetic: the
phosphorus tested count rose 246 → 344, phosphorus hit rate 5.3% → 3.8%
(the 13 significant genes unchanged). Nitrogen is untouched. Same Lin
change applied in the sibling walkthrough analysis.

## Results

**Hit rates** (`data/02_hit_rate_results.csv`):

| Group | Distinct COGs | Distinct loci | Tests | Up | Down | % significant |
|---|---|---|---|---|---|---|
| Nitrogen — all 41 pooled | 37 | ~230 | 438 | 37 | 31 | **15.5%** |
| Nitrogen — Direction N (13) | 10 | 55 | 99 | 12 | 11 | 23.2% |
| Nitrogen — Direction mixed (28) | 27 | ~175 | 339 | 25 | 20 | 13.3% |
| Phosphorus — all 41 pooled | 38 | 308 | 344 | 8 | 5 | **3.8%** |
| Phosphorus — Direction N (13) | 10 | 53 | 53 | 2 | 0 | 3.8% |
| Phosphorus — Direction mixed (28) | 28 | 259 | 291 | 6 | 5 | 3.8% |

**Nitrogen background comparison** (`data/02_h_summary.txt`; background
overall rate **19.1%**, 1,274/6,664):

| Group | Observed | Bootstrap null mean (max) | Bootstrap p | Fisher OR | Fisher p |
|---|---|---|---|---|---|
| Nitrogen pooled (41) | 15.5% | 18.3% (25.2%) | 0.930 | 0.78 | 0.067 |
| Nitrogen Direction N (13) | 23.2% | 22.9% (39.8%) | 0.443 | 1.28 | 0.30 |
| Nitrogen Direction mixed (28) | 13.3% | 16.8% (27.4%) | 0.969 | 0.65 | **0.0067** |

**Phosphorus** is descriptive only (no valid background): **3.8%**
significant (13/344), 8 up / 5 down — essentially flat. (Was 5.3% / 246
tests before the 2026-09-09 Lin reclassification, which added ~98 Lin
NATL2A tested-not-significant calls; the 13 significant genes are the
same.)

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
- **Phosphorus really is flat.** With the Martiny + Lin reclassifications
  the phosphorus side now covers 38 of 41 COGs (308 loci, 344 tests), and
  only 13 are significant, split 8 up / 5 down (3.8%). The candidate list
  shows no phosphorus-starvation response.

- **11 nitrogen cells have paralogs moving opposite ways within one COG.**
  `COG0477` (prop/MFS permeases) is the clearest — 3 up / 6 down under
  Weissberg RNA-seq, 2 up / 3 down under Read, 1 up / 4 down under Tolonen.
  `COG0697` (RhaT) is 4 up / 3 down under Weissberg. `COG0188` (GyrA)
  3 up / 2 down; `COG0399` (WecE) 1/1 in two experiments; `COG0443`
  (DnaK), `COG1596` (Wza), `COG1793` (ligA) each 1/1 under Weissberg
  RNA-seq. These are broad transporter/repair families, not coherent
  single-direction responses — consistent with the overall negative
  result. Full list: `data/03_within_cell_direction_splits.csv`.

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

**2026-09-09 — Lin uninfected P RNA-seq added to
`RECLASSIFY_ABSENT_AS_NOT_SIG`.** Researcher-directed, overriding the
walkthrough analysis's 2026-09-08 "not applied to Lin" note. Grounds: the
live KG marks the Lin uninfected experiment
(`..._plimited_natl2a_rnaseq_uninfected`) `table_scope = significant_only`,
and RNA-seq measures the whole transcriptome — so a COG gene with a NATL2A
locus but no row in Lin's reported set was measured and not significant,
the same reasoning already applied to the two Martiny microarrays. Effect:
phosphorus tested count 246 → 344, phosphorus hit rate 5.3% → 3.8%
(13 significant genes unchanged, 8 up / 5 down). Nitrogen unaffected (Lin
is phosphorus). Cascades to `03_figures.py`, `05_within_cell_direction_
splits.py` (no new split cells — all 11 are nitrogen), and step 6's
sensitivity checks (nitrogen-only — unchanged).

## Decide-gate checklist

- **Outputs produced:** `scripts/01_extract_target_de.py` (phosphorus +
  Lin 59h + Martiny **and Lin** table-absent reclassification),
  `scripts/02_compute_hit_rates.py`
  (nutrient split, N/mixed subsets, N-subset background test),
  `scripts/03_figures.py` (flipped, fraction-based, functional-category
  brackets, `k/n/m` cell text, orange split cells, Arial/no-bold style
  pass, citation row labels + N/P bracket, darker grey, single-line
  absent), `scripts/04_heatmap_pptx.py` (editable PowerPoint, same),
  `scripts/05_within_cell_direction_splits.py` (2026-09-09 refinement);
  `data/01_target_gene_experiment_matrix.csv` (1,179 rows),
  `data/02_hit_rate_results.csv`, `data/02_bootstrap_null_distribution.csv`,
  `data/02_h_summary.txt`, `data/03_within_cell_direction_splits.csv`
  (11 rows); `figures/01_gene_experiment_heatmap.png` + `.pptx`,
  `figures/02_pct_significant_vs_background.png`.
- **Results presented:** hit-rate table, background comparison table,
  phosphorus descriptive result, both figures (all above).
- **QC gate:** matrix row count still 1,179 (Lin reclassification moves
  ~98 phosphorus rows from `no_data_at_timepoint` to `not_significant`, no
  rows added/removed). `status × nutrient` crosstab checked: nitrogen
  counts identical to the pre-Lin run; phosphorus `not_significant`
  331 (was ~233), `no_data_at_timepoint` 183 (was ~281). Bootstrap reused
  seed 42 / 10,000 iterations; nitrogen bootstrap/Fisher numbers identical
  to the pre-Lin run. Split-cell rendering cross-checked: the 11 orange
  cells in `03_figures.py` match
  `data/03_within_cell_direction_splits.csv` gene-for-gene (GyrA/Weissberg
  `↑3↓2`, COG0477/Weissberg `↑3↓6`, etc.). Figures re-rendered and visually
  checked (Arial found; brackets, citation labels, 3-stroke absent cells,
  distinct no-data beige all render).
- **Decisions made this step:** redo on reopened inputs; heatmap
  fraction-based; N-subset background test folded in; split-direction cells
  drawn orange with `↑n↓n` arrows; Arial/no-bold style pass + citation row
  labels + N/P bracket; Lin uninfected P RNA-seq added to the table-absent
  reclassification (all 2026-09-09).
- **Advance rationale:** both nutrients computed on the full COG gene set,
  every hypothesis has concrete numbers and a figure, the negative pooled
  result is stable across the resolution change. Ready for step 6.
