# Step 5 — Analyze

## Context

Step 4 verified the methods. This step runs them on the full 61-gene,
10-experiment dataset, resolves several visualization judgment calls
raised by the researcher (gene inclusion, timepoint handling, duplicate
measurements), and produces the main supplementary figure.

## What I did

**Timepoint selection.** Time-course experiments cannot be pooled across
timepoints if the goal is one clean "starvation vs. replete" comparison
per experiment (co-defined with the researcher). For each of the 10
experiments, identified a single representative starvation timepoint,
using the KG's `growth_phase` field where it distinguishes a starved
state, and the researcher's own knowledge of the source papers where it
doesn't (2 cases — see table below and Decisions). This required
**reopening step 3**: the background pool had been built pooling all
timepoints; it was rebuilt restricted to the same single-timepoint rule
(`3_analysis_framing/notebook.md` — Reopened section).

**Extraction** (`scripts/01_extract_target_de.py`): for all 61 named
target genes (all resolved genes from step 2, not just the 60 with
evidence — the 61st, unresolved-for-evidence gene still gets a row,
correctly showing "no data" everywhere), pulled DE data from each of the
10 experiments at its chosen timepoint, tagging each (gene, experiment)
as "matched" (gene's N/P annotation equals the experiment's nutrient) or
"cross" (it doesn't).

**Hypothesis computation** (`scripts/02_compute_hypotheses.py`): ran
`hit_rate()` on matched-N, matched-P, cross (both directions), and the 8
positive controls; ran `bootstrap_pvalue()` (10,000 iterations) and
Fisher's exact test on matched-N vs. the nitrogen background pool (H3,
nitrogen only per step 3).

**Main figure** (`scripts/03_heatmap_figure.py`): gene x experiment
heatmap, all 61 genes (rows) x all 10 experiments (columns, not
collapsed/averaged across studies — co-defined with the researcher after
checking concordance, see Surprises), colored by up / down / tested-not-
significant / no-data.

## Results

**Starvation timepoint per experiment** (locked with the researcher):

| # | Experiment | Chosen timepoint | Basis |
|---|---|---|---|
| 1 | Weissberg proteomics (MED4, N) | day 14 | `growth_phase=nutrient_limited`, before the later "death"-phase timepoints |
| 2 | Weissberg RNA-seq (MED4, N) | day 14 (only point) | — |
| 3 | Read RNA-seq (MED4, N) | 24h | only timepoint marked `nutrient_limited` (3h/12h are `acute_stress`) |
| 4 | Tolonen microarray (MED4, N) | 12h | `[user-provided]`: cultures started to decline at 12h per the source paper (KG marks all 6 timepoints `acute_stress` — no `nutrient_limited` signal exists in the KG for this rapid-onset design) |
| 5 | Tolonen microarray (MIT9313, N) | 12h | same as #4 |
| 6 | Lin RNA-seq (NATL2A, P) | **59h** | latest `nutrient_limited` timepoint in the arm that never had phosphate re-added; the `50h/59h (P added)` recovery timepoints are excluded. Corrected from 46h on 2026-09-08 — see Decisions ("Lin timepoint"). At 46h only 4 of Lin's 34 genes are significant; at 59h, 18 |
| 7 | Fuszard proteomics (MIT9312, P) | only point (stationary) | — |
| 8 | Fuszard proteomics (NATL2A, P) | only point (exponential) | flagged: not the same physiological state as #7 (culture collapsed before reaching stationary phase) |
| 9 | Martiny microarray (MED4, P) | 48h | `growth_phase=nutrient_limited`, matches the paper's own stated q<0.05-at-48h criterion |
| 10 | Martiny microarray (MIT9313, P) | 24h (latest available) | flagged: paper's stated "q<0.05 at t=48h" criterion has no matching 48h timepoint for this strain in the KG — included with a visible caveat (marked `*` in the figure) |

**Hit-rate results** (after the 2026-09-08 changes — step-2 gene-identity
corrections, Lin at 59h, Martiny table-absent reclassification; see
Decisions):

| Comparison | Distinct genes | Tests | Up | Down | Not sig. | % significant | % up | % down |
|---|---|---|---|---|---|---|---|---|
| H1: N-genes in N-experiments (matched) | 29 | 104 | 46 | 0 | 58 | 44.2% | 44.2% | 0.0% |
| H1: P-genes in P-experiments (matched) | 31 | 66 | 38 | 2 | 26 | 60.6% | 57.6% | 3.0% |
| H2: N-genes in P-experiments (cross) | 29 | 48 | 0 | 2 | 46 | 4.2% | 0.0% | 4.2% |
| H2: P-genes in N-experiments (cross) | 30 | 111 | 10 | 13 | 88 | 20.7% | 9.0% | 11.7% |
| Positive controls — N (ntcA, glnA, amtB/amt1, ureA) | 4 | 20 | 16 | 0 | 4 | 80.0% | 80.0% | 0.0% |
| Positive controls — P (pstS, phoA, phoB, phoR) | 4 | 14 | 12 | 1 | 1 | 92.9% | 85.7% | 7.1% |

Prior values (pre-2026-09-08): H1 N matched 41.3% (43/104); H1 P matched
75.7% (28/37); H2 N-in-P 100% (1/1); P controls 84.6% (11/13). The
H2 N-in-P jump from 1 test to 48 is the Martiny reclassification finally
making the phosphorus-experiment cross-tests visible; the H1 N rise
(41.3% → 44.2%) is the real `urtA` locus (step 2 had `urtE`'s).

**H3 — nitrogen noise/bootstrap check:** observed matched-N hit rate
**44.2%** vs. a background/null distribution with mean 18.8% (std 3.8%)
over 10,000 bootstrap iterations — **empirical p-value < 0.0001** (0 of
10,000 random size-matched draws reached it). Fisher's exact test on the
same comparison (46/104 significant in target genes vs. 1,337/7,047 in
the background pool): **odds ratio 3.39, p = 4.98e-09**. Both tests agree
— slightly firmer than the pre-2026-09-08 values (OR 3.01, p=2.15e-07).
The nitrogen background pool is unchanged by every 2026-09-08 fix (it is
built by gene *name* exclusion, and none of the phosphorus-side changes
touch the nitrogen `all_detected_genes` experiments).

**Full data:** `data/01_target_gene_experiment_matrix.csv` (610 rows, one
per gene x experiment), `data/02_hypothesis_results.csv`,
`data/02_bootstrap_null_distribution.csv`, `data/02_h3_summary.txt`.

**Main figure:** `figures/01_gene_experiment_heatmap.png` — revised after
researcher review: transposed (experiments as rows, genes as columns),
larger fonts, and restricted to the 40 of 61 genes that showed an
up/down response in at least 1 experiment (genes silent everywhere are
dropped from the figure but remain in the full data table). N and P
experiments in separate row blocks, N and P genes in separate column
blocks, the MIT9313-phosphorus significance-criterion caveat marked with
`*`. Revised again 2026-09-08 (see Decisions): (a) the single "no data"
category split into "no data (gene present, not in this experiment's
table)" and "gene absent from strain genome" (`no_locus_in_strain`), the
latter hatched; (b) the underlying matrix re-extracted after the step-2
gene-identity corrections, the Lin 59h timepoint, and the Martiny
table-absent reclassification. `scripts/07_heatmap_pptx.py` also emits
`figures/01_gene_experiment_heatmap.pptx` — the same figure with every
cell as an editable rectangle for hand-annotation in PowerPoint.

**Percentage-upregulated table and figure**
(`scripts/04_pct_upregulated.py`): for every gene, the fraction of its
**matched-nutrient** tested experiments (N genes in N experiments, P
genes in P experiments; cross-nutrient tests excluded — decided
2026-09-08, see Decisions) that came back significantly upregulated. Full
table: `data/04_pct_upregulated.csv`. Figure: `figures/02_pct_upregulated.png`.
7 nitrogen genes reach 100% under N starvation (`ntcA`, `glnA`, `cynA`,
`cynD`, `focA`, `nirA`, `nirX` — the last 3 tested in only 1 experiment
each). Several phosphorus genes reach 100% under P starvation, including
`phoA`, `phoB`, `phoE`, `pstA`, `pstC`, `gap3` (all 3/3) — `phoE` now
that its locus is corrected (step 2). Mean per-gene matched-nutrient
upregulation rate: 43.1% (N target genes), 50.4% (P target genes).

**Percentage strip plot, target vs. background** (`scripts/05_pct_scatter.py`,
`figures/03_pct_scatter.png`): the same per-gene matched-nutrient
percentages as one dot each, grouped as N target genes (mean **43.1%**),
N background pool (mean 6.1%, n=4,019 background genes by locus tag), and
P target genes (mean **50.4%**, no P background exists per step 3).
Restricting to matched-nutrient tests (2026-09-08) realigns the N
target-gene group with the matched-nitrogen H3 bootstrap population — the
43.1% per-gene mean now tracks the 44.2% pooled hit rate. N target-gene
dots sit far above the N background dots as a whole distribution, not
just in a summary statistic. Caveat: background genes are grouped by
locus tag, not a
cross-strain gene name like the target genes are (`04_pct_upregulated.csv`
pools e.g. `ntcA` across all 4 strains into one rate) — background genes
have no cross-strain identity to pool by, since each strain's background
sample is drawn independently. A MED4 background gene can be tested in up
to 3 of the 4 background experiments (proteomics/RNA-seq/microarray all
being MED4); a MIT9313 background gene only in 1 — this is why the
background dots visibly cluster at a few discrete values (0%, 33%, 50%,
67%, 100%) rather than spreading continuously like the target genes.

**Nutrient-collapsed heatmap** (`scripts/06_pct_heatmap_by_nutrient.py`,
`figures/04_pct_heatmap_by_nutrient.png`, added after researcher review of
the first two figures, then revised again 2026-08-12 after a second
researcher review — see Decisions): two stacked panels, one per gene
group — N-acquisition genes on top, P-acquisition genes below — each with
its own pair of rows ("N starvation experiments", "P starvation
experiments"). Cell = % of that gene's *tested* experiments (cells with no
data for that gene/nutrient combination are excluded from the denominator,
not counted as 0%) that came back significant in whichever direction
(up/down) had more hits; colored red (up) / blue (down) / gray (tested,
not significant) / hatched (no data), color intensity scaled continuously
by the percentage. The percentage is only printed as a number inside the
cell when it is >=55% — a readability threshold for white-on-color
legibility, not a data cutoff; cells below it are still colored by their
exact percentage. This makes the matched-vs-cross pattern (H1 vs. H2)
visible gene-by-gene in one image: N-acquisition genes are red/upregulated
in the top panel's "N starvation" row (matched) and mostly gray/hatched in
its "P starvation" row (cross), and the mirror pattern holds for
P-acquisition genes in the bottom panel — the same conclusion as the
H1/H2 hit-rate table, now visible per gene rather than only as an
aggregate rate. Full data: `data/06_pct_by_nutrient.csv`.

## Surprises

- **Duplicate-measurement concordance check** (co-defined with the
  researcher before deciding how to display repeated measurements):
  comparing the 9 MED4 nitrogen genes measured independently by both
  Weissberg (day14) and Read (24h), direction never disagreed between the
  two studies, and 6 of 9 genes agreed fully on both direction and
  significance. The 2 disagreements (`amt1`, `ureG`) were both borderline
  cases — same direction, broadly similar magnitude, but one study's
  significance call crossed the threshold and the other's didn't. This
  directly motivated **not** collapsing repeated measurements into an
  average (see Decisions) — an average would have hidden exactly this
  kind of informative disagreement.
- **`urtA`/`urtE` — RESOLVED 2026-09-08.** The original note here read:
  "`urtA`/`urtE` share a Cyanorak ID (`CK_00008074`), a source-list
  data-entry issue; both genes now appear with identical data, the honest
  representation given they are indistinguishable in this KG build." That
  was wrong — they are *not* indistinguishable. `CK_00008074` is `urtE`'s
  group; the spreadsheet gave `urtA` that ID by mistake. The real `urtA`
  group is `CK_00000076` (`PMM0970` in MED4, "substrate-binding
  component"), a distinct gene with distinct data. Caught by the step-2
  name-vs-Cyanorak cross-check (`2_kg_selection/scripts/06_name_vs_cyanorak_crosscheck.py`).
  After the fix, `urtA` (real) and `urtE` are separate columns with
  different values, and real `urtA` is more N-responsive — H1 nitrogen
  rate 41.3% → 44.2%. See `gaps_and_friction.md` (2026-09-08).
- **H2 cross-nutrient evidence for N-genes-in-P-experiments was nearly
  empty** (1 test) *before* the Martiny reclassification. As of 2026-09-08
  it is 48 tests (2 significant), because the Martiny microarray's
  genome-present-but-table-absent nitrogen genes are now counted as
  "tested, not significant" — so "N genes do not respond to P starvation"
  now rests on real evidence, not a data vacuum.

## Decisions

**2026-08-10 — One starvation timepoint per experiment, not pooled.**
Co-defined with the researcher (see Context); reopens step 3's background
pool. Full table above.

**2026-08-10 — Exclude phosphate-repletion ("recovery") timepoints
entirely**, not just from the figure. The Lin et al. NATL2A experiment's
2 `recovery` timepoints (P added back) measure reversal, not starvation —
counting them as starvation evidence would have diluted the real signal.
Confirmed by the researcher specifying timepoint 46h (last pre-repletion
point) as the chosen comparison.

**2026-08-10 — All 61 genes shown, not just significant ones.** A
"significant only" figure would hide real findings — particularly for H2,
where "this gene did not respond" is itself the answer.

**2026-08-10 — Repeated measurements (same gene, strain, platform,
different study) shown as separate columns, not averaged.** Motivated by
the concordance check above: averaging would hide genuine
study-to-study disagreement on borderline calls.

**2026-08-12 — Nutrient-collapsed heatmap redesigned into two stacked
panels, gene labels recolored black, "Gene" corner label added.**
Co-defined with the researcher during a walkthrough of the figure: the
original single-block layout (2 rows shared across N- and P-acquisition
gene columns side by side) was replaced with two independent panels — N
genes on top with their own 2-row block, P genes below with their own
2-row block — so each gene group's matched/cross rows sit together rather
than sharing rows with the other group. Gene-name column labels changed
from nutrient-colored (blue/red) to black for readability; the "N/P
acquisition genes" panel titles moved closer to their panels and changed
to black; a small "Gene" label was added above the row-label column in
each panel to identify what the columns are. No underlying data changed —
`data/06_pct_by_nutrient.csv` has the same computation, only the figure's
layout differs. This redoes part of step 5 (`scripts/06_pct_heatmap_by_nutrient.py`,
`figures/04_pct_heatmap_by_nutrient.png`); recorded here per the redo path
rather than as a silent edit.

**2026-09-08 — Figure 1: split the "no data" cell category into
"gene present but not in this experiment's table" and "gene absent from
the strain's genome", the latter hatched.** Researcher-directed, after
this analysis was closed. The extraction script already recorded these as
two distinct states (`no_data_at_timepoint` vs. `no_locus_in_strain` —
see `scripts/01_extract_target_de.py`); `scripts/03_heatmap_figure.py`
had been collapsing them to one color. The genome-absent cells the
researcher called out — `cynA/B/D/S`, `PMM707/719/721`, `psiP1`, `phoA`,
`ptrA`, `unkP5` absent from MIT9313; `nirA`/`nirX`, `focA` absent from
MED4 — were checked against `data/01_target_gene_experiment_matrix.csv`
and all were already `no_locus_in_strain`, consistent with the step-2
resolution. No data or hypothesis result changes: both non-tested states
were already excluded from `hit_rate()` (`scripts/02_compute_hypotheses.py`
counts only `significant_up`/`significant_down`/`not_significant`). Only
`scripts/03_heatmap_figure.py` and `figures/01_gene_experiment_heatmap.png`
change; recorded here per the redo path rather than as a silent edit.
(`scripts/07_heatmap_pptx.py` added the same day emits an editable
PowerPoint of the figure — `figures/01_gene_experiment_heatmap.pptx`,
every cell a named rectangle.)

**2026-09-08 — Re-extracted the whole matrix after the step-2 gene-identity
reopen.** Step 2 corrected 6 gene loci (4× `phoE`, `unkP2`, `urtA`; plus
`ptrA`-NATL2A). `scripts/01_extract_target_de.py` re-run picks up the new
resolution. Result changes: H1 N matched 41.3% → 44.2% (real `urtA`); H1 P
matched moved (see next two decisions for the other drivers); H3 p<0.0001
throughout (Fisher OR 3.01 → 3.39). Cascade approved by the researcher.

**2026-09-08 — Lin timepoint corrected 46h → 59h.** The 46h choice
("last `nutrient_limited` before phosphate re-added") misread the
experiment's two arms: `59h` is also `nutrient_limited` and is the latest
point in the arm that never got phosphate; `50h/59h (P added)` are the
recovery arm (still excluded). At 46h only 4 of Lin's 34 genes clear
significance; at 59h, 18 — the `pst` operon and `phoR` do not engage until
59h. KG `expression_status` verified as a genuine per-timepoint DESeq2
call. Researcher-directed. `gaps_and_friction.md` 2026-09-08.

**2026-09-08 — Martiny table-absent genes reclassified as "tested, not
significant" (Martiny only).** The 2 Martiny phosphorus microarrays have
`table_scope = filtered_subset` with a pure significance filter (q<0.05 at
48h) on a genome-wide microarray — so a genome-present target gene not in
the table was measured and did not pass. `RECLASSIFY_ABSENT_AS_NOT_SIG` in
`scripts/01_extract_target_de.py`. **Not** applied to Fuszard (iTRAQ — a
fixed detected-protein list; absence = not detected, verified against the
live KG) or Lin (curated operon-level table). Co-defined with the
researcher (chose "Martiny only" over "Martiny + Fuszard"). Effect: lowers
the phosphorus matched rate (removes inflated all-significant
denominators) and lifts H2 N-in-P from 1 test to 48.

**2026-09-08 — Figures 2 and 3 restricted to matched-nutrient tests.**
Both pool a gene's rate across every experiment it appears in; after the
Martiny reclassification, nitrogen genes picked up their Martiny
phosphorus cross-tests as "not significant," dragging the Figure 3 "N
target genes" group down to ~27% and out of step with the matched-nitrogen
H3 bootstrap it exists to illustrate. Restricting `scripts/04_pct_upregulated.py`
(and therefore `05_pct_scatter.py`) to `match_type == "matched"` puts it
back at 43.1%, tracking the 44.2% pooled hit rate. Cross-nutrient
behaviour stays visible in Figure 4, which splits by nutrient. Researcher
was offered pooled-vs-matched and did not object to the recommended
matched-only.

## Decide-gate checklist

- **Outputs produced:** `scripts/01_extract_target_de.py` (Lin 59h +
  Martiny reclassification), `scripts/02_compute_hypotheses.py`,
  `scripts/03_heatmap_figure.py` (5-category legend + hatching),
  `scripts/04_pct_upregulated.py` (matched-nutrient only),
  `scripts/05_pct_scatter.py` (matched-nutrient only),
  `scripts/06_pct_heatmap_by_nutrient.py`, `scripts/07_heatmap_pptx.py`
  (new — editable PowerPoint); `data/01_target_gene_experiment_matrix.csv`,
  `data/02_hypothesis_results.csv`,
  `data/02_bootstrap_null_distribution.csv`, `data/02_h3_summary.txt`,
  `data/04_pct_upregulated.csv`, `data/05_pct_scatter_background.csv`,
  `data/06_pct_by_nutrient.csv`;
  `figures/01_gene_experiment_heatmap.png` + `.pptx`,
  `figures/02_pct_upregulated.png`, `figures/03_pct_scatter.png`,
  `figures/04_pct_heatmap_by_nutrient.png`. Step 3 reopened (2026-08-10):
  `3_analysis_framing/scripts/02_build_controls.py` updated,
  `3_analysis_framing/data/02_negative_background_pool.csv` regenerated
  (unaffected by the 2026-09-08 changes — nitrogen `all_detected_genes`
  only, excluded by gene name).
- **Results presented:** timepoint table, hit-rate table, H3
  bootstrap/Fisher results, all four figures (above).
- **QC gate (original):** gene-count cross-check surfaced the
  `urtA`/`urtE` locus collision → fixed. **QC gate (2026-09-08 reopen):**
  step-2 paralog audit + name cross-check (see `2_kg_selection/notebook.md`);
  H3 nitrogen background confirmed unchanged by all phosphorus-side fixes;
  matrix diff after the Martiny reclassification confirmed to touch only
  the 2 Martiny experiments' rows (line-ending artifact ruled out).
- **Advance rationale:** all three hypotheses recomputed on the
  identity-corrected matrix with the corrected Lin timepoint; every figure
  regenerated; the nitrogen headline result is firmer, not weaker; caveats
  carried to `gaps_and_friction.md`. Ready for step 6 re-evaluation.
