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
| 6 | Lin RNA-seq (NATL2A, P) | 46h | last `nutrient_limited` timepoint before phosphate was added back (`recovery` timepoints excluded entirely — see Decisions) |
| 7 | Fuszard proteomics (MIT9312, P) | only point (stationary) | — |
| 8 | Fuszard proteomics (NATL2A, P) | only point (exponential) | flagged: not the same physiological state as #7 (culture collapsed before reaching stationary phase) |
| 9 | Martiny microarray (MED4, P) | 48h | `growth_phase=nutrient_limited`, matches the paper's own stated q<0.05-at-48h criterion |
| 10 | Martiny microarray (MIT9313, P) | 24h (latest available) | flagged: paper's stated "q<0.05 at t=48h" criterion has no matching 48h timepoint for this strain in the KG — included with a visible caveat (marked `*` in the figure) |

**Hit-rate results:**

| Comparison | Distinct genes | Tests | Up | Down | Not sig. | % significant | % up | % down |
|---|---|---|---|---|---|---|---|---|
| H1: N-genes in N-experiments (matched) | 29 | 104 | 43 | 0 | 61 | 41.3% | 41.3% | 0.0% |
| H1: P-genes in P-experiments (matched) | 18 | 37 | 26 | 2 | 9 | 75.7% | 70.3% | 5.4% |
| H2: N-genes in P-experiments (cross) | 1 | 1 | 0 | 1 | 0 | 100% | 0% | 100% |
| H2: P-genes in N-experiments (cross) | 30 | 109 | 9 | 13 | 87 | 20.2% | 8.3% | 11.9% |
| Positive controls — N (ntcA, glnA, amtB/amt1, ureA) | 4 | 20 | 16 | 0 | 4 | 80.0% | 80.0% | 0.0% |
| Positive controls — P (pstS, phoA, phoB, phoR) | 4 | 13 | 10 | 1 | 2 | 84.6% | 76.9% | 7.7% |

**H3 — nitrogen noise/bootstrap check:** observed matched-N hit rate
41.3% vs. a background/null distribution with mean 18.8% (std 3.8%) over
10,000 bootstrap iterations — **empirical p-value < 0.0001** (0 of 10,000
random size-matched draws reached 41.3%). Bootstrap wall-clock time:
26.7 seconds (as expected — no KG calls in the loop). Fisher's exact test
on the same comparison (43/104 significant in target genes vs. 1,337/7,047
in the background pool): odds ratio 3.01, p = 2.15e-07. Both tests agree:
the nitrogen matched-response rate is not explainable by chance.

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
`*`.

**Percentage-upregulated table and figure**
(`scripts/04_pct_upregulated.py`): for every gene, the fraction of its
*tested* experiments (excluding no-data cells) that came back
significantly upregulated. Full table: `data/04_pct_upregulated.csv`.
Figure: `figures/02_pct_upregulated.png`. 6 nitrogen genes reach 100%
(upregulated in every experiment that tested them): `ntcA`, `cynA`,
`cynD`, `focA`, `nirA`, `nirX` — though `focA`, `nirA`, `nirX` were each
only tested in 1 experiment, so "100%" there means "up the one time it
was checked," not a robust rate. No phosphorus gene reaches 100%; the
highest is `PMM719` at 80% (4 of 5 tests).

**Percentage strip plot, target vs. background** (`scripts/05_pct_scatter.py`,
`figures/03_pct_scatter.png`): the same per-gene percentages as one dot
each, grouped as N target genes (mean 40.5%), N background pool (mean
6.1%, n=4,019 background genes by locus tag), and P target genes (mean
19.4%, no P background exists per step 3). Gives a visual complement to
the H3 bootstrap p-value: N target-gene dots sit visibly higher than the
N background dots as a whole distribution, not just in a single summary
statistic. Caveat: background genes are grouped by locus tag, not a
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
the first two figures): 2 rows ("N starvation experiments", "P starvation
experiments") x all 61 genes (columns, grouped N-acquisition then
P-acquisition), cell = that gene's % of tests significant under that
row's nutrient specifically, colored red (up) / blue (down) / gray
(tested, not significant) / hatched (no data), color intensity scaled by
the percentage. This makes the matched-vs-cross pattern (H1 vs. H2)
visible gene-by-gene in one image: the N-acquisition columns are
red/upregulated in the "N starvation" row and mostly gray/hatched in the
"P starvation" row, and the mirror pattern holds for P-acquisition
columns — the same conclusion as the H1/H2 hit-rate table, now visible
per gene rather than only as an aggregate rate. Full data:
`data/06_pct_by_nutrient.csv`.

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
- **`urtA`/`urtE` share a Cyanorak ortholog-group ID in the researcher's
  source spreadsheet** (`CK_00008074` for both), so they resolve to the
  identical locus tag in every strain — a data-entry issue in the source
  list, not a KG problem (urtA and urtE are different proteins
  biologically: the urea ABC transporter's substrate-binding component
  vs. its ATPase component). First draft of the extraction script
  silently dropped one of the two genes because of this collision (a dict
  keyed by locus tag, overwritten when two gene names mapped to the same
  locus); caught by a gene-count mismatch (60 vs. the expected 61) and
  fixed by keying on (gene, locus) pairs instead of locus alone. Both
  genes now appear in the figure with identical data, which is the
  honest representation given they are indistinguishable in this KG
  build.
- **H2 cross-nutrient evidence for N-genes-in-P-experiments is nearly
  empty** (1 test total) — consistent with the P experiments' narrow,
  pre-filtered tables (step 3) leaving almost no room for genes outside
  their own significant/filtered gene lists, let alone off-topic
  (nitrogen) genes.

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

## Decide-gate checklist

- **Outputs produced:** `scripts/01_extract_target_de.py`,
  `scripts/02_compute_hypotheses.py`, `scripts/03_heatmap_figure.py`,
  `scripts/04_pct_upregulated.py`; `data/01_target_gene_experiment_matrix.csv`,
  `data/02_hypothesis_results.csv`,
  `data/02_bootstrap_null_distribution.csv`, `data/02_h3_summary.txt`,
  `data/04_pct_upregulated.csv`; `scripts/05_pct_scatter.py`,
  `data/05_pct_scatter_background.csv`; `scripts/06_pct_heatmap_by_nutrient.py`,
  `data/06_pct_by_nutrient.csv`; `figures/01_gene_experiment_heatmap.png`,
  `figures/02_pct_upregulated.png`, `figures/03_pct_scatter.png`,
  `figures/04_pct_heatmap_by_nutrient.png`. Step 3 reopened:
  `3_analysis_framing/scripts/02_build_controls.py` updated,
  `3_analysis_framing/data/02_negative_background_pool.csv` regenerated.
- **Results presented:** timepoint table, hit-rate table, H3
  bootstrap/Fisher results, heatmap (all above).
- **QC gate:** gene count cross-checked (61 expected vs. 60 initially
  produced) → traced to the urtA/urtE Cyanorak-ID collision → fixed and
  re-verified (61 genes in final output). Bootstrap runtime measured and
  reported as promised to the researcher (26.7s, matches the
  no-KG-calls-in-the-loop expectation from step 4).
- **Advance rationale:** all three hypotheses are computed with concrete
  numbers, the main figure is built and reviewed, and every visualization
  judgment call the researcher raised (gene inclusion, timepoints,
  duplicates) has been resolved with evidence, not assumption. Ready for
  step 6 to evaluate against the step-1 framing and harvest caveats.
