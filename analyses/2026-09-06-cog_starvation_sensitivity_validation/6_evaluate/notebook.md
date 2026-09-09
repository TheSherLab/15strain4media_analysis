# Step 6 — Evaluate

## Context

Step 5 (redo 2026-09-09) computed the hit rates on the full 41-COG gene
set across all 10 experiments. This step assesses the result against the
step-3 framing, runs the robustness checks the "use every gene per COG"
decision raised, and harvests caveats.

## What I did

**Sensitivity checks** (`scripts/01_sensitivity_checks.py`): the
nitrogen headline (15.5%, below background) rerun three ways —
(a) excluding the 3 COGs with >8 genes in a strain (COG0477, COG0697,
COG0845); (b) with each COG collapsed to one value per experiment
(a COG "responds" if ≥50% of its genes are significant), so a 30-gene COG
counts once like a 1-gene COG; and (c) the per-COG nitrogen
significant-gene fractions for the caveats.

## Results

**Nitrogen robustness** (`data/01_sensitivity_summary.txt`; background
overall rate **19.1%**):

| Cut | Significant / tests | Rate | Bootstrap p | Fisher OR (p) |
|---|---|---|---|---|
| all genes, all 41 COGs (step 5 headline) | 68 / 438 | 15.5% | 0.930 | 0.78 (0.067) |
| excluding the 3 big COGs (COG0477/0697/0845) | 37 / 267 | 13.9% | 0.965 | 0.68 (0.031) |
| per-COG unit (≥50%-of-genes threshold) | 18 / 121 | 14.9% | — | — |

Every cut lands **at or below** the 19.1% background. The result does not
depend on the big COGs and does not depend on counting genes vs. counting
COGs.

**Top nitrogen "responsive" COGs** (`data/02_per_cog_nitrogen_fractions.csv`):
COG1596 (Wza, 40%, 2/5), COG1787 (Mrr/RM, 33%, 1/3), COG2192
(carbamoyltransferase, 29%, 2/7), COG0188 (gyrase/topo, 27%, 7/26) — all
`N`-tagged, all small, all with roughly balanced up/down.

## Evaluation against the step-3 framing

**Nitrogen response elevated above background — NOT supported.** The
pooled 41-COG rate (15.5%) is below the random background (19.1%);
bootstrap p = 0.93, Fisher's exact p = 0.067 (trend below, not
significant). The `N`-tagged subset (23.2%) trends above background but is
squarely inside the bootstrap null (p = 0.44, Fisher OR 1.28, p = 0.30).
The `mixed` subset (13.3%) is significantly *below* background (Fisher
p = 0.0067). `[interpretation]` These are general cellular-maintenance
genes — chaperones, DNA repair, transporters, LPS. As a group they respond
to nitrogen starvation *no more than*, and the `mixed` half *less than*, a
random gene. This is the same conclusion as the 2026-09-06 run (16.2% vs
19.0%), now on the full COG gene set and robust to every sensitivity cut.

**Phosphorus response — flat.** 5.3% significant (13/246), 8 up / 5 down,
descriptive only. `[interpretation]` The candidate list shows no
phosphorus-starvation response in the one genome-wide phosphorus dataset
(Martiny); Lin and Fuszard add almost nothing. Consistent with — and
weaker than — the already-weak nitrogen signal.

**What this means for the comparative-genomics finding.**
`[interpretation]` The 41 COGs were selected because their cross-strain
*copy number* tracks the starvation-sensitivity phenotype. This analysis
asked whether that genomic signal has in-vivo expression support: does
starvation actually change how much these genes are transcribed/translated?
The answer is **no** — not above what any random gene shows. This does not
refute the copy-number association (gene *dosage* differences between
strains can matter without the genes being starvation-*regulated*), but it
means the expression data give the copy-number finding no independent
support, and the two lines of evidence should be reported as separate,
not mutually reinforcing.

No preregistration was written (step 3 stated the hypothesis in prose).
This evaluation compares the pooled and subset hit rates, and their
background tests, to what step 3 said would count as support.

## Surprises

- **The result got slightly *more* negative on the full gene set**, not
  less. Adding dnaK1/2/3 and the other real family members (the point of
  the step-2 reopen) did not lift the hit rate — those chaperones are not
  strongly nitrogen-starvation-regulated in these experiments either.
- **COG0188 (gyrase/topoisomerase) at 27% (7/26)** is the only big-ish COG
  above background, and it is `N`-tagged. DNA topology management under
  starvation is plausible `[interpretation]`, but 7/26 with a balanced
  3 up / 4 down split is not a clean induction signal.

## Caveats

- **The pooled rate mixes very unequal COGs.** COG0477 alone contributes
  ~112 of the 438 nitrogen tests (25–34 MFS permeases per strain). The
  "use every gene" choice matches the copy-number test's gene set but
  means a few broad families weight the pooled number. The sensitivity
  checks (drop-big-COGs, per-COG unit) are why the headline is trustworthy
  despite this.
- **Nitrogen strain coverage is 2 strains** (MED4, MIT9313); phosphorus
  effectively 2 (MED4, MIT9313 via Martiny), with Lin/Fuszard contributing
  a handful of table genes. Narrow relative to the 13–15 strains the
  copy-number test spanned.
- **A few genes carry 2–3 COG numbers** (multi-domain proteins like zipN);
  their expression is counted once per COG. Small non-independence, carried
  from how the copy-number test counted them.
- **Individual (gene × experiment) tests are not independent** — the same
  gene is tested across experiments/strains. The bootstrap and Fisher's
  test address this at the aggregate for the nitrogen comparison; the
  per-COG fractions are descriptive.
- **No phosphorus background test is possible** — the phosphorus 5.3% is
  descriptive, not a confirmed "below background".
- **Significance criteria differ by platform** (padj+fold-change vs
  padj-only vs fold-change-only), carried unchanged from the walkthrough;
  the bootstrap/Fisher correct for this for nitrogen only.

## Decisions

**2026-09-09 — Report the copy-number finding and the expression result
as separate, non-reinforcing lines of evidence.** The expression data do
not support the comparative-genomics selection; they also do not refute it
(dosage ≠ regulation). Co-defined framing for the paper's Discussion.

## Decide-gate checklist

- **Outputs produced:** `scripts/01_sensitivity_checks.py`;
  `data/01_sensitivity_summary.txt`, `data/02_per_cog_nitrogen_fractions.csv`.
- **Results presented:** robustness table, top-responsive-COG list,
  evaluation against framing, caveats (all above).
- **QC gate:** all three nitrogen cuts land at/below background —
  conclusion is not an artifact of the big COGs or the per-gene unit.
  Bootstrap reused seed 42 / 10,000 iterations.
- **Decisions made this step:** report the two evidence lines separately
  (2026-09-09).
- **Advance rationale:** all hypotheses evaluated, the negative result is
  robust to the sensitivity checks, caveats harvested, `paper.md`
  Discussion written. Analysis complete.
