# Gaps and friction log

Append-only. Each entry: date, short name, what happened, workaround/impact if any.

**2026-08-10 — `growth_phase` doesn't always identify the "starvation" timepoint, and one experiment's stated significance criterion doesn't match its available data.** While selecting a single representative starvation timepoint per experiment (step 5), two gaps surfaced: (1) Tolonen et al. 2006's rapid-onset nitrogen-deprivation time courses (MED4 and MIT9313 microarray) never reach `growth_phase = nutrient_limited` at any of their 6 timepoints (0-48h) — every timepoint is labeled `acute_stress`, so the KG alone gives no signal for "when did starvation take hold." Resolved with researcher-supplied knowledge from the source paper (cultures began declining at 12h), not a KG field — tagged `[user-provided]` in `5_analyze/notebook.md`. (2) Martiny et al. 2006's MIT9313 phosphorus microarray experiment carries the description "genes with q<0.05 at t=48h" (copied from the same paper's MED4 entry), but the KG's actual timepoint data for MIT9313 stops at 24h — no 48h row exists. Cannot determine from the KG whether this is a missing-data gap or whether the description was never updated for MIT9313. Included in the analysis with a visible caveat rather than resolved. Impact: any future work with `growth_phase` should not assume it marks "starvation achieved" without checking whether the experiment is a rapid-onset design; any use of a `table_scope_detail` string should be verified against the experiment's actual timepoint list before trusting it strain-by-strain within a multi-strain publication.

**2026-08-12 — The nitrogen "background/noise rate" is not one number; it
ranges from 3.0% to 47.2% across the 4 background experiments.** While
evaluating H3 in step 6, breaking the pooled 18.8% background mean out
per experiment showed a ~16-fold spread (Tolonen MIT9313 microarray 3.0%,
Weissberg proteomics 9.1%, Tolonen MED4 microarray 17.6%, Weissberg
RNA-seq 47.2% — this last one higher than the observed matched-nitrogen
target rate itself, 41.3%). The bootstrap in `4_methods/np_response.py`
already handles this correctly (it draws each iteration's random sample
from the same per-experiment pools the real target genes were tested in,
not from one pooled rate), so the H3 result is not compromised — but the
single pooled mean reported in step 5's notebook understates how much the
"noise" floor varies by platform/experiment. Impact: future analyses
using this bootstrap pattern should report (or at least check) the
per-experiment background rate, not just the pooled mean/std, since a
single "background rate" can be misleading on its own even when the
underlying resampling is done correctly.

**2026-08-10 — MIT9312 phosphorus evidence collapses under `significant_only` table_scope.** Step 2 resolved 44 of the researcher's 92 genes to a locus tag in MIT9312, but only 2 have actual DE evidence. MIT9312's sole in-scope experiment (Fuszard et al. 2012 proteomics, `10.1186/2046-9063-8-7`) has `table_scope = significant_only` — the source publication's supplementary table only lists genes that passed its own fold-change cutoff (>1.6 or <0.6), so every gene not already significant in that paper is invisible to this KG build, not "tested and not significant." Impact: MIT9312 will contribute almost no data to steps 4-6 regardless of framing choices, and any background/noise gene set for MIT9312 would be a biased (pre-filtered) sample, not a random draw from the tested genome — same caveat the KG's `table_scope` field is designed to surface. Carry this into step 3 framing when deciding whether MIT9312 supports a background/noise comparison.
