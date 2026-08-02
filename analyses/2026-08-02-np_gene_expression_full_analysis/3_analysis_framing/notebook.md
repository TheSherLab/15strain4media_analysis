# Step 3 — Analysis framing

## Context

Step 2 enumerated the KG entries this analysis runs against. This step
turns the locked question into three testable hypotheses, selects the
positive and negative control gene sets, and states what result would
count as support or non-support for each hypothesis, in KG-operational
terms.

## Co-define

Researcher proposed three hypotheses directly (rather than Claude
proposing framing options): a matched-direction hypothesis, a
cross-nutrient hypothesis, and a noise/method-limitation check using a
background gene set. Claude proposed operationalizing the noise check
with a random background sample and confirmed two judgment calls with
the researcher before building it:
- Background set size and "not starvation-related" filter: **100 genes
  per strain**, sampled from genes actually tested in this analysis's
  N/P experiments, excluding the 54 target genes and any gene whose
  product/name text matches an N/P-metabolism keyword (see What I did).
- Background genes need not share identity across strains (independent
  random draw per strain), since the comparison is "does this strain's
  background rate look like this strain's target-gene rate," not a
  cross-strain matched design.

A second data-availability gap surfaced while building the background
set (see Surprises) and reopened step 2 a second time.

## What I did

- `01_select_background_genes.py` — for each strain, pulled every tested
  gene (`differential_expression_by_gene(significant_only=False,
  verbose=True)`) across that strain's N/P experiments, excluded the 54
  target loci and any N/P-metabolism-keyword-matching gene, and randomly
  sampled up to 100 (seed=42).
- Selected **positive controls** by hand from the 54-gene list:
  canonical, well-established nutrient-starvation marker genes present
  in this analysis's own target set — `[interpretation]`, domain
  knowledge, not derived from the KG, but every gene named is a KG
  entity already in the target set with resolved loci (verified in step
  2).
  - Phosphorus: `pstS`, `phoA`, `phoB`, `phoR` (core Pho-regulon genes —
    phosphate-binding transporter subunit, alkaline phosphatase, and
    the PhoB/PhoR two-component regulatory pair).
  - Nitrogen: `ntcA` (master nitrogen-regulatory transcription factor),
    `glnA` (glutamine synthetase, the canonical N-assimilation enzyme),
    `amtB/amt1` (ammonium transporter), `ureA` (urea-uptake operon).
  - All 8 resolve to a locus in all 4 relevant strains except `phoA`,
    which resolves in MED4, MIT9312, NATL2A (not MIT9313 — consistent
    with step 2's resolution).

## Hypotheses

**H1 — matched, directional.** Nitrogen-annotated genes are
upregulated (higher expression under starvation than under replete
control) under nitrogen starvation; phosphorus-annotated genes are
upregulated under phosphorus starvation. Operationalized: among
(gene x strain x matched-nutrient experiment x timepoint) tests, what
fraction return `expression_status == significant_up`?

**H2 — cross-nutrient.** Under nitrogen starvation, phosphorus-annotated
genes are also differentially expressed (either direction); under
phosphorus starvation, nitrogen-annotated genes are also differentially
expressed. Operationalized: among (gene x strain x cross-nutrient
experiment x timepoint) tests, what fraction return
`expression_status` in `{significant_up, significant_down}`? No
directional prediction — this is deliberately open, since a shared
general-stress response and a nutrient-specific response would both
produce "differentially expressed" but for different reasons.

**H3 — noise / method-limitation check.** The apparent signal in H1/H2
could be an artifact of the significance-calling method rather than
biology specific to N/P-acquisition genes, if a similarly-sized random
set of unrelated genes shows a comparably high rate of "significant"
hits in the same experiments. Operationalized: compare the target
genes' significant-hit rate (H1/H2) against the same rate computed for
the 303-gene background set (up to 100 per strain, `data/01_background_genes.csv`)
in the same experiments.

## Target

- 53/54 genes with a resolved locus in >=1 of the 4 relevant strains
  (`ptxB/phnD2` has none — see Surprises); 29 nitrogen-annotated, 24
  phosphorus-annotated with usable loci.
- 4 strains: MED4, MIT9312, MIT9313, NATL2A.
- 10 experiments (5 nitrogen, 5 phosphorus — narrower than step 2's
  first pass; see Surprises) from `2_kg_selection/data/02_np_experiments.csv`.
- Gene→locus mapping from `2_kg_selection/data/03_gene_loci.csv`.

## Positive control

8 genes (4 N, 4 P) from the target set itself, selected as canonical
literature markers (see What I did). Expectation: these should show
`significant_up` under their matched nutrient starvation at a high rate
— if they don't, the DE-detection method or experiment selection has a
problem, independent of whether the broader 54-gene hypothesis holds.

## Negative control

303 genes (100 MED4, 32 MIT9312, 100 MIT9313, 71 NATL2A — MIT9312 and
NATL2A capped below 100 by a small tested-gene pool; see Surprises),
not annotated as N/P-related, drawn from genes actually tested in the
same 10 experiments. Expectation: these should show a **low**
significant-hit rate relative to the target genes' H1 rate — if the
background rate is comparably high, H1/H2 findings should be read as
weak evidence of nutrient-specific regulation, not strong evidence.

## Expected outcome (operational)

Computed in step 5 via a shared method module (step 4):
1. Matched-nutrient significant-hit rate for target genes (H1), split
   up/down.
2. Cross-nutrient significant-hit rate for target genes (H2), split
   up/down.
3. Same two rates for the background set (H3 comparison).

No numeric threshold is preregistered for "supported" vs "not
supported" — step 6 will describe the actual rates and directions
first, then interpret against these three hypotheses in prose.

## Surprises

- **MIT9301 dropped entirely.** While building the background gene
  pool, `differential_expression_by_gene` returned zero rows for
  MIT9301's two phosphorus experiments — both are `METABOLOMICS`
  (metabolite abundance), which carries no gene-level DE edges in the
  KG. Since MIT9301 had no nitrogen experiments either, this drops it
  from the analysis entirely. Raised with the researcher; decided to
  restrict all experiments (not just MIT9301's) to gene-expression/
  proteomics omics types (RNASEQ, PROTEOMICS, MICROARRAY), which is a
  second redo of step 2 (see that step's notebook). Strain count is now
  4 (MED4, MIT9312, MIT9313, NATL2A), and experiment count is 10 (5 N /
  5 P — one fewer P experiment than the prior 7, since MIT9301's 2
  metabolomics rows are gone alongside a phosphorus proteomics
  experiment that had also counted MIT9301... actually verify against
  `02_np_experiments.csv`: the 2 dropped P rows are exactly MIT9301's
  2 metabolomics experiments).
- **`ptxB/phnD2` has no locus in any of the 4 remaining strains.** Its
  only prior evidence (from the original 91-gene triage) was in
  MIT9301, which is now out of scope. It is excluded from the target
  set for this analysis, not "unresolved" — a scope consequence, not a
  resolution failure.
- **MIT9312 and NATL2A have small tested-gene pools** (77 and 127
  genes respectively, vs. ~1,900–2,200 for MED4/MIT9313) because their
  only phosphorus experiment is a targeted iTRAQ proteomics panel,
  which quantifies far fewer proteins than RNA-seq/microarray. This
  caps their background-set size below 100 and means their H1/H2 rates
  for these two strains are estimated from a much smaller universe of
  measurable genes — a caveat to carry into step 6, not a defect in the
  sampling.

## Decisions

- **2026-08-02** — Analysis restricted to gene-expression/proteomics
  omics types only (RNASEQ, PROTEOMICS, MICROARRAY); metabolomics
  experiments excluded because they carry no gene-level DE data. This
  is a second redo of step 2's experiment list and gene-locus
  resolution (first redo dropped SS120; this one drops MIT9301 and its
  2 metabolomics experiments).
- **2026-08-02** — Positive/negative control gene sets fixed as
  described above; no numeric significance threshold preregistered for
  H1/H2/H3 — step 6 describes rates first, interprets second.

## Decide-gate checklist

- **Outputs produced** — `scripts/01_select_background_genes.py` →
  `data/01_background_genes.csv` (303 rows). Positive-control gene list
  recorded in this notebook (not a separate data file — it's 8 rows
  drawn directly from step 2's `03_gene_loci.csv`).
- **Results presented** — hypotheses, target/control definitions, and
  background-pool-size table shown inline above, matching what was
  shown to the researcher in chat.
- **QC gate** — confirmed all 8 positive-control genes resolve to a
  locus in their relevant strains (7/8 in all 4, `phoA` in 3/4,
  consistent with step 2); confirmed background genes are drawn only
  from each strain's actually-tested gene pool (not the full genome) so
  the H1-vs-H3 comparison is apples-to-apples; confirmed 0 target loci
  leaked into the background sample (explicit exclusion check in the
  script).
- **Decisions made this step** — see Decisions above.
- **Advance rationale** — three hypotheses are stated in operational,
  KG-measurable terms; positive and negative control gene sets are
  selected and validated; ready for step 4 (methods: the shared
  significant-hit-rate computation used across target, positive-control,
  and background sets).
