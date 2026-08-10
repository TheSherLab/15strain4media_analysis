# Step 1 — Research question

## Context

Researcher supplied a 92-gene list (nitrogen/phosphorus acquisition genes,
compiled from prior literature review; `Dataset 2.xlsx`, sheet "N and P
acquisition genes", the "15 strains 4 media" paper revision) and asked
whether these genes' expression has been tested in vivo under N or P
starvation vs. replete media, across differential expression, proteomics,
and microarray data. A prior analysis
(`analyses/2026-08-02-np_gene_expression_full_analysis/`) already answered
a closely related question; the researcher chose to redo the process from
scratch as a hands-on walkthrough (not because the prior analysis's
decisions were considered wrong), so this notebook makes its own
independent calls rather than importing the prior ones.

## KG context

Grounding queries run before locking scope:

- `kg_release_info` — KG version `0.1.0-alpha.6` (production build,
  2026-06-16), verdict `ok`. 124,751 genes, 197 experiments, 43
  publications, 47 organisms in this build.
- `list_experiments(summary=true, organism="Prochlorococcus")` — 91
  Prochlorococcus experiments total; by treatment type, 12 nitrogen and 9
  phosphorus (next largest: 18 carbon, 10 light).
- `list_experiments(organism="Prochlorococcus", treatment_type=["nitrogen"], verbose=true)` — 12 experiments, 4 publications (MED4, MIT9313, SS120).
- `list_experiments(organism="Prochlorococcus", treatment_type=["phosphorus"], verbose=true)` — 9 experiments, 4 publications (MED4, MIT9312, MIT9313, NATL2A, SS120, MIT9301).

Structural surprise: not every "nitrogen" or "phosphorus" treatment-type
experiment is a literal starvation-vs-replete contrast. The nitrogen set
includes 3 experiments comparing growth on an alternate N source (cyanate,
urea, or nitrite as sole N source) vs. replete medium — a substrate switch,
not deprivation — and 1 using a chemical inhibitor (azaserine) as an
N-limitation proxy rather than actual media manipulation. The phosphorus
set includes 2 metabolomics-only experiments with no gene-level
differential expression, and one RNA-seq pair confounded with phage
infection status. Background factors also vary: several experiments run in
coculture with a bacterial partner (Alteromonas) rather than axenic
monoculture.

## What I did

Reviewed the researcher's shared gene-reference table (92 genes, N/P
acquisition, with prior KG-evidence resolution status) and the existing
`2026-08-02` analysis's `paper.md` to understand what had already been
attempted. Surfaced this to the researcher before doing any new work, since
redoing it would duplicate existing effort; researcher chose to redo as an
independent walkthrough.

Ran the KG grounding queries above, then worked through the following
clarifying questions with the researcher:

1. **Hypotheses to test** — researcher selected: (a) matched response
   (N-annotated genes respond to N starvation, P-annotated genes to P
   starvation), (b) cross-nutrient response (do genes respond to the
   non-matching nutrient's starvation too — is the response
   nutrient-specific or general stress?), and (c) a noise/robustness check:
   is any observed pattern distinguishable from measurement noise or
   background rate, evaluated with significance testing and bootstrap
   resampling (not just point-estimate rates).
2. **Confounded backgrounds** (coculture, phage infection) — researcher
   chose to restrict the main comparison to axenic, uninfected experiments
   only, for the cleanest causal read (any expression change attributable
   to N/P status alone, not a partner organism or infection).
3. **Alternate-N-source experiments** (cyanate/urea/nitrite vs. replete) —
   researcher chose to exclude these from "starvation vs. replete": they
   test substrate preference, not nutrient deprivation.

## Decisions

**2026-08-10 — Organism scope restricted to the researcher's 15 study
strains.** The researcher's source spreadsheet covers a "15 strains 4
media" study; the researcher asked to restrict this analysis to those same
15 Prochlorococcus strains: MED4, MIT9312, MIT9313, MIT1327, MIT0604,
NATL2A, MIT9515, MIT9215, AS9601, PAC1, MIT9202, SB, MIT9301, MIT1314,
NATL1A. Checked all 15 against the KG (`list_organisms`,
`organism_names=[...]`): all 15 exist as `genome_strain` entries, but only
7 have any experiments at all in this KG build (the other 8 — MIT1327,
MIT0604, MIT9202, MIT9215, MIT9515, PAC1, MIT1314, SB — are genome-only,
zero experiments, zero publications). Of the 7 with experiments, treatment
coverage is uneven: MED4 and MIT9313 have both nitrogen- and
phosphorus-treatment experiments; MIT9312, MIT9301, and NATL2A have
phosphorus-treatment experiments only (no nitrogen); AS9601 and NATL1A
have experiments but only for salt stress, not N or P. Net effect: the
nitrogen side of this analysis can only ever draw on MED4 and MIT9313; the
phosphorus side can draw on MED4, MIT9312, MIT9313, MIT9301, and NATL2A.
This also resolves the SS120 azaserine experiment question from the KG
grounding above — SS120 is not one of the researcher's 15 strains, so it is
excluded by this scope regardless of the proxy-vs-literal-starvation
question.

**2026-08-10 — Redo as independent walkthrough, not reuse.** The
researcher wants to experience each step's judgment calls directly (for
understanding and for the supplementary figure's methods narrative) rather
than inherit the prior analysis's already-made decisions. This analysis
does not read from or depend on `2026-08-02-np_gene_expression_full_analysis/`;
any similarity in conclusions is a cross-check, not a dependency.

**2026-08-10 — Scope restricted to axenic, uninfected, literal
starvation/limitation contrasts.** Coculture and phage-infected background
experiments, and alternate-N-source substrate-switch experiments, are
excluded from the main comparison. This is a deliberate precision-over-recall
choice: it shrinks the usable experiment pool but keeps every included
contrast interpretable as "N or P withheld/limited vs. replete" with no
second variable riding along. (Step 3 will confirm how many experiments and
genes survive this filter — if the phosphorus side becomes too thin to
support the noise/bootstrap check, that will be revisited.)

## Locked research question

For the 92 nitrogen- and phosphorus-acquisition genes in the researcher's
reference list, do genes annotated for a given nutrient show a
differential-expression response when axenic, uninfected Prochlorococcus —
restricted to the researcher's 15 study strains (MED4, MIT9312, MIT9313,
MIT1327, MIT0604, NATL2A, MIT9515, MIT9215, AS9601, PAC1, MIT9202, SB,
MIT9301, MIT1314, NATL1A) — is grown under literal starvation/limitation of
that nutrient (N or P withheld or limited in the medium) versus
nutrient-replete medium — and is that response nutrient-specific (matched:
N-genes respond to N starvation, P-genes to P starvation) or general
(cross: genes respond regardless of which nutrient is withheld)? Evidence
sources: RNA-seq, proteomics, and microarray differential expression in the
KG. A third question asks whether any observed matched/cross pattern is
statistically distinguishable from a background/noise rate, using
significance testing and bootstrap resampling.

## Decide-gate checklist

- **Outputs produced:** this notebook; analysis scaffold (`paper.md`,
  `gaps_and_friction.md`, `.gitignore`).
- **Results presented:** KG grounding counts (above); 15-strain data
  availability table; locked question.
- **QC gate:** KG release verdict checked → `ok` (0.1.0-alpha.6, 16/16
  schema asserts pass). All 15 researcher strains checked against
  `list_organisms` → all exist as genome entries; 7/15 have any expression
  data, 2/15 have nitrogen-treatment data, 5/15 have phosphorus-treatment
  data.
- **Decisions made this step:** redo-as-independent-walkthrough;
  axenic/uninfected/literal-starvation scope restriction; 15-strain scope
  restriction (all above, dated 2026-08-10).
- **Advance rationale:** question and scope are locked with the
  researcher's explicit input on four genuine judgment calls (hypotheses,
  confound handling, alt-N-source handling, strain scope); ready to move to
  step 2 (KG entries) to enumerate the actual publications/experiments that
  satisfy this scope.
