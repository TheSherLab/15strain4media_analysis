# N/P acquisition gene expression response under starvation — walkthrough analysis

## Question

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

## Background

The KG holds 12 Prochlorococcus nitrogen-treatment and 9 phosphorus-treatment
experiments. Applying the step-1 scope (axenic, uninfected, literal
starvation/limitation, restricted to the researcher's 15 study strains)
narrows this to 10 usable experiments across 4 strains (MED4, MIT9312,
MIT9313, NATL2A) and 6 publications: 5 nitrogen (Weissberg et al. 2025 —
proteomics and RNA-seq, axenic arm only; Read et al. 2017 — RNA-seq;
Tolonen et al. 2006 — microarray, MED4 and MIT9313) and 5 phosphorus (Lin,
Ding, Zeng 2015 — RNA-seq, uninfected arm only; Fuszard et al. 2012 —
proteomics, MIT9312 and NATL2A; Martiny, Coleman, Chisholm 2006 —
microarray, MED4 and MIT9313). Excluded: coculture- and
phage-infection-background arms, 2 strains not in the researcher's 15
(SS120), 4 alternate-N-source substrate-switch experiments, and 2
metabolomics-only experiments with no gene-level differential expression.

Of the researcher's 92-gene reference list, 91 carry a Cyanorak
ortholog-group ID; resolving each via `genes_by_homolog_group` against the
4 in-scope strains finds a locus tag for 61 genes in at least one strain
(MED4 53, NATL2A 48, MIT9312 44, MIT9313 42). A locus tag only means the
gene is present in that strain's genome, not that it was measured;
checking each resolved gene against the 10 in-scope experiments'
differential-expression rows narrows this to **60 genes (29
nitrogen-annotated, 31 phosphorus-annotated) with actual expression
evidence** in at least one in-scope experiment. The remaining genes have
no locus tag in any of the 4 strains — consistent with known
strain/ecotype-level gene content variation in Prochlorococcus rather than
a resolution gap (one exception, `phnW`, has no Cyanorak ID at all and
resolves only in *Pseudomonas putida* on direct name lookup — genuinely
absent from Prochlorococcus in this KG build). One strain, MIT9312,
resolves 44 genes but has evidence for only 2, because its sole in-scope
experiment's source table only reports genes that were already
significant in the original publication (see `gaps_and_friction.md`). See
`2_kg_selection/notebook.md` for the full experiment and gene resolution
tables.

## Methods

## Results

## Discussion

## References

1. Weissberg O, Aharonovich D, Sher D (2025). Transcriptomic and Proteomic
   Analysis Reveals Nitrogen Recycling as a Core Mechanism for
   Prochlorococcus Prolonged Survival. *bioRxiv*.
   https://doi.org/10.1101/2025.11.24.690089
2. Read RW, Berube PM, Biller SJ, Neveux I, Cubillos-Ruiz A, Chisholm SW,
   Grzymski JJ (2017). Nitrogen cost minimization is promoted by
   structural changes in the transcriptome of N-deprived Prochlorococcus
   cells. *The ISME Journal*. https://doi.org/10.1038/ismej.2017.88
3. Tolonen AC, Aach J, Lindell D, Johnson ZI, Rector T, Steen R, Church
   GM, Chisholm SW (2006). Global gene expression of Prochlorococcus
   ecotypes in response to changes in nitrogen availability. *Molecular
   Systems Biology*. https://doi.org/10.1038/msb4100087
4. Lin X, Ding H, Zeng Q (2015). Transcriptomic response during phage
   infection of a marine cyanobacterium under phosphorus-limited
   conditions. *Environmental Microbiology*.
   https://doi.org/10.1111/1462-2920.13104
5. Fuszard MA, Wright PC, Biggs CA (2012). Comparative quantitative
   proteomics of Prochlorococcus ecotypes to a decrease in environmental
   phosphate concentrations. *Aquatic Biosystems*.
   https://doi.org/10.1186/2046-9063-8-7
6. Martiny AC, Coleman ML, Chisholm SW (2006). Phosphate acquisition
   genes in Prochlorococcus ecotypes: Evidence for genome-wide adaptation.
   *PNAS*. https://doi.org/10.1073/pnas.0601301103
