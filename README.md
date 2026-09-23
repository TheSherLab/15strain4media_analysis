# 15 Prochlorococcus Strains Across Four Media

This repository contains the data, exploratory analyses, figures, and
knowledge-graph (KG) analyses supporting the paper **"15 strains 4 media"**.
The study compares growth and starvation-related phenotypes across 15
*Prochlorococcus* strains grown in four media conditions.

The repository is intended to keep the paper's analysis materials together and
to make the computational work traceable. It includes the original Excel
datasets, analysis notebooks, genome annotations, generated analysis artifacts,
and two follow-up analyses using the laboratory's multi-omics KG.

## Repository contents

| Path | Contents |
|---|---|
| [`Datasets/`](Datasets/) | The three source Excel workbooks used for the paper (`Dataset 1.xlsx`, `Dataset 2.xlsx`, and `Dataset 3.xlsx`). |
| [`analyses/Analysis and figures/`](analyses/Analysis%20and%20figures/) | Jupyter notebooks for growth rates, figures, nutrient profiles, and supplementary figures. |
| [`analyses/Analysis and figures/gffs/`](analyses/Analysis%20and%20figures/gffs/) | GFF/GenBank genome annotations for the study strains. |
| [`analyses/KG analysis/`](analyses/KG%20analysis/) | Reproducible KG-backed analyses, staged data, scripts, figures, and research notes. |
| [`analyses/README.md`](analyses/README.md) | Notes on the analysis artifact structure and six-step research workflow. |

## Study scope

The paper covers these 15 *Prochlorococcus* strains:

`MED4`, `MIT9312`, `MIT9313`, `MIT1327`, `MIT0604`, `NATL2A`, `MIT9515`,
`MIT9215`, `AS9601`, `PAC1`, `MIT9202`, `SB`, `MIT9301`, `MIT1314`, and
`NATL1A`.

The KG analyses use a narrower evidence scope where the KG has matching
experiments: axenic, uninfected *Prochlorococcus* under literal nitrogen or
phosphorus starvation/limitation compared with nutrient-replete conditions.
The final scope contains 10 usable experiments across four strains and six
publications, using RNA-seq, proteomics, and microarray evidence.

## KG analyses

### Nitrogen and phosphorus acquisition genes

[`2026-08-10-np_starvation_expression_walkthrough/`](analyses/KG%20analysis/2026-08-10-np_starvation_expression_walkthrough/)

This analysis tests whether a 92-gene nitrogen/phosphorus acquisition list
responds to matched or non-matched nutrient starvation. It includes gene
identity checks, paralog audits, experiment selection, significance criteria,
hit-rate calculations, bootstrap tests, figures, and evaluation notes.

The final evidence set contains 60 target genes with expression evidence: 29
nitrogen-annotated genes and 31 phosphorus-annotated genes. Important limits of
the evidence, including filtered source tables and the lack of an unbiased
phosphorus background population, are recorded in the analysis notebooks.

### COG starvation-sensitivity validation

[`2026-09-06-cog_starvation_sensitivity_validation/`](analyses/KG%20analysis/2026-09-06-cog_starvation_sensitivity_validation/)

This analysis tests whether the 41 COGs identified by the paper's
cross-strain copy-number analysis also show an in-vivo expression response to
nitrogen or phosphorus starvation. Gene resolution was rebuilt from the full
per-strain COG membership rather than a single representative locus.

Headline results recorded in the analysis:

- Nitrogen: 15.5% significant tests for all 41 COGs, compared with a 19.1%
  nitrogen background; bootstrap p = 0.93.
- The nitrogen-tagged subset reached 23.2%, but was not distinguishable from
  its background (bootstrap p = 0.44).
- Phosphorus: 3.8% significant tests, reported descriptively because a valid
  phosphorus background was not available.
- Sensitivity checks gave the same overall conclusion after removing the
  largest COGs or collapsing each COG to one value per experiment.

## Reproducing the analyses

The notebooks in `analyses/Analysis and figures/` contain the paper's original
data processing and figure workflows. The KG analyses contain their own
`notebook.md` research records, scripts, staged CSV data, figures, and methods.
Start with the notebook in each analysis directory, then follow the scripts in
the numbered subdirectories.

The KG-backed work requires access to the laboratory multi-omics KG and the
corresponding Python environment. The KG is used for experiment metadata,
gene/locus resolution, orthology, and differential-expression evidence; the
repository records the extracted data and analysis decisions used for the
reported results.

The project declares its Python dependencies in [`pyproject.toml`](pyproject.toml).
For a KG-enabled environment, install the dependencies with:

```bash
uv sync
```

Do not commit credentials or private KG connection details. They should be
configured locally according to the laboratory's KG connection instructions.

## Manuscript files

The current manuscript revision is **15 strains 4 media paper revised version
92026**. The working paper folder contains the manuscript and supplemental
PDFs, editable Word files, figures, and response-to-reviewers materials. The
GitHub repository contains the computational materials and source datasets;
the manuscript working files are maintained separately from this repository.

## Data and interpretation notes

- Gene names can map to multiple loci or paralogs. The KG analyses therefore
  record locus-level resolution and explicitly audit paralogs.
- A locus present in a genome is not automatically evidence that it was
  measured in an experiment.
- KG results are tagged and documented with their evidence source, analysis
  decisions, and known gaps in the relevant `notebook.md` and
  `gaps_and_friction.md` files.
- The KG analyses are validation analyses of the paper's strain-level and
  COG-level findings; they do not replace the paper's original experimental
  datasets or statistical workflows.

## License and citation

Please cite the associated paper when reusing the datasets, notebooks, or
figures. Add the final publication citation here when the paper is published.
