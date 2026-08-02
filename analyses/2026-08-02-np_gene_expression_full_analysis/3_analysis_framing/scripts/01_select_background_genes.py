"""Select the negative-control ("noise check") gene set for hypothesis 3:
100 genes per strain, not annotated as nitrogen/phosphorus-related, drawn
from the genes actually tested (any expression status) in this analysis's
12 N/P starvation experiments -- so the comparison is apples-to-apples
with the 54 target genes (same experiments, same detection sensitivity).

Method:
1. Pull every (gene x experiment x timepoint) DE row, significant or not
   (significant_only=False), for each strain's N/P experiment set from
   step 2 -- this is the tested-gene universe per strain.
2. Exclude the 54 target loci for that strain (data/03_gene_loci.csv from
   step 2).
3. Exclude any gene whose product or gene_name text matches an N/P
   metabolism keyword (nitrogen, ammonium, urea, nitrate, nitrite,
   cyanate, phosphate, phosphonate, "phospho", glutamine/glutamate
   synthase/synthetase, starvation) -- deliberately broad (over-excludes
   some unrelated "phospho-" central-metabolism genes) since the goal is
   a clean "definitely not N/P-acquisition-related" pool, not maximum
   recall.
4. Randomly sample up to 100 remaining genes per strain (seed=42, for
   reproducibility).

Inputs:  ../../2_kg_selection/data/02_np_experiments.csv,
         ../../2_kg_selection/data/03_gene_loci.csv
Outputs: data/01_background_genes.csv (up to 100 genes x 5 strains)

Usage:
  uv run python analyses/2026-08-02-np_gene_expression_full_analysis/3_analysis_framing/scripts/01_select_background_genes.py
"""
from pathlib import Path

import pandas as pd
from multiomics_explorer import differential_expression_by_gene

STEP_DIR = Path(__file__).resolve().parents[1]
DATA_DIR = STEP_DIR / "data"
KG_SELECTION_DATA = STEP_DIR.parent / "2_kg_selection" / "data"

N_PER_STRAIN = 100
SEED = 42

EXCLUDE_KEYWORDS = [
    "nitrogen", "ammoni", "urea", "nitrate", "nitrite", "cyanat", "cyanase",
    "phosphat", "phosphon", "phospho", "glutamine synthetase",
    "glutamate synthase", "starv",
]


def is_np_related(product: str, gene_name: str) -> bool:
    text = f"{product or ''} {gene_name or ''}".lower()
    return any(kw in text for kw in EXCLUDE_KEYWORDS)


def main() -> None:
    experiments = pd.read_csv(KG_SELECTION_DATA / "02_np_experiments.csv")
    gene_loci = pd.read_csv(KG_SELECTION_DATA / "03_gene_loci.csv")

    records = []
    for organism, exp_sub in experiments.groupby("organism_name"):
        exp_ids = exp_sub["experiment_id"].unique().tolist()
        result = differential_expression_by_gene(
            organism=organism,
            experiment_ids=exp_ids,
            significant_only=False,
            verbose=True,
            limit=None,
        )
        pool = pd.DataFrame(result["results"])
        pool = pool.drop_duplicates(subset="locus_tag")[
            ["locus_tag", "gene_name", "product", "gene_category"]
        ]

        target_loci = set(gene_loci.loc[gene_loci["organism_name"] == organism, "locus_tag"])
        pool = pool[~pool["locus_tag"].isin(target_loci)]

        not_related = pool[
            ~pool.apply(lambda r: is_np_related(r["product"], r["gene_name"]), axis=1)
        ]

        n = min(N_PER_STRAIN, len(not_related))
        sample = not_related.sample(n=n, random_state=SEED).sort_values("locus_tag")

        print(
            f"{organism}: tested pool {len(pool) + len(target_loci)} -> "
            f"{len(pool)} after excluding {len(target_loci)} target loci -> "
            f"{len(not_related)} after N/P-keyword exclusion -> sampled {len(sample)}"
        )

        for _, r in sample.iterrows():
            records.append(
                {
                    "organism_name": organism,
                    "locus_tag": r["locus_tag"],
                    "gene_name": r["gene_name"],
                    "product": r["product"],
                    "gene_category": r["gene_category"],
                }
            )

    out = pd.DataFrame(records)
    out_path = DATA_DIR / "01_background_genes.csv"
    out.to_csv(out_path, index=False)
    print(f"\nWrote {len(out)} background genes to {out_path}")
    print(out.groupby("organism_name").size().to_string())


if __name__ == "__main__":
    main()
