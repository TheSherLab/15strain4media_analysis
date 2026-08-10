"""
Step 3 -- materialize positive and negative control gene sets.

Purpose: build the actual gene lists (not just describe the rule) for:
  - positive controls: 4 canonical nitrogen-starvation-response genes
    (ntcA, glnA, amtB/amt1, ureA) and 4 canonical phosphorus-starvation-
    response genes (pstS, phoA, phoB, phoR) -- all confirmed present in
    the 60-gene evidence set from step 2.
  - negative/background: the full population of genes actually tested in
    each of the 4 nitrogen all_detected_genes experiments (the only ones
    with an unbiased gene population -- see 01_significance_criteria.py),
    excluding the 92 target genes and any gene whose product or gene name
    contains a nitrogen-metabolism keyword. This is nitrogen-only, per the
    step-3 framing decision that hypothesis 3 (noise/bootstrap check) does
    not run on phosphorus (no unbiased phosphorus background exists).

No pre-set background size: the full excluded-filtered population is kept
per experiment; the bootstrap (step 4/5) draws size-matched random samples
from it.

Inputs: ../../2_kg_selection/data/01_np_experiments.csv,
        ../../2_kg_selection/data/03_gene_evidence.csv (target gene list)
Outputs:
  data/02_positive_controls.csv
  data/02_negative_background_pool.csv -- one row per (background gene, experiment)

Usage: uv run analyses/2026-08-10-np_starvation_expression_walkthrough/3_analysis_framing/scripts/02_build_controls.py
"""

import re
from pathlib import Path

import pandas as pd
from multiomics_explorer import differential_expression_by_gene, GraphConnection

BASE = Path(__file__).resolve().parent.parent
DATA_DIR = BASE / "data"
STEP2_DATA = BASE.parent / "2_kg_selection" / "data"

POSITIVE_CONTROLS = {
    "N": ["ntcA", "glnA", "amtB/amt1", "ureA"],
    "P": ["pstS", "phoA", "phoB", "phoR"],
}

# Nitrogen-all_detected_genes experiments only (unbiased gene population).
N_UNBIASED_EXPERIMENTS = [
    "10.1101/2025.11.24.690089_growth_state_pro99lown_nutrient_starvation_med4_proteomics_axenic",
    "10.1101/2025.11.24.690089_growth_state_pro99lown_nutrient_starvation_med4_rnaseq_axenic",
    "10.1038/msb4100087_nitrogen_nitrogen_deprivation_med4_med4_microarray",
    "10.1038/msb4100087_nitrogen_nitrogen_deprivation_mit9313_mit9313_microarray",
]

N_KEYWORDS = [
    "nitrogen", "nitrate", "nitrite", "ammonium", "ammonia", "urea", "urease",
    "cyanate", "glutamine synthetase", "glutamate synthase", "gogat",
    "nitrogenase", "amino acid transport",
]
N_KEYWORD_PATTERN = re.compile("|".join(re.escape(k) for k in N_KEYWORDS), re.IGNORECASE)


def main() -> None:
    evidence = pd.read_csv(STEP2_DATA / "03_gene_evidence.csv")
    target_gene_names = set(evidence["gene_name"].dropna().unique())

    # -- positive controls: pull their locus tags from the resolved-genes table --
    resolution = pd.read_csv(STEP2_DATA / "02_gene_locus_resolution.csv")
    pos_rows = []
    for nutrient, genes in POSITIVE_CONTROLS.items():
        for g in genes:
            hits = resolution[(resolution["gene_name"] == g) & (resolution["status"] == "resolved")]
            for _, h in hits.iterrows():
                pos_rows.append({"nutrient": nutrient, "gene_name": g,
                                  "organism_name": h["organism_name"], "locus_tag": h["locus_tag"]})
    pos_df = pd.DataFrame(pos_rows)
    pos_df.to_csv(DATA_DIR / "02_positive_controls.csv", index=False)
    print(f"Positive controls: {len(pos_df)} (gene, strain) rows, "
          f"{pos_df.groupby('nutrient')['gene_name'].nunique().to_dict()}")

    # -- negative/background: full tested population per N-unbiased experiment,
    #    minus target genes, minus N-keyword products --
    exp_meta = pd.read_csv(STEP2_DATA / "01_np_experiments.csv").set_index("experiment_id")
    bg_rows = []
    with GraphConnection() as conn:
        for exp_id in N_UNBIASED_EXPERIMENTS:
            organism = exp_meta.loc[exp_id, "organism_name"].replace("Prochlorococcus ", "")
            result = differential_expression_by_gene(
                organism=organism, experiment_ids=[exp_id], limit=None, verbose=True, conn=conn,
            )
            de = pd.DataFrame(result["results"])
            de["gene_name_str"] = de["gene_name"].fillna("")
            de["product_str"] = de["product"].fillna("")

            is_target = de["gene_name_str"].isin(target_gene_names)
            is_n_keyword = (de["gene_name_str"].str.contains(N_KEYWORD_PATTERN)
                             | de["product_str"].str.contains(N_KEYWORD_PATTERN))
            background = de[~is_target & ~is_n_keyword]

            distinct_genes = background["locus_tag"].nunique()
            print(f"{exp_id}: {de['locus_tag'].nunique()} distinct tested genes -> "
                  f"{distinct_genes} in background pool after excluding "
                  f"{is_target.sum()} target-list rows and {is_n_keyword.sum()} N-keyword rows")

            bg_rows.append(background[["locus_tag", "gene_name", "product", "experiment_id",
                                        "expression_status", "log2fc", "padj"]])

    bg_df = pd.concat(bg_rows, ignore_index=True)
    bg_df.to_csv(DATA_DIR / "02_negative_background_pool.csv", index=False)
    print(f"\nWrote {len(bg_df)} rows to data/02_negative_background_pool.csv")
    print("Distinct background genes per experiment:")
    print(bg_df.groupby("experiment_id")["locus_tag"].nunique())


if __name__ == "__main__":
    main()
