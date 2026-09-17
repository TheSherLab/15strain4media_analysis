"""
Step 3 -- build the nitrogen negative/background pool for the bootstrap
and Fisher's exact checks: every gene actually tested in each of the 4
unfiltered ("all_detected_genes") nitrogen experiments, at the same
single chosen starvation timepoint the prior walkthrough analysis locked
(reused directly -- no new timepoint judgment call here), excluding:

  1. This analysis's own 41-COG target list (by locus tag, per strain --
     not by gene_name text, since several of these genes have no
     gene_name populated in the KG; locus tag is unambiguous).
  2. The prior walkthrough analysis's 92-gene N/P-acquisition list (by
     locus tag, per strain) -- co-defined with the researcher: those
     genes are also already hypothesized to respond to starvation, so
     leaving them in the "random" pool would inflate the null baseline.
  3. Nitrogen-keyword product/gene-name text matches (same keyword list
     the prior analysis used), to catch other nitrogen-associated genes
     not on either curated list.

Inputs:
  ../../2_kg_selection/data/02_gene_locus_resolution.csv (this analysis's 41-COG resolution)
  ../../../2026-08-10-np_starvation_expression_walkthrough/2_kg_selection/data/02_gene_locus_resolution.csv (92-gene resolution)
  ../../2_kg_selection/data/01_np_experiments.csv
Outputs: data/02_negative_background_pool.csv

Usage: uv run analyses/2026-09-06-cog_starvation_sensitivity_validation/3_analysis_framing/scripts/02_build_background_pool.py
"""

import re
from pathlib import Path

import pandas as pd
from multiomics_explorer import differential_expression_by_gene, GraphConnection

BASE = Path(__file__).resolve().parent.parent
DATA_DIR = BASE / "data"
STEP2_DATA = BASE.parent / "2_kg_selection" / "data"
PRIOR_RESOLUTION_CSV = (
    BASE.parent.parent
    / "2026-08-10-np_starvation_expression_walkthrough"
    / "2_kg_selection" / "data" / "02_gene_locus_resolution.csv"
)

N_UNBIASED_EXPERIMENTS = [
    "10.1101/2025.11.24.690089_growth_state_pro99lown_nutrient_starvation_med4_proteomics_axenic",
    "10.1101/2025.11.24.690089_growth_state_pro99lown_nutrient_starvation_med4_rnaseq_axenic",
    "10.1038/msb4100087_nitrogen_nitrogen_deprivation_med4_med4_microarray",
    "10.1038/msb4100087_nitrogen_nitrogen_deprivation_mit9313_mit9313_microarray",
]

# Same chosen starvation timepoint per experiment as the prior analysis
# (5_analyze/notebook.md there; reused directly, no new judgment call).
N_STARVATION_TIMEPOINT = {
    "10.1101/2025.11.24.690089_growth_state_pro99lown_nutrient_starvation_med4_proteomics_axenic": "day 14",
    "10.1101/2025.11.24.690089_growth_state_pro99lown_nutrient_starvation_med4_rnaseq_axenic": None,
    "10.1038/msb4100087_nitrogen_nitrogen_deprivation_med4_med4_microarray": "12h",
    "10.1038/msb4100087_nitrogen_nitrogen_deprivation_mit9313_mit9313_microarray": "12h",
}

N_KEYWORDS = [
    "nitrogen", "nitrate", "nitrite", "ammonium", "ammonia", "urea", "urease",
    "cyanate", "glutamine synthetase", "glutamate synthase", "gogat",
    "nitrogenase", "amino acid transport",
]
N_KEYWORD_PATTERN = re.compile("|".join(re.escape(k) for k in N_KEYWORDS), re.IGNORECASE)


def main() -> None:
    cog_resolution = pd.read_csv(STEP2_DATA / "02_gene_locus_resolution.csv")
    cog_loci_resolved = cog_resolution[cog_resolution["status"] == "resolved"]
    cog_exclude = set(
        zip(cog_loci_resolved["organism_name"], cog_loci_resolved["locus_tag"])
    )

    prior_resolution = pd.read_csv(PRIOR_RESOLUTION_CSV)
    prior_resolved = prior_resolution[prior_resolution["status"] == "resolved"]
    prior_exclude = set(
        zip(prior_resolved["organism_name"], prior_resolved["locus_tag"])
    )

    exclude_loci = cog_exclude | prior_exclude
    print(f"Exclusion set: {len(cog_exclude)} (organism, locus) pairs from the 41-COG "
          f"list, {len(prior_exclude)} from the 92-gene list, {len(exclude_loci)} combined (union).")

    exp_meta = pd.read_csv(STEP2_DATA / "01_np_experiments.csv").set_index("experiment_id")
    bg_rows = []
    with GraphConnection() as conn:
        for exp_id in N_UNBIASED_EXPERIMENTS:
            organism_full = exp_meta.loc[exp_id, "organism_name"]
            organism = organism_full.replace("Prochlorococcus ", "")
            result = differential_expression_by_gene(
                organism=organism, experiment_ids=[exp_id], limit=None, verbose=True, conn=conn,
            )
            de = pd.DataFrame(result["results"])

            chosen_tp = N_STARVATION_TIMEPOINT[exp_id]
            if chosen_tp is not None:
                before = len(de)
                de = de[de["timepoint"] == chosen_tp]
                print(f"{exp_id}: restricted to timepoint '{chosen_tp}' -> {len(de)} of {before} rows")

            de["gene_name_str"] = de["gene_name"].fillna("")
            de["product_str"] = de["product"].fillna("")

            is_target_locus = de["locus_tag"].apply(lambda lt: (organism_full, lt) in exclude_loci)
            is_n_keyword = (de["gene_name_str"].str.contains(N_KEYWORD_PATTERN)
                             | de["product_str"].str.contains(N_KEYWORD_PATTERN))
            background = de[~is_target_locus & ~is_n_keyword]

            distinct_genes = background["locus_tag"].nunique()
            print(f"{exp_id}: {de['locus_tag'].nunique()} distinct tested genes (at chosen timepoint) -> "
                  f"{distinct_genes} in background pool after excluding "
                  f"{is_target_locus.sum()} target-locus rows (41-COG + 92-gene lists) and "
                  f"{is_n_keyword.sum()} N-keyword rows")

            bg_rows.append(background[["locus_tag", "gene_name", "product", "experiment_id",
                                        "expression_status", "log2fc", "padj"]])

    bg_df = pd.concat(bg_rows, ignore_index=True)
    out_path = DATA_DIR / "02_negative_background_pool.csv"
    bg_df.to_csv(out_path, index=False)
    print(f"\nWrote {len(bg_df)} rows to {out_path}")
    print("Distinct background genes per experiment:")
    print(bg_df.groupby("experiment_id")["locus_tag"].nunique())


if __name__ == "__main__":
    main()
