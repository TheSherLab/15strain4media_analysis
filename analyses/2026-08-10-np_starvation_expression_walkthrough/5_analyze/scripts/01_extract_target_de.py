"""
Step 5 -- extract DE data for all 60 target genes across all 10 in-scope
experiments, restricted to each experiment's single chosen starvation
timepoint (step 5 timepoint table, see notebook.md).

For every (gene, experiment) pair, produces exactly one of 4 states:
  - "no_locus_in_strain": the gene has no locus tag in this experiment's
    strain (genuine strain-level gene absence, from step 2 resolution).
  - "no_data_at_timepoint": the gene has a locus tag but no DE row at the
    chosen timepoint in this experiment (not reported/tested there).
  - "significant_up" / "significant_down" / "not_significant": from the
    KG's own expression_status at the chosen timepoint.

Also tags each (gene, experiment) pair as "matched" (gene's N/P annotation
equals the experiment's nutrient) or "cross" (it doesn't), for hypotheses
1 and 2.

Inputs: ../../2_kg_selection/data/01_np_experiments.csv,
        ../../2_kg_selection/data/02_gene_locus_resolution.csv,
        ../../2_kg_selection/data/03_gene_evidence.csv
Outputs: data/01_target_gene_experiment_matrix.csv

Usage: uv run analyses/2026-08-10-np_starvation_expression_walkthrough/5_analyze/scripts/01_extract_target_de.py
"""

from pathlib import Path

import pandas as pd
from multiomics_explorer import differential_expression_by_gene, GraphConnection

BASE = Path(__file__).resolve().parent.parent
DATA_DIR = BASE / "data"
STEP2_DATA = BASE.parent / "2_kg_selection" / "data"

# Chosen starvation timepoint per experiment (step 5 timepoint table,
# notebook.md). None = single-timepoint experiment (no filter needed).
STARVATION_TIMEPOINT = {
    "10.1101/2025.11.24.690089_growth_state_pro99lown_nutrient_starvation_med4_proteomics_axenic": "day 14",
    "10.1101/2025.11.24.690089_growth_state_pro99lown_nutrient_starvation_med4_rnaseq_axenic": None,
    "10.1038/ismej.2017.88_nitrogen_stress_ndepleted_pro99_medium_med4_rnaseq": "24h",
    "10.1038/msb4100087_nitrogen_nitrogen_deprivation_med4_med4_microarray": "12h",
    "10.1038/msb4100087_nitrogen_nitrogen_deprivation_mit9313_mit9313_microarray": "12h",
    "10.1111/1462-2920.13104_phosphorus_plimited_natl2a_rnaseq_uninfected": "46h",
    "10.1186/2046-9063-8-7_pi_limitation_mit9312_itraq": None,
    "10.1186/2046-9063-8-7_pi_limitation_natl2a_itraq": None,
    "10.1073/pnas.0601301103_phosphorus_phosphate_starvation_med4_microarray": "48h",
    "10.1073/pnas.0601301103_phosphorus_phosphate_starvation_mit9313_microarray": "24h",  # flagged: see notebook.md
}


def main() -> None:
    experiments = pd.read_csv(STEP2_DATA / "01_np_experiments.csv")
    included = experiments[experiments["include"]].copy()
    included["organism_short"] = included["organism_name"].str.replace("Prochlorococcus ", "", regex=False)

    resolution = pd.read_csv(STEP2_DATA / "02_gene_locus_resolution.csv")
    evidence = pd.read_csv(STEP2_DATA / "03_gene_evidence.csv")
    target_genes = evidence[["n_or_p", "gene_name"]].drop_duplicates()

    rows = []
    with GraphConnection() as conn:
        for _, exp in included.iterrows():
            exp_id = exp["experiment_id"]
            organism_full = exp["organism_name"]
            organism_short = exp["organism_short"]
            chosen_tp = STARVATION_TIMEPOINT[exp_id]

            gene_locus = resolution[
                (resolution["organism_name"] == organism_full) & (resolution["status"] == "resolved")
            ].set_index("gene_name")["locus_tag"].to_dict()

            # List of (gene_name, n_or_p, locus_tag) -- NOT a dict keyed by locus_tag,
            # since two list entries can share a locus tag (e.g. urtA/urtE share a
            # Cyanorak ID in the source spreadsheet -- a genuine data-quality issue,
            # not a bug here; both must still appear as separate rows).
            gene_locus_pairs = []
            for _, g in target_genes.iterrows():
                gene_name = g["gene_name"]
                if gene_name not in gene_locus:
                    rows.append({
                        "gene_name": gene_name, "n_or_p": g["n_or_p"], "experiment_id": exp_id,
                        "nutrient": exp["nutrient"], "organism_name": organism_full,
                        "match_type": "matched" if g["n_or_p"] == exp["nutrient"][0].upper() else "cross",
                        "locus_tag": None, "status": "no_locus_in_strain",
                        "log2fc": None, "padj": None, "timepoint": chosen_tp,
                    })
                    continue
                gene_locus_pairs.append((gene_name, g["n_or_p"], gene_locus[gene_name]))

            distinct_loci = list({locus for _, _, locus in gene_locus_pairs})
            if distinct_loci:
                result = differential_expression_by_gene(
                    organism=organism_short, locus_tags=distinct_loci,
                    experiment_ids=[exp_id], limit=None, conn=conn,
                )
                de = pd.DataFrame(result["results"])
                if chosen_tp is not None and len(de):
                    de = de[de["timepoint"] == chosen_tp]
                de_by_locus = {r["locus_tag"]: r for _, r in de.iterrows()} if len(de) else {}

                for gene_name, n_or_p, locus in gene_locus_pairs:
                    match_type = "matched" if n_or_p == exp["nutrient"][0].upper() else "cross"
                    if locus in de_by_locus:
                        row = de_by_locus[locus]
                        rows.append({
                            "gene_name": gene_name, "n_or_p": n_or_p, "experiment_id": exp_id,
                            "nutrient": exp["nutrient"], "organism_name": organism_full,
                            "match_type": match_type, "locus_tag": locus,
                            "status": row["expression_status"], "log2fc": row["log2fc"],
                            "padj": row["padj"], "timepoint": chosen_tp,
                        })
                    else:
                        rows.append({
                            "gene_name": gene_name, "n_or_p": n_or_p, "experiment_id": exp_id,
                            "nutrient": exp["nutrient"], "organism_name": organism_full,
                            "match_type": match_type, "locus_tag": locus,
                            "status": "no_data_at_timepoint", "log2fc": None, "padj": None,
                            "timepoint": chosen_tp,
                        })

            print(f"{exp_id} ({organism_short}, timepoint={chosen_tp}): "
                  f"{len(gene_locus_pairs)} genes queried ({len(distinct_loci)} distinct loci)")

    df = pd.DataFrame(rows)
    out_path = DATA_DIR / "01_target_gene_experiment_matrix.csv"
    df.to_csv(out_path, index=False)
    print(f"\nWrote {len(df)} rows to {out_path}")
    print("\nStatus counts:")
    print(df["status"].value_counts())
    print("\nBy match_type x status:")
    print(pd.crosstab(df["match_type"], df["status"]))


if __name__ == "__main__":
    main()
