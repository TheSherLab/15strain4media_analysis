"""
Step 2 -- a locus tag only means the gene exists in a strain's genome
(step 2 resolution). This checks each resolved (COG, strain, locus) row
against the 10 in-scope experiments' differential_expression_by_gene rows
to find which genes were actually measured.

Inputs: data/01_np_experiments.csv, data/02_gene_locus_resolution.csv
Outputs: data/03_gene_evidence.csv

Usage: uv run analyses/2026-09-06-cog_starvation_sensitivity_validation/2_kg_selection/scripts/03_check_evidence.py
"""

from pathlib import Path

import pandas as pd
from multiomics_explorer import differential_expression_by_gene, GraphConnection

BASE = Path(__file__).resolve().parent.parent
DATA_DIR = BASE / "data"


def main() -> None:
    experiments = pd.read_csv(DATA_DIR / "01_np_experiments.csv")
    resolution = pd.read_csv(DATA_DIR / "02_gene_locus_resolution.csv")
    resolved = resolution[resolution["status"] == "resolved"].copy()

    rows = []
    with GraphConnection() as conn:
        for organism_full, group in resolved.groupby("organism_name"):
            organism_short = organism_full.replace("Prochlorococcus ", "")
            exp_ids = experiments[
                experiments["organism_name"] == organism_full
            ]["experiment_id"].tolist()
            if not exp_ids:
                for _, g in group.iterrows():
                    rows.append({**g.to_dict(), "has_evidence": False, "n_de_rows": 0})
                continue

            loci = group["locus_tag"].unique().tolist()
            result = differential_expression_by_gene(
                organism=organism_short, locus_tags=loci, experiment_ids=exp_ids,
                limit=None, conn=conn,
            )
            de = pd.DataFrame(result["results"])
            counts = de["locus_tag"].value_counts() if len(de) else pd.Series(dtype=int)

            for _, g in group.iterrows():
                n = int(counts.get(g["locus_tag"], 0))
                rows.append({**g.to_dict(), "has_evidence": n > 0, "n_de_rows": n})

            print(f"{organism_short}: {len(loci)} distinct loci queried across "
                  f"{len(exp_ids)} experiments, {(counts > 0).sum()} with >=1 DE row")

    df = pd.DataFrame(rows)
    out_path = DATA_DIR / "03_gene_evidence.csv"
    df.to_csv(out_path, index=False)
    print(f"\nWrote {len(df)} rows to {out_path}")

    with_evidence = df[df["has_evidence"]]
    n_cogs_with_evidence = with_evidence["cog_number"].nunique()
    print(f"\n{n_cogs_with_evidence} of 41 COGs have >=1 DE row in >=1 in-scope experiment.")

    by_direction = (
        with_evidence.drop_duplicates(["cog_number", "direction"])
        .groupby("direction")["cog_number"].nunique()
    )
    print("\nCOGs with evidence, by Direction tag:")
    print(by_direction.to_string())

    per_strain = with_evidence.groupby("organism_name")["cog_number"].nunique().sort_values(ascending=False)
    print("\nPer-strain COGs with evidence (of resolved):")
    print(per_strain.to_string())


if __name__ == "__main__":
    main()
