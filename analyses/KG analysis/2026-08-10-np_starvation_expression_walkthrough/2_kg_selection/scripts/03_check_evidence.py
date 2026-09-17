"""
Step 2 -- check which resolved genes actually have expression evidence.

Purpose: a locus tag (02_resolve_genes.py) only means the gene exists in
a strain's genome. This script checks, for each resolved (gene, strain)
pair, whether the gene has any differential-expression row (any
significance call, any direction) in the 10 in-scope experiments
(01_select_experiments.py) for that strain -- i.e. whether it was actually
measured, not just annotated.

Inputs: ../data/01_np_experiments.csv, ../data/02_gene_locus_resolution.csv
Outputs: data/03_gene_evidence.csv -- one row per resolved (gene, strain),
  with has_evidence (bool) and the count of DE rows found.

Usage: uv run analyses/2026-08-10-np_starvation_expression_walkthrough/2_kg_selection/scripts/03_check_evidence.py
"""

from pathlib import Path

import pandas as pd
from multiomics_explorer import differential_expression_by_gene, GraphConnection

DATA_DIR = Path(__file__).resolve().parent.parent / "data"


def main() -> None:
    experiments = pd.read_csv(DATA_DIR / "01_np_experiments.csv")
    included_exp = experiments[experiments["include"]]

    resolution = pd.read_csv(DATA_DIR / "02_gene_locus_resolution.csv")
    resolved = resolution[resolution["status"] == "resolved"].copy()
    print(f"Checking evidence for {len(resolved)} (gene, strain) pairs "
          f"({resolved['gene_name'].nunique()} distinct genes)")

    rows = []
    with GraphConnection() as conn:
        for strain, strain_resolved in resolved.groupby("organism_name"):
            strain_exp_ids = included_exp[included_exp["organism_name"] == strain]["experiment_id"].tolist()
            if not strain_exp_ids:
                for _, r in strain_resolved.iterrows():
                    rows.append({**r.to_dict(), "has_evidence": False, "de_row_count": 0,
                                 "note": "no in-scope experiments for this strain"})
                continue

            locus_tags = strain_resolved["locus_tag"].tolist()
            result = differential_expression_by_gene(
                organism=strain,
                locus_tags=locus_tags,
                experiment_ids=strain_exp_ids,
                verbose=False,
                limit=None,
                conn=conn,
            )
            counts = pd.Series([r["locus_tag"] for r in result["results"]]).value_counts()

            for _, r in strain_resolved.iterrows():
                n = int(counts.get(r["locus_tag"], 0))
                rows.append({**r.to_dict(), "has_evidence": n > 0, "de_row_count": n, "note": ""})

    df = pd.DataFrame(rows)
    out_path = DATA_DIR / "03_gene_evidence.csv"
    df.to_csv(out_path, index=False)
    print(f"\nWrote {len(df)} rows to {out_path}")

    print("\nEvidence summary per strain:")
    print(df.groupby("organism_name")["has_evidence"].agg(["sum", "count"]))

    any_evidence_genes = df[df["has_evidence"]]["gene_name"].nunique()
    total_genes = 92
    print(f"\nDistinct genes with actual DE evidence in >=1 in-scope experiment: "
          f"{any_evidence_genes} / {total_genes}")

    print("\nBy N/P annotation:")
    with_evidence = df[df["has_evidence"]].drop_duplicates("gene_name")
    print(with_evidence.groupby("n_or_p")["gene_name"].nunique())


if __name__ == "__main__":
    main()
