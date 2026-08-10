"""
Step 2 -- resolve the 92 researcher-supplied genes to locus tags.

Purpose: for each gene in the researcher's source list (identified by a
Cyanorak ortholog-group ID, not a bare gene name -- avoids paralog
ambiguity, see gene-identity.md Rule 2), find its locus tag in each of
the 4 strains that survived experiment selection (01_select_experiments.py):
MED4, MIT9312, MIT9313, NATL2A.

Uses genes_by_homolog_group(group_ids=["cyanorak:<CK_id>"], organisms=[...])
-- the ortholog-group route -- rather than name-based resolve_gene, since
the source list already carries the unambiguous Cyanorak group ID for all
but 2 rows (unnamed hypothetical genes with no gene name and no Cyanorak ID
in the source spreadsheet; these are recorded as unresolved).

Inputs: ../data/00_source_gene_list.csv (copied from the researcher's
  source spreadsheet extract).
Outputs: data/02_gene_locus_resolution.csv -- one row per (gene x strain),
  with the resolved locus_tag or a not_found/no_cyanorak_id flag.

Usage: uv run analyses/2026-08-10-np_starvation_expression_walkthrough/2_kg_selection/scripts/02_resolve_genes.py
"""

from pathlib import Path

import pandas as pd
from multiomics_explorer import genes_by_homolog_group, GraphConnection

DATA_DIR = Path(__file__).resolve().parent.parent / "data"
SOURCE_CSV = DATA_DIR / "00_source_gene_list.csv"

STRAINS = [
    "Prochlorococcus MED4",
    "Prochlorococcus MIT9312",
    "Prochlorococcus MIT9313",
    "Prochlorococcus NATL2A",
]


def main() -> None:
    genes = pd.read_csv(SOURCE_CSV)
    genes.columns = [c.strip() for c in genes.columns]
    genes["Cyanorak Ids"] = genes["Cyanorak Ids"].str.strip()
    print(f"Loaded {len(genes)} genes from source list")
    print(f"Genes with a Cyanorak ID: {genes['Cyanorak Ids'].notna().sum()}")
    print(f"Genes with no Cyanorak ID (unresolvable via ortholog group): "
          f"{genes['Cyanorak Ids'].isna().sum()}")

    rows = []
    with GraphConnection() as conn:
        for _, gene_row in genes.iterrows():
            ck_id = gene_row["Cyanorak Ids"]
            gene_name = gene_row["Gene_name"] if pd.notna(gene_row["Gene_name"]) else "(unnamed)"
            n_or_p = gene_row["N or P gene?"]

            if pd.isna(ck_id) or ck_id == "":
                for strain in STRAINS:
                    rows.append({
                        "n_or_p": n_or_p, "gene_name": gene_name, "cyanorak_id": None,
                        "organism_name": strain, "locus_tag": None,
                        "status": "no_cyanorak_id_in_source",
                    })
                continue

            group_id = f"cyanorak:{ck_id}"
            result = genes_by_homolog_group(
                group_ids=[group_id], organisms=STRAINS, verbose=False, conn=conn,
            )
            found_by_strain = {r["organism_name"]: r["locus_tag"] for r in result["results"]}
            not_found_group = group_id in (result.get("not_found") or [])

            for strain in STRAINS:
                if not_found_group:
                    status = "group_not_found_in_kg"
                    locus_tag = None
                elif strain in found_by_strain:
                    status = "resolved"
                    locus_tag = found_by_strain[strain]
                else:
                    status = "not_matched_in_strain"
                    locus_tag = None
                rows.append({
                    "n_or_p": n_or_p, "gene_name": gene_name, "cyanorak_id": ck_id,
                    "organism_name": strain, "locus_tag": locus_tag, "status": status,
                })

    df = pd.DataFrame(rows)
    out_path = DATA_DIR / "02_gene_locus_resolution.csv"
    df.to_csv(out_path, index=False)
    print(f"\nWrote {len(df)} rows to {out_path}")

    print("\nResolution status counts:")
    print(df["status"].value_counts())

    print("\nPer-strain resolved gene counts:")
    print(df[df["status"] == "resolved"].groupby("organism_name")["gene_name"].nunique())

    genes_with_any_resolution = df[df["status"] == "resolved"]["gene_name"].nunique()
    print(f"\nDistinct genes resolved to >=1 locus tag in >=1 in-scope strain: "
          f"{genes_with_any_resolution} / {genes['Gene_name'].notna().sum() + genes['Gene_name'].isna().sum()}")


if __name__ == "__main__":
    main()
