"""
Step 2 (redo, 2026-09-08) -- audit every target Cyanorak group for
multi-paralog membership in the 4 in-scope strains.

Motivation: 02_resolve_genes.py built `{organism_name: locus_tag}` from
genes_by_homolog_group results -- a dict keyed by strain, so when a
Cyanorak group has >1 member in a strain the LAST one silently wins.
This surfaced for phoE (CK_00002330), a 3-6 member "outer membrane porin"
family per strain: the resolver picked a non-pho-operon paralog in MED4,
MIT9312 and NATL2A. This script finds every group with the same problem.

Inputs: ../data/00_source_gene_list.csv
Outputs: ../data/05_paralog_audit.csv -- one row per (gene x strain) with
  the full member list where >1, plus a flag column.

Usage: uv run analyses/2026-08-10-np_starvation_expression_walkthrough/2_kg_selection/scripts/05_paralog_audit.py
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
    genes["Cyanorak Ids"] = genes["Cyanorak Ids"].astype("string").str.strip()

    rows = []
    with GraphConnection() as conn:
        for _, g in genes.iterrows():
            ck_id = g["Cyanorak Ids"]
            gene_name = g["Gene_name"] if pd.notna(g["Gene_name"]) else "(unnamed)"
            n_or_p = g["N or P gene?"]
            if pd.isna(ck_id) or ck_id == "":
                continue

            result = genes_by_homolog_group(
                group_ids=[f"cyanorak:{ck_id}"], organisms=STRAINS,
                verbose=True, conn=conn, limit=None,
            )
            by_strain: dict[str, list[dict]] = {s: [] for s in STRAINS}
            for r in result["results"]:
                by_strain.setdefault(r["organism_name"], []).append(r)

            for strain in STRAINS:
                members = by_strain.get(strain, [])
                rows.append({
                    "n_or_p": n_or_p,
                    "gene_name": gene_name,
                    "cyanorak_id": ck_id,
                    "organism_name": strain,
                    "n_members": len(members),
                    "member_loci": "; ".join(m["locus_tag"] for m in members),
                    "member_products": " | ".join(
                        f'{m["locus_tag"]}={m.get("product") or m.get("consensus_product") or "?"}'
                        for m in members
                    ),
                    "ambiguous": len(members) > 1,
                })

    df = pd.DataFrame(rows)
    out_path = DATA_DIR / "05_paralog_audit.csv"
    df.to_csv(out_path, index=False)
    print(f"Wrote {len(df)} rows to {out_path}")

    amb = df[df["ambiguous"]]
    print(f"\n{amb['gene_name'].nunique()} genes have >1 group member in >=1 in-scope strain:")
    for gene_name, sub in amb.groupby("gene_name"):
        strains = ", ".join(
            f'{s.replace("Prochlorococcus ", "")}(n={n})'
            for s, n in zip(sub["organism_name"], sub["n_members"])
        )
        print(f"  {gene_name:12} [{sub['cyanorak_id'].iloc[0]}]  {strains}")


if __name__ == "__main__":
    main()
