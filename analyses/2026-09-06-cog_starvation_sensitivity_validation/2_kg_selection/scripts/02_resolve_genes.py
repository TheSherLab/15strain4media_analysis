"""
Step 2 -- resolve the 41 COGs to locus tags in the 15 study strains, using
the researcher's own COG-to-Cyanorak-ID mapping (verified in step 1: 41/41
COGs covered, 364 gene rows across the 15 strains, zero unmatched
groups/strains).

Inputs: data/00_source_normixed_cogs.csv, data/00_cog_to_ck_id_mapping.csv
Outputs: data/02_gene_locus_resolution.csv -- one row per (COG, strain)

Usage: uv run analyses/2026-09-06-cog_starvation_sensitivity_validation/2_kg_selection/scripts/02_resolve_genes.py
"""

from pathlib import Path

import pandas as pd
from multiomics_explorer import genes_by_homolog_group, GraphConnection

BASE = Path(__file__).resolve().parent.parent
DATA_DIR = BASE / "data"

STRAINS = [
    "MED4", "MIT9312", "MIT9313", "MIT1327", "MIT0604", "NATL2A", "MIT9515",
    "MIT9215", "AS9601", "PAC1", "MIT9202", "SB", "MIT9301", "MIT1314", "NATL1A",
]
STRAINS_FULL = [f"Prochlorococcus {s}" for s in STRAINS]


def main() -> None:
    source = pd.read_csv(DATA_DIR / "00_source_normixed_cogs.csv")
    mapping = pd.read_csv(DATA_DIR / "00_cog_to_ck_id_mapping.csv")
    merged = source.merge(mapping, left_on="COGs", right_on="cog_numbers")
    assert len(merged) == 41, f"expected 41 merged rows, got {len(merged)}"

    group_ids = [f"cyanorak:{ck}" for ck in merged["CK_ID"].unique()]

    with GraphConnection() as conn:
        result = genes_by_homolog_group(
            group_ids=group_ids, organisms=STRAINS_FULL, verbose=True, limit=None, conn=conn,
        )
        print(f"total_matching (gene rows): {result['total_matching']}")
        print(f"not_found_groups: {result['not_found_groups']}")
        print(f"not_matched_groups: {result['not_matched_groups']}")
        gene_rows = pd.DataFrame(result["results"])

    # gene_rows has columns: locus_tag, gene_name, product, organism_name, gene_category, group_id (=cyanorak:CK_...)
    gene_rows["CK_ID"] = gene_rows["group_id"].str.replace("cyanorak:", "", regex=False)

    # join back to COG numbers (a CK_ID can map to >1 COG, e.g. CK_00002384 -> COG4252 + COG2114).
    # A (CK_ID, strain) pair can also match >1 gene -- a genuine in-genome gene
    # duplication (verified: e.g. MED4 has 2 separate carbamoyltransferase-family
    # loci both assigned to CK_00002335), not a data artifact. Keep every match
    # as its own row rather than collapsing to one (Rule 2: locus tags, not gene
    # names -- paralogs are real biology).
    rows = []
    for strain_full in STRAINS_FULL:
        strain_genes = gene_rows[gene_rows["organism_name"] == strain_full]
        for _, cog_row in merged.iterrows():
            ck = cog_row["CK_ID"]
            matches = strain_genes[strain_genes["CK_ID"] == ck]
            if len(matches):
                for _, g in matches.iterrows():
                    rows.append({
                        "cog_number": cog_row["COGs"], "direction": cog_row["Direction"],
                        "protein_name": cog_row["Protien"], "CK_ID": ck, "organism_name": strain_full,
                        "locus_tag": g["locus_tag"], "gene_name": g["gene_name"], "product": g["product"],
                        "status": "resolved",
                    })
            else:
                rows.append({
                    "cog_number": cog_row["COGs"], "direction": cog_row["Direction"],
                    "protein_name": cog_row["Protien"], "CK_ID": ck, "organism_name": strain_full,
                    "locus_tag": None, "gene_name": None, "product": None,
                    "status": "not_matched_in_strain",
                })

    df = pd.DataFrame(rows)
    out_path = DATA_DIR / "02_gene_locus_resolution.csv"
    df.to_csv(out_path, index=False)
    print(f"\nWrote {len(df)} rows to {out_path}")

    # per-strain resolved-COG counts (distinct COG numbers with >=1 locus tag in that strain)
    per_strain = (
        df[df["status"] == "resolved"]
        .groupby("organism_name")["cog_number"].nunique()
        .sort_values(ascending=False)
    )
    print("\nPer-strain resolved COG counts (of 41):")
    print(per_strain.to_string())

    # multi-gene-per-strain check (a CK_ID matching >1 gene in the same strain)
    dupe_check = (
        gene_rows.groupby(["CK_ID", "organism_name"]).size().reset_index(name="n")
    )
    dupes = dupe_check[dupe_check["n"] > 1]
    if len(dupes):
        print("\nCK_ID with >1 gene in the same strain (genuine in-genome gene "
              "duplication -- both kept as separate rows, not collapsed):")
        print(dupes.to_string(index=False))
    else:
        print("\nNo CK_ID resolves to >1 gene within the same strain.")


if __name__ == "__main__":
    main()
