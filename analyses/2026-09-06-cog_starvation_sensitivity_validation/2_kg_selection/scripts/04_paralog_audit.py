"""
Step 2 (reopened 2026-09-09) -- member-count audit of the researcher's
COG -> Cyanorak-ID mapping, mirroring the prior walkthrough analysis's
`2_kg_selection/scripts/05_paralog_audit.py`.

For every supplied Cyanorak group ID, count how many member genes it has
in each of the 15 study strains. A group with >1 member in a strain is a
family, not a single gene -- resolving it to "the ortholog" for that
strain is a guess unless the members are a genuine in-genome duplication
of the same gene. `02_resolve_genes.py` here keeps every member as its own
row (it never silently picks one, unlike the bug the walkthrough hit), so
this audit's job is to CONFIRM each multi-member case is a real
duplication, not to unbreak a picker.

Inputs:  data/00_cog_to_ck_id_mapping.csv, data/00_source_normixed_cogs.csv
Outputs: data/05_paralog_audit.csv

Usage: uv run analyses/2026-09-06-cog_starvation_sensitivity_validation/2_kg_selection/scripts/04_paralog_audit.py
"""

from pathlib import Path

import pandas as pd
from multiomics_explorer import genes_by_homolog_group, GraphConnection

DATA_DIR = Path(__file__).resolve().parent.parent / "data"

STRAINS = [
    "MED4", "MIT9312", "MIT9313", "MIT1327", "MIT0604", "NATL2A", "MIT9515",
    "MIT9215", "AS9601", "PAC1", "MIT9202", "SB", "MIT9301", "MIT1314", "NATL1A",
]
STRAINS_FULL = [f"Prochlorococcus {s}" for s in STRAINS]


def main() -> None:
    mapping = pd.read_csv(DATA_DIR / "00_cog_to_ck_id_mapping.csv")
    source = pd.read_csv(DATA_DIR / "00_source_normixed_cogs.csv")
    cog_meta = source.set_index("COGs")[["Protien", "Direction"]].to_dict("index")

    ck_to_cogs: dict[str, list[str]] = {}
    for _, r in mapping.iterrows():
        ck_to_cogs.setdefault(r["CK_ID"], []).append(r["cog_numbers"])

    rows = []
    with GraphConnection() as conn:
        for ck, cogs in ck_to_cogs.items():
            result = genes_by_homolog_group(
                group_ids=[f"cyanorak:{ck}"], organisms=STRAINS_FULL,
                verbose=True, limit=None, conn=conn,
            )
            members = pd.DataFrame(result["results"])
            for strain_full in STRAINS_FULL:
                if len(members):
                    hit = members[members["organism_name"] == strain_full]
                else:
                    hit = members
                loci = sorted(hit["locus_tag"]) if len(hit) else []
                products = sorted(set(hit["product"].dropna())) if len(hit) else []
                rows.append({
                    "CK_ID": ck,
                    "cog_numbers": "|".join(cogs),
                    "protein_name": "|".join(str(cog_meta.get(c, {}).get("Protien", "")) for c in cogs),
                    "direction": "|".join(str(cog_meta.get(c, {}).get("Direction", "")) for c in cogs),
                    "organism_name": strain_full,
                    "n_members_in_strain": len(loci),
                    "member_loci": "; ".join(loci),
                    "member_products": " :: ".join(products),
                    "multi_member": len(loci) > 1,
                })

    df = pd.DataFrame(rows)
    out = DATA_DIR / "05_paralog_audit.csv"
    df.to_csv(out, index=False)
    print(f"Wrote {len(df)} rows to {out}\n")

    multi = df[df["multi_member"]]
    print(f"{multi['CK_ID'].nunique()} Cyanorak IDs have >1 member in >=1 study strain, "
          f"across {len(multi)} (CK_ID, strain) pairs:\n")
    if len(multi):
        print(multi[["CK_ID", "cog_numbers", "protein_name", "organism_name",
                     "n_members_in_strain", "member_loci", "member_products"]].to_string(index=False))


if __name__ == "__main__":
    main()
