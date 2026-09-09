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

Multi-member fix (2026-09-08): the Cyanorak ID is NOT always a 1:1 ortholog
group. `05_paralog_audit.py` found 4 target groups with >1 member in an
in-scope strain -- worst is phoE (CK_00002330), a 3-6 member "outer membrane
porin" family per strain. The first version of this script built a dict
keyed by organism_name, so the LAST member silently won; for phoE it picked
a non-pho-operon porin in MED4/MIT9312/MIT9313/NATL2A. Now: a group with >1
member in a strain is resolved via MANUAL_LOCUS_OVERRIDE (curated by genomic
synteny to the pho/N operon -- see gaps_and_friction.md 2026-09-08); with no
override it is flagged `ambiguous_multi_member` and left unresolved rather
than guessed.

Inputs: ../data/00_source_gene_list.csv (copied from the researcher's
  source spreadsheet extract).
Outputs: data/02_gene_locus_resolution.csv -- one row per (gene x strain),
  with the resolved locus_tag or a not_found/ambiguous/no_cyanorak_id flag.

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

# Curated picks for Cyanorak groups with >1 member in a strain
# (05_paralog_audit.py). Key: (cyanorak_id, organism_name) -> locus_tag.
# Basis: genomic synteny to the pho gene island (phoB-phoR-phoA-porin-pstS),
# confirmed with gene_neighbors + strong P-induction of the chosen locus
# (e.g. phoE PMM0709: log2FC 117 at Martiny MED4 48h). Full rationale in
# gaps_and_friction.md (2026-09-08).
MANUAL_LOCUS_OVERRIDE = {
    ("CK_00002330", "Prochlorococcus MED4"): "PMM0709",        # phoE -- between phoA (PMM0708) & pstS (PMM0710)
    ("CK_00002330", "Prochlorococcus MIT9312"): "PMT9312_0721",  # phoE -- in phoA (0720)/pstS (0722) cluster
    ("CK_00002330", "Prochlorococcus MIT9313"): "PMT0998",      # phoE -- next to phoB (PMT0994)/phoR/pstS
    ("CK_00002330", "Prochlorococcus NATL2A"): "PMN2A_0440",    # phoE -- between phoA (PMN2A_0439) & pstS (PMN2A_0441)
    ("CK_00003432", "Prochlorococcus MED4"): "PMM0715",         # unkP2 -- pho island (arsR/arsB neighbours), sig_up Martiny 48h
    ("CK_00043821", "Prochlorococcus MIT9313"): "PMT0993",      # pstS -- adjacent to phoB (PMT0994); PMT0508 is the non-island paralog
    # PMM719 / CK_00044628 / NATL2A has 2 members (PMN2A_0477, PMN2A_0585);
    # neither has any DE row in an in-scope phosphorus experiment, so the
    # choice is immaterial to every downstream result. Left to fall through
    # to the ambiguous flag rather than guessed.
}

# The Cyanorak ID in the source spreadsheet points at the WRONG group.
# Found by 06_name_vs_cyanorak_crosscheck.py (name route vs Cyanorak route).
CYANORAK_ID_OVERRIDE = {
    # urtA was given CK_00008074, which is urtE's group ("ATPase component
    # UrtE"). Both the urtA and urtE spreadsheet rows carried it, so both
    # resolved to urtE's locus -- the analysis had urtE twice and no urtA.
    # CK_00000076 is the real urtA group ("substrate binding component"),
    # clean 1:1 across all 4 strains.
    "urtA": "CK_00000076",
}

# Cyanorak ID resolves in SOME strains but misses a real ortholog in another.
# Applied after the group lookup. Key: (gene_name, organism_name) -> locus.
FORCE_LOCUS = {
    # ptrA's spreadsheet ID CK_00056804 is a MED4 singleton group (1 member,
    # PMM0718). NATL2A's ptrA (PMN2A_0435, gene_name "ptrA", in the pho
    # island next to phoB) is in a different Cyanorak group (CK_00001606) and
    # was flagged not_matched. MIT9312/MIT9313 genuinely lack the pho-island
    # Crp regulator (checked by synteny -- no Crp gene in either pho island).
    ("ptrA", "Prochlorococcus NATL2A"): "PMN2A_0435",
}


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
            gene_name = gene_row["Gene_name"] if pd.notna(gene_row["Gene_name"]) else "(unnamed)"
            n_or_p = gene_row["N or P gene?"]
            ck_id = CYANORAK_ID_OVERRIDE.get(gene_name, gene_row["Cyanorak Ids"])

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
                group_ids=[group_id], organisms=STRAINS, verbose=False, conn=conn, limit=None,
            )
            # list of loci per strain -- NOT a dict (a strain can have >1 member)
            members_by_strain: dict[str, list[str]] = {s: [] for s in STRAINS}
            for r in result["results"]:
                members_by_strain.setdefault(r["organism_name"], []).append(r["locus_tag"])
            not_found_group = group_id in (result.get("not_found") or [])

            for strain in STRAINS:
                members = members_by_strain.get(strain, [])
                override = MANUAL_LOCUS_OVERRIDE.get((ck_id, strain))
                forced = FORCE_LOCUS.get((gene_name, strain))
                if forced is not None:
                    status, locus_tag = "resolved", forced
                elif not_found_group:
                    status, locus_tag = "group_not_found_in_kg", None
                elif not members:
                    status, locus_tag = "not_matched_in_strain", None
                elif len(members) == 1:
                    status, locus_tag = "resolved", members[0]
                elif override is not None:
                    assert override in members, f"{override} not in {ck_id} members for {strain}: {members}"
                    status, locus_tag = "resolved", override
                else:
                    status, locus_tag = "ambiguous_multi_member", None
                rows.append({
                    "n_or_p": n_or_p, "gene_name": gene_name, "cyanorak_id": ck_id,
                    "organism_name": strain, "locus_tag": locus_tag, "status": status,
                    "n_group_members_in_strain": len(members),
                    "all_members_in_strain": "; ".join(members),
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
