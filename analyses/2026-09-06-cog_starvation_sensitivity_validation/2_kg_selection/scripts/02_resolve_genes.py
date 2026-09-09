"""
Step 2 (reopened 2026-09-09) -- resolve the full per-strain COG membership
(`00_wholegenome_cog_members.csv`, from the researcher's whole-genome
annotation) to multiomics-KG locus tags, so expression can be pulled.

Each gene in the whole-genome file carries a Cyanorak cluster number. We
map (cluster, strain) -> KG locus tag(s) with a single batched
`genes_by_homolog_group` call over all 165 clusters x 13 strains, then
join back to the whole-genome gene rows.

  - 1 whole-genome gene, 1 KG gene for that (cluster, strain)  -> resolved
  - N whole-genome genes, N KG genes (in-genome duplication)   -> resolved, all kept
  - whole-genome gene, 0 KG genes for that (cluster, strain)   -> unmapped_in_kg
  - count mismatch                                             -> flagged, all KG loci kept

Output is driven by the whole-genome file (the researcher's COG assignment
is authoritative); KG is used only to attach an expression-queryable locus.

Inputs:  data/00_wholegenome_cog_members.csv
Outputs: data/02_gene_locus_resolution.csv  -- one row per (cog, strain, KG locus)

Usage: uv run analyses/2026-09-06-cog_starvation_sensitivity_validation/2_kg_selection/scripts/02_resolve_genes.py
"""

from pathlib import Path

import pandas as pd
from multiomics_explorer import genes_by_homolog_group, GraphConnection

BASE = Path(__file__).resolve().parent.parent
DATA_DIR = BASE / "data"


def main() -> None:
    members = pd.read_csv(DATA_DIR / "00_wholegenome_cog_members.csv")
    strains = sorted(members["strain"].unique())
    strains_full = [f"Prochlorococcus {s}" for s in strains]
    clusters = sorted(members["cluster_number"].dropna().unique())
    print(f"{len(members)} whole-genome gene rows, {len(clusters)} clusters, {len(strains)} strains")

    with GraphConnection() as conn:
        res = genes_by_homolog_group(
            group_ids=[f"cyanorak:{c}" for c in clusters],
            organisms=strains_full, verbose=True, limit=None, conn=conn,
        )
    kg = pd.DataFrame(res["results"])
    kg["cluster_number"] = kg["group_id"].str.replace("cyanorak:", "", regex=False)
    kg["strain"] = kg["organism_name"].str.replace("Prochlorococcus ", "", regex=False)
    print(f"KG returned {len(kg)} (cluster x gene x strain) rows; "
          f"not_found_groups={res.get('not_found_groups')}")

    # {(cluster, strain): [kg rows]}
    kg_by_key: dict[tuple, list[dict]] = {}
    for _, r in kg.iterrows():
        kg_by_key.setdefault((r["cluster_number"], r["strain"]), []).append(r.to_dict())

    rows = []
    for (cog, strain, cluster), grp in members.groupby(["cog_numbers", "strain", "cluster_number"], dropna=False):
        wg_rows = grp.to_dict("records")
        kg_hits = kg_by_key.get((cluster, strain), [])
        direction = wg_rows[0]["direction"]

        if not kg_hits:
            for w in wg_rows:
                rows.append({
                    "cog_number": cog, "direction": direction, "organism_name": f"Prochlorococcus {strain}",
                    "cluster_number": cluster, "wg_locus_name": w["locus_name"], "wg_id": w["ckpro_id"],
                    "locus_tag": None, "kg_gene_name": None, "kg_product": None,
                    "wg_product": w["product"], "status": "unmapped_in_kg",
                })
            continue

        status = "resolved" if len(kg_hits) == len(wg_rows) else "count_mismatch_kept_all"
        for h in kg_hits:
            rows.append({
                "cog_number": cog, "direction": direction, "organism_name": h["organism_name"],
                "cluster_number": cluster, "wg_locus_name": wg_rows[0]["locus_name"] if len(wg_rows) == 1 else "(multiple)",
                "wg_id": wg_rows[0]["ckpro_id"] if len(wg_rows) == 1 else "(multiple)",
                "locus_tag": h["locus_tag"], "kg_gene_name": h.get("gene_name"),
                "kg_product": h.get("product"), "wg_product": wg_rows[0]["product"],
                "status": status,
            })

    df = pd.DataFrame(rows).drop_duplicates(subset=["cog_number", "organism_name", "locus_tag"])
    out = DATA_DIR / "02_gene_locus_resolution.csv"
    df.to_csv(out, index=False)
    print(f"\nWrote {len(df)} rows to {out}")
    print("\nstatus counts:")
    print(df["status"].value_counts().to_string())

    res_df = df[df["status"] != "unmapped_in_kg"]
    print(f"\n{res_df['cog_number'].nunique()} / 41 COGs have >=1 KG-resolved locus in >=1 strain")
    ev = ["Prochlorococcus MED4", "Prochlorococcus MIT9313",
          "Prochlorococcus NATL2A", "Prochlorococcus MIT9312"]
    ev_df = res_df[res_df["organism_name"].isin(ev)]
    print(f"evidence strains: {len(ev_df)} (cog x strain x locus) rows, "
          f"{ev_df['cog_number'].nunique()} COGs, {ev_df['locus_tag'].nunique()} distinct loci")
    print("\nunmapped_in_kg detail:")
    um = df[df["status"] == "unmapped_in_kg"]
    if len(um):
        print(um.groupby(["cog_number"])["cluster_number"].apply(lambda s: sorted(set(s))).to_string())


if __name__ == "__main__":
    main()
