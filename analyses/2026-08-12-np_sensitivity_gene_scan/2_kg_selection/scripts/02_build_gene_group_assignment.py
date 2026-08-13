"""
Step 2 -- consolidate the raw gene + ortholog-edge pulls into one
per-gene assignment table: which ortholog group (if any) represents
this gene for cross-strain comparison, and which system supplied it.

Assignment rule (co-defined with the researcher):
  1. Cyanorak group (specificity_rank=0) if the gene has one -- primary.
  2. Else, Eggnog rank=1 (tightest/family-level) group if it has one --
     fallback, used only for genes Cyanorak didn't curate.
  3. Else, "true orphan" -- no group in either system; excluded from the
     per-ortholog-group scan (step 5), tracked separately as a per-strain
     count (see 03_orphan_summary.csv).

Inputs: data/01_genes_all_strains.csv, data/02_ortholog_edges_all_strains.csv
Outputs: data/03_gene_group_assignment.csv (one row per gene)
         data/04_group_landscape_summary.csv (one row per source/rank)
         data/05_orphan_summary_by_strain.csv

Usage: uv run analyses/2026-08-12-np_sensitivity_gene_scan/2_kg_selection/scripts/02_build_gene_group_assignment.py
"""

from pathlib import Path

import pandas as pd

BASE = Path(__file__).resolve().parent.parent
DATA_DIR = BASE / "data"


def main() -> None:
    genes = pd.read_csv(DATA_DIR / "01_genes_all_strains.csv")
    edges = pd.read_csv(DATA_DIR / "02_ortholog_edges_all_strains.csv")
    genes["key"] = genes["organism_name"] + "::" + genes["locus_tag"]
    edges["key"] = edges["organism_name"] + "::" + edges["locus_tag"]

    cyanorak = edges[edges["source"] == "cyanorak"][["key", "ortholog_group_id"]]
    cyanorak = cyanorak.rename(columns={"ortholog_group_id": "cyanorak_group"})
    assert cyanorak["key"].is_unique, "expected exactly 0 or 1 cyanorak group per gene"

    eggnog1 = edges[(edges["source"] == "eggnog") & (edges["specificity_rank"] == 1)][
        ["key", "ortholog_group_id"]
    ]
    eggnog1 = eggnog1.rename(columns={"ortholog_group_id": "eggnog1_group"})
    eggnog1 = eggnog1.drop_duplicates("key")  # a gene could in principle map to >1 rank-1 group; keep first

    merged = genes.merge(cyanorak, on="key", how="left").merge(eggnog1, on="key", how="left")

    def assign(row):
        if pd.notna(row["cyanorak_group"]):
            return row["cyanorak_group"], "cyanorak"
        if pd.notna(row["eggnog1_group"]):
            return row["eggnog1_group"], "eggnog_rank1"
        return None, "true_orphan"

    assigned = merged.apply(assign, axis=1, result_type="expand")
    merged["final_group_id"], merged["group_source"] = assigned[0], assigned[1]

    out_cols = [
        "organism_name", "locus_tag", "gene_name", "gene_category", "product",
        "group", "final_group_id", "group_source",
    ]
    merged[out_cols].to_csv(DATA_DIR / "03_gene_group_assignment.csv", index=False)
    print(f"Wrote {len(merged)} gene rows to data/03_gene_group_assignment.csv")
    print(merged["group_source"].value_counts().to_string())

    # Group-landscape summary (distinct groups + avg strain coverage), for the record.
    landscape_rows = []
    for (source, rank), grp in edges.groupby(["source", "specificity_rank"]):
        landscape_rows.append({
            "source": source,
            "specificity_rank": rank,
            "distinct_groups": grp["ortholog_group_id"].nunique(),
            "edges": len(grp),
            "avg_strains_per_group": grp.groupby("ortholog_group_id")["organism_name"].nunique().mean(),
        })
    landscape_df = pd.DataFrame(landscape_rows)
    landscape_df.to_csv(DATA_DIR / "04_group_landscape_summary.csv", index=False)
    print(f"\nWrote {len(landscape_df)} rows to data/04_group_landscape_summary.csv")

    # True-orphan count per strain.
    orphans = merged[merged["group_source"] == "true_orphan"]
    orphan_summary = orphans.groupby(["organism_name", "group"]).size().reset_index(name="true_orphan_count")
    orphan_unknown = orphans[orphans["gene_category"] == "Unknown"].groupby("organism_name").size()
    orphan_summary["unknown_function_count"] = orphan_summary["organism_name"].map(orphan_unknown).fillna(0).astype(int)
    orphan_summary.to_csv(DATA_DIR / "05_orphan_summary_by_strain.csv", index=False)
    print(f"\nWrote {len(orphan_summary)} rows to data/05_orphan_summary_by_strain.csv")
    print(orphan_summary.sort_values("true_orphan_count", ascending=False).to_string(index=False))


if __name__ == "__main__":
    main()
