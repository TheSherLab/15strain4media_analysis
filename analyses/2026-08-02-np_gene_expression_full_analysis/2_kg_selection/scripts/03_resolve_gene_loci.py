"""Resolve the 54 N/P-acquisition genes (from 01_export_evidence_gene_list.py)
to locus tags, restricted to the 6 Prochlorococcus strains that have an
axenic nitrogen or phosphorus starvation experiment (from
02_list_np_experiments.py) -- MED4, MIT9301, MIT9312, MIT9313, NATL2A,
and SS120 (CCMP1375).

Rule 2 (gene identity): locus tags, not gene names -- this re-derives the
mapping fresh from the KG rather than reusing the prior triage analysis's
cached results (which recorded only strain names, not locus tags, for
"has evidence").

Resolution method mirrors the prior triage analysis (steps 02-04 there):
- 53 of 54 genes resolve directly by name via resolve_gene.
- phoA is the one ontology-routed exception: its Gene.gene_name field is
  empty in the KG, so it was originally found via KEGG orthology term
  K01077, which is shared with phoB -- candidates are filtered to (a)
  exclude loci already resolved by name as phoB, and (b) keep only the
  "probable alkaline phosphatase" product (the K01077 Pfam family also
  contains unrelated DedA-family genes, which this filter excludes).
  No genes in this 54-gene list use the narX1->narM name remap from the
  prior analysis (narX1/narM are not in the evidence list).
- PMM707/PMM719/PMM721/PMM722 are not gene symbols -- they are unpadded
  MED4 locus tags used as gene labels in the source spreadsheet (known
  issue from the prior analysis: real format is PMM0707 etc.). These are
  looked up directly via gene_overview on the zero-padded locus tag,
  restricted to MED4 (the only strain they identify), not via
  resolve_gene name matching.

Inputs:  data/01_evidence_gene_list.csv, data/02_np_experiments.csv
Outputs: data/03_gene_loci.csv (one row per gene x strain x locus_tag)

Usage:
  uv run python analyses/2026-08-02-np_gene_expression_full_analysis/2_kg_selection/scripts/03_resolve_gene_loci.py
"""
from pathlib import Path

import pandas as pd
from multiomics_explorer import GraphConnection, gene_overview, genes_by_ontology, resolve_gene

DATA_DIR = Path(__file__).resolve().parents[1] / "data"

PHOA_KEGG_KO = "kegg.orthology:K01077"
PHOA_PRODUCT_SUBSTR = "phosphatase"

# gene_name_raw -> zero-padded MED4 locus tag (see module docstring).
UNPADDED_MED4_LOCUS_TAGS = {
    "PMM707": "PMM0707",
    "PMM719": "PMM0719",
    "PMM721": "PMM0721",
    "PMM722": "PMM0722",
}
MED4 = "Prochlorococcus MED4"


def resolve_by_name(gene_name: str, strains: list[str], conn: GraphConnection) -> list[dict]:
    alt_names = [n for n in gene_name.split("/") if n.strip()]
    seen = {}
    for name in alt_names:
        result = resolve_gene(identifier=name, organism="Prochlorococcus", limit=200, conn=conn)
        for r in result["results"]:
            if r["organism_name"] in strains:
                seen[(r["organism_name"], r["locus_tag"])] = r
    return list(seen.values())


def resolve_phoa(strains: list[str], phob_loci: set[str], conn: GraphConnection) -> list[dict]:
    rows = []
    for strain in strains:
        hits = genes_by_ontology(
            ontology="kegg",
            term_ids=[PHOA_KEGG_KO],
            organism=strain,
            min_gene_set_size=0,
            max_gene_set_size=100000,
            limit=50,
            conn=conn,
        )
        candidate_loci = [r["locus_tag"] for r in hits["results"] if r["locus_tag"] not in phob_loci]
        if not candidate_loci:
            continue
        ov = gene_overview(locus_tags=candidate_loci, limit=len(candidate_loci), conn=conn)
        for row in ov["results"]:
            if PHOA_PRODUCT_SUBSTR in (row.get("product") or "").lower():
                rows.append(
                    {
                        "organism_name": strain,
                        "locus_tag": row["locus_tag"],
                        "gene_name_in_kg": row.get("gene_name") or "",
                        "product": row.get("product") or "",
                    }
                )
    return rows


def main() -> None:
    genes = pd.read_csv(DATA_DIR / "01_evidence_gene_list.csv")
    experiments = pd.read_csv(DATA_DIR / "02_np_experiments.csv")
    strains = sorted(experiments["organism_name"].unique())
    print(f"{len(genes)} genes to resolve across {len(strains)} strains: {strains}")

    records = []
    with GraphConnection() as conn:
        phob = resolve_by_name("phoB", strains, conn)
        phob_loci = {r["locus_tag"] for r in phob}

        for _, row in genes.iterrows():
            gene_name = row["gene_name_raw"]
            if gene_name == "phoA":
                hits = resolve_phoa(strains, phob_loci, conn)
            elif gene_name in UNPADDED_MED4_LOCUS_TAGS:
                if MED4 not in strains:
                    hits = []
                else:
                    padded = UNPADDED_MED4_LOCUS_TAGS[gene_name]
                    ov = gene_overview(locus_tags=[padded], limit=1, conn=conn)
                    hits = [
                        {
                            "organism_name": MED4,
                            "locus_tag": r["locus_tag"],
                            "gene_name_in_kg": r.get("gene_name") or "",
                            "product": r.get("product") or "",
                        }
                        for r in ov["results"]
                    ]
            else:
                hits = [
                    {
                        "organism_name": r["organism_name"],
                        "locus_tag": r["locus_tag"],
                        "gene_name_in_kg": r.get("gene_name") or "",
                        "product": r.get("product") or "",
                    }
                    for r in resolve_by_name(gene_name, strains, conn)
                ]
            for h in hits:
                records.append(
                    {
                        "n_or_p": row["n_or_p"],
                        "gene_name_raw": gene_name,
                        "organism_name": h["organism_name"],
                        "locus_tag": h["locus_tag"],
                        "gene_name_in_kg": h["gene_name_in_kg"],
                        "product": h["product"],
                    }
                )

    out = pd.DataFrame(records).drop_duplicates(subset=["gene_name_raw", "organism_name", "locus_tag"])
    out_path = DATA_DIR / "03_gene_loci.csv"
    out.to_csv(out_path, index=False)
    print(f"\nWrote {len(out)} (gene x strain x locus) rows to {out_path}")

    print("\nGenes resolved in >=1 relevant strain:", out["gene_name_raw"].nunique(), "of", len(genes))
    missing = set(genes["gene_name_raw"]) - set(out["gene_name_raw"])
    if missing:
        print("Genes with NO locus in any of the 6 relevant strains:", sorted(missing))

    print("\nLoci per (n_or_p, strain):")
    print(out.groupby(["n_or_p", "organism_name"])["locus_tag"].nunique().to_string())


if __name__ == "__main__":
    main()
