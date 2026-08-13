"""
Step 2 -- bulk pull: every gene for the 15 study strains, plus every
Gene_in_ortholog_group edge for those genes (both Cyanorak and eggnog
sources). This is the raw material the genome-wide association scan
(step 5) will run on.

No higher-level MCP/package tool covers "every gene x its ortholog
group(s), scoped to a set of organisms" -- genes_by_homolog_group takes
group_ids as input (wrong direction), gene_homologs takes locus_tags
(one gene at a time). run_cypher is the documented escape hatch for
exactly this shape of question (docs://guide/start_here).

Inputs: KG (via multiomics_explorer.run_cypher)
Outputs: data/01_genes_all_strains.csv
         data/02_ortholog_edges_all_strains.csv

Usage: uv run analyses/2026-08-12-np_sensitivity_gene_scan/2_kg_selection/scripts/01_pull_gene_ortholog_data.py
"""

from pathlib import Path

import pandas as pd
from multiomics_explorer import run_cypher

BASE = Path(__file__).resolve().parent.parent
DATA_DIR = BASE / "data"

STRAINS = [
    "Prochlorococcus MED4", "Prochlorococcus MIT9515", "Prochlorococcus MIT9202",
    "Prochlorococcus MIT9215", "Prochlorococcus MIT0604", "Prochlorococcus AS9601",
    "Prochlorococcus MIT9312", "Prochlorococcus MIT1314", "Prochlorococcus NATL1A",
    "Prochlorococcus MIT9301", "Prochlorococcus SB", "Prochlorococcus NATL2A",
    "Prochlorococcus PAC1", "Prochlorococcus MIT9313", "Prochlorococcus MIT1327",
]

N_SENSITIVE = {
    "Prochlorococcus MED4", "Prochlorococcus MIT9515", "Prochlorococcus MIT9202",
    "Prochlorococcus MIT9215", "Prochlorococcus MIT0604", "Prochlorococcus AS9601",
    "Prochlorococcus MIT9312", "Prochlorococcus MIT1314", "Prochlorococcus NATL1A",
}
MIXED = {
    "Prochlorococcus MIT9301", "Prochlorococcus SB", "Prochlorococcus NATL2A",
    "Prochlorococcus PAC1", "Prochlorococcus MIT9313", "Prochlorococcus MIT1327",
}


def strain_list_literal(strains: list[str]) -> str:
    quoted = ", ".join(f'"{s}"' for s in strains)
    return f"[{quoted}]"


def main() -> None:
    assert N_SENSITIVE | MIXED == set(STRAINS)
    assert not (N_SENSITIVE & MIXED)

    strains_literal = strain_list_literal(STRAINS)

    genes_query = f"""
    MATCH (g:Gene)
    WHERE g.organism_name IN {strains_literal}
    RETURN g.organism_name AS organism_name, g.locus_tag AS locus_tag,
           g.gene_name AS gene_name, g.gene_category AS gene_category,
           g.product AS product
    """
    genes_result = run_cypher(query=genes_query, limit=100_000)
    genes_df = pd.DataFrame(genes_result["results"])
    genes_df["group"] = genes_df["organism_name"].map(
        lambda o: "N_sensitive" if o in N_SENSITIVE else "mixed"
    )
    genes_df.to_csv(DATA_DIR / "01_genes_all_strains.csv", index=False)
    print(f"Genes pulled: {len(genes_df)} (truncated={genes_result.get('truncated')})")
    print(genes_df.groupby("organism_name").size().to_string())

    edges_query = f"""
    MATCH (g:Gene)-[:Gene_in_ortholog_group]->(og:OrthologGroup)
    WHERE g.organism_name IN {strains_literal}
    RETURN g.organism_name AS organism_name, g.locus_tag AS locus_tag,
           og.id AS ortholog_group_id, og.source AS source,
           og.specificity_rank AS specificity_rank,
           og.taxonomic_level AS taxonomic_level,
           og.consensus_gene_name AS consensus_gene_name
    """
    edges_result = run_cypher(query=edges_query, limit=300_000)
    edges_df = pd.DataFrame(edges_result["results"])
    edges_df.to_csv(DATA_DIR / "02_ortholog_edges_all_strains.csv", index=False)
    print(f"\nOrtholog-group edges pulled: {len(edges_df)} (truncated={edges_result.get('truncated')})")
    print(edges_df.groupby(["source", "specificity_rank"]).size().to_string())


if __name__ == "__main__":
    main()
