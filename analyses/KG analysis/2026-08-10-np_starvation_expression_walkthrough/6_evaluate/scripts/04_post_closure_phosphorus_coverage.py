"""
Post-closure QC (2026-09-06) -- quantify what fraction of each strain's
genome each of the 5 in-scope phosphorus experiments' KG table actually
covers, to check whether the narrow coverage already flagged for Martiny
(step 3/6) holds across all 5 phosphorus experiments regardless of method
(RNA-seq, proteomics, microarray), or is specific to Martiny/microarray.

See gaps_and_friction.md (2026-09-06 entry) for the write-up.

Usage: uv run analyses/2026-08-10-np_starvation_expression_walkthrough/6_evaluate/scripts/04_post_closure_phosphorus_coverage.py
"""

from pathlib import Path

import pandas as pd
from multiomics_explorer import differential_expression_by_gene, GraphConnection, list_organisms

BASE = Path(__file__).resolve().parent.parent
DATA_DIR = BASE / "data"

EXPERIMENTS = [
    ("Lin RNA-seq (NATL2A)", "10.1111/1462-2920.13104_phosphorus_plimited_natl2a_rnaseq_uninfected", "NATL2A"),
    ("Fuszard proteomics (MIT9312)", "10.1186/2046-9063-8-7_pi_limitation_mit9312_itraq", "MIT9312"),
    ("Fuszard proteomics (NATL2A)", "10.1186/2046-9063-8-7_pi_limitation_natl2a_itraq", "NATL2A"),
    ("Martiny microarray (MED4)", "10.1073/pnas.0601301103_phosphorus_phosphate_starvation_med4_microarray", "MED4"),
    ("Martiny microarray (MIT9313)", "10.1073/pnas.0601301103_phosphorus_phosphate_starvation_mit9313_microarray", "MIT9313"),
]


def main() -> None:
    rows = []
    with GraphConnection() as conn:
        orgs = list_organisms(conn=conn)
        orgs_df = pd.DataFrame(orgs["results"])
        genome_size = dict(zip(
            orgs_df["organism_name"].str.replace("Prochlorococcus ", "", regex=False),
            orgs_df["gene_count"],
        ))

        for label, exp_id, organism in EXPERIMENTS:
            r = differential_expression_by_gene(organism=organism, experiment_ids=[exp_id], limit=None, conn=conn)
            d = pd.DataFrame(r["results"])
            n_genes = d["locus_tag"].nunique() if len(d) else 0
            n_rows = len(d)
            gsize = genome_size.get(organism)
            pct = round(100 * n_genes / gsize, 1) if gsize else None
            rows.append({
                "experiment": label, "experiment_id": exp_id, "strain": organism,
                "n_rows": n_rows, "distinct_genes": n_genes, "genome_size": gsize,
                "pct_genome_covered": pct,
            })

    df = pd.DataFrame(rows)
    out_path = DATA_DIR / "04_phosphorus_genome_coverage.csv"
    df.to_csv(out_path, index=False)
    print(df.to_string(index=False))
    print(f"\nWrote {len(df)} rows to {out_path}")


if __name__ == "__main__":
    main()
