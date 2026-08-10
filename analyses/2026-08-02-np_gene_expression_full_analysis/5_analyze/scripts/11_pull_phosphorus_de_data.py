"""Mirror of 06_pull_nitrogen_volcano_data.py, for the phosphorus side:
pull raw DE data for all 53 usable genes (both N- and P-annotated),
restricted to the 5 phosphorus starvation experiments.

This is H1_P (P genes x P experiments) + H2_Ngenes_Pexp (N genes x P
experiments) combined into one gene-level table with log2fc/padj kept.
N genes here will be sparse: H2_Ngenes_Pexp was only 6 rows / 2 genes
in step 5 (phosphorus tables are narrow and pre-filtered to
already-significant genes -- see 2_kg_selection/notebook.md), so most
N-annotated genes have zero rows in this file. That sparsity is a
finding in itself, not a bug -- kept as-is, not padded or worked around.

Input:  2_kg_selection/data/02_np_experiments.csv, 03_gene_loci.csv
Output: data/11_phosphorus_de_data.csv

Usage:
  uv run python analyses/2026-08-02-np_gene_expression_full_analysis/5_analyze/scripts/11_pull_phosphorus_de_data.py
"""
from collections import defaultdict
from pathlib import Path

import pandas as pd
from multiomics_explorer import GraphConnection, differential_expression_by_gene, list_experiments

STEP_DIR = Path(__file__).resolve().parents[1]
ANALYSIS_ROOT = STEP_DIR.parent
KG_SELECTION_DATA = ANALYSIS_ROOT / "2_kg_selection" / "data"
DATA_DIR = STEP_DIR / "data"


def main() -> None:
    experiments = pd.read_csv(KG_SELECTION_DATA / "02_np_experiments.csv")
    gene_loci = pd.read_csv(KG_SELECTION_DATA / "03_gene_loci.csv")
    p_exp_ids = experiments.loc[experiments["n_or_p"] == "P", "experiment_id"].tolist()

    gene_meta = gene_loci.set_index(["organism_name", "locus_tag"])[["n_or_p", "gene_name_raw"]]

    by_organism: dict[str, list[str]] = defaultdict(list)
    for organism, locus_tag in gene_loci[["organism_name", "locus_tag"]].itertuples(index=False, name=None):
        by_organism[organism].append(locus_tag)

    rows = []
    with GraphConnection() as conn:
        for organism, loci in by_organism.items():
            loci = sorted(set(loci))
            exp_lookup = list_experiments(experiment_ids=p_exp_ids, organism=organism, limit=None, conn=conn)
            org_experiment_ids = [r["experiment_id"] for r in exp_lookup["results"]]
            if not org_experiment_ids:
                continue
            result = differential_expression_by_gene(
                organism=organism,
                locus_tags=loci,
                experiment_ids=org_experiment_ids,
                significant_only=False,
                verbose=True,
                limit=None,
                conn=conn,
            )
            for r in result["results"]:
                meta = gene_meta.loc[(organism, r["locus_tag"])]
                rows.append(
                    {
                        "n_or_p": meta["n_or_p"],
                        "gene_name_raw": meta["gene_name_raw"],
                        "organism_name": organism,
                        "locus_tag": r["locus_tag"],
                        "experiment_id": r["experiment_id"],
                        "omics_type": r["omics_type"],
                        "timepoint": r.get("timepoint"),
                        "log2fc": r.get("log2fc"),
                        "padj": r.get("padj"),
                        "expression_status": r["expression_status"],
                    }
                )

    out = pd.DataFrame(rows)
    out_path = DATA_DIR / "11_phosphorus_de_data.csv"
    out.to_csv(out_path, index=False)
    print(f"Wrote {len(out)} rows to {out_path}")
    print(out.groupby(["n_or_p", "omics_type"]).size().to_string())
    print(f"\nDistinct genes: {out.groupby('n_or_p')['gene_name_raw'].nunique().to_dict()}")


if __name__ == "__main__":
    main()
