"""Freeze the 54-gene evidence list from the prior triage analysis
(analyses/2026-07-13-n_p_genes_in_vivo_experiments/) as this analysis's
local input.

That analysis resolved 91 N/P acquisition genes against the KG and found
54 with in vivo expression evidence in some Prochlorococcus strain
(has_expression_evidence == True). This script copies just those 54 rows
into this analysis's own data/ so this analysis does not depend on
another analysis directory's files at run time.

Input:  ../../2026-07-13-n_p_genes_in_vivo_experiments/data/n_p_genes_experimental_evidence_final.csv
Output: data/01_evidence_gene_list.csv (54 rows: n_or_p, gene_name_raw,
        cyanorak_id, resolution_status, n_strains_with_evidence)

Usage:
  uv run python 01_export_evidence_gene_list.py
"""
from pathlib import Path

import pandas as pd

DATA_DIR = Path(__file__).resolve().parents[1] / "data"
SOURCE = (
    Path(__file__).resolve().parents[3]
    / "2026-07-13-n_p_genes_in_vivo_experiments"
    / "data"
    / "n_p_genes_experimental_evidence_final.csv"
)


def main() -> None:
    df = pd.read_csv(SOURCE)
    evidence = df[df["has_expression_evidence"] == True]  # noqa: E712
    keep_cols = [
        "n_or_p",
        "gene_name_raw",
        "cyanorak_id",
        "resolution_status",
        "n_strains_with_evidence",
    ]
    out = evidence[keep_cols].reset_index(drop=True)

    out_path = DATA_DIR / "01_evidence_gene_list.csv"
    out.to_csv(out_path, index=False)
    print(f"Wrote {len(out)} genes to {out_path}")
    print(out["n_or_p"].value_counts().to_string())
    print()
    print("resolution_status breakdown:")
    print(out["resolution_status"].value_counts().to_string())


if __name__ == "__main__":
    main()
