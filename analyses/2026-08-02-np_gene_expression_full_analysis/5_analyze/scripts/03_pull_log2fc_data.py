"""Pull raw log2FC values (not just the up/down/not-significant call
used in 01_compute_hit_rates.py) for the 5 groups that have an
interpretable comparison: H1_N, H1_P, positive_control_N,
positive_control_P, background_N. (background_P is excluded throughout
this analysis -- see 5_analyze/notebook.md.)

hit_rate() (step 4) only keeps expression_status; this script re-queries
the same (gene, experiment) pairs but keeps log2fc and omics_type too,
so effect size (not just hit/no-hit) can be visualized, and so
cross-platform pooling can be avoided per this methodology's statistical
rigor rule (RNA-seq/proteomics/microarray log2FC magnitudes are not
directly comparable -- compare within platform only).

Inputs:  2_kg_selection/data/02_np_experiments.csv, 03_gene_loci.csv
         3_analysis_framing/data/01_background_genes.csv
Outputs: data/03_log2fc_raw.csv (one row per gene x experiment x timepoint,
         all 5 groups, with log2fc, padj, omics_type, expression_status)

Usage:
  uv run python analyses/2026-08-02-np_gene_expression_full_analysis/5_analyze/scripts/03_pull_log2fc_data.py
"""
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd
from multiomics_explorer import GraphConnection, differential_expression_by_gene, list_experiments

STEP_DIR = Path(__file__).resolve().parents[1]
ANALYSIS_ROOT = STEP_DIR.parent
KG_SELECTION_DATA = ANALYSIS_ROOT / "2_kg_selection" / "data"
FRAMING_DATA = ANALYSIS_ROOT / "3_analysis_framing" / "data"
DATA_DIR = STEP_DIR / "data"

POSITIVE_CONTROL_N = ["ntcA", "glnA", "amtB/amt1", "ureA"]
POSITIVE_CONTROL_P = ["pstS", "phoA", "phoB", "phoR"]


def pairs_from(df: pd.DataFrame) -> list[tuple[str, str]]:
    return list(df[["organism_name", "locus_tag"]].itertuples(index=False, name=None))


def pull(group: str, pairs: list[tuple[str, str]], experiment_ids: list[str], conn: GraphConnection) -> list[dict]:
    by_organism: dict[str, list[str]] = defaultdict(list)
    for organism, locus_tag in pairs:
        by_organism[organism].append(locus_tag)

    rows = []
    for organism, loci in by_organism.items():
        loci = sorted(set(loci))
        exp_lookup = list_experiments(experiment_ids=experiment_ids, organism=organism, limit=None, conn=conn)
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
            rows.append(
                {
                    "group": group,
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
    return rows


def main() -> None:
    experiments = pd.read_csv(KG_SELECTION_DATA / "02_np_experiments.csv")
    gene_loci = pd.read_csv(KG_SELECTION_DATA / "03_gene_loci.csv")
    background = pd.read_csv(FRAMING_DATA / "01_background_genes.csv")

    n_exp_ids = experiments.loc[experiments["n_or_p"] == "N", "experiment_id"].tolist()
    p_exp_ids = experiments.loc[experiments["n_or_p"] == "P", "experiment_id"].tolist()

    n_target = gene_loci[gene_loci["n_or_p"] == "N"]
    p_target = gene_loci[gene_loci["n_or_p"] == "P"]
    pos_n = gene_loci[gene_loci["gene_name_raw"].isin(POSITIVE_CONTROL_N)]
    pos_p = gene_loci[gene_loci["gene_name_raw"].isin(POSITIVE_CONTROL_P)]

    groups = {
        "H1_N": (pairs_from(n_target), n_exp_ids),
        "H1_P": (pairs_from(p_target), p_exp_ids),
        "positive_control_N": (pairs_from(pos_n), n_exp_ids),
        "positive_control_P": (pairs_from(pos_p), p_exp_ids),
        "background_N": (pairs_from(background), n_exp_ids),
    }

    all_rows = []
    with GraphConnection() as conn:
        for group, (pairs, exp_ids) in groups.items():
            rows = pull(group, pairs, exp_ids, conn)
            print(f"{group}: {len(rows)} rows")
            all_rows.extend(rows)

    out = pd.DataFrame(all_rows)
    out_path = DATA_DIR / "03_log2fc_raw.csv"
    out.to_csv(out_path, index=False)
    print(f"\nWrote {len(out)} rows to {out_path}")
    print("\nBy group x omics_type (row counts):")
    print(out.groupby(["group", "omics_type"]).size().to_string())


if __name__ == "__main__":
    main()
