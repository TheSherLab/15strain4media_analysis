"""Run hit_rate() (step 4) across all group x experiment-set combinations
needed for H1 (matched), H2 (cross), H3 (noise check), and the
positive-control sanity check.

Groups:
- H1_N              : nitrogen-annotated target genes  x nitrogen experiments
- H1_P              : phosphorus-annotated target genes x phosphorus experiments
- H2_Ngenes_Pexp    : nitrogen-annotated target genes  x phosphorus experiments (cross)
- H2_Pgenes_Nexp    : phosphorus-annotated target genes x nitrogen experiments (cross)
- positive_control_N: N-side positive-control genes (ntcA, glnA, amtB/amt1, ureA) x nitrogen experiments
- positive_control_P: P-side positive-control genes (pstS, phoA, phoB, phoR) x phosphorus experiments
- background_N      : background gene set x nitrogen experiments (H3 -- valid noise check)

background_P (H3 for phosphorus) is deliberately NOT computed: all 5
phosphorus experiments are pre-filtered by the source publications to
already-significant genes (table_scope significant_only/filtered_subset
-- see 2_kg_selection/notebook.md's experiment table), so a "background"
sample drawn from those same tables is not a random sample of the
genome -- it inherits the same selection bias as the target genes, which
would make any H3-style comparison uninterpretable rather than merely
weak. Researcher decision (2026-08-02): skip it rather than report and
caveat it.

Inputs:  2_kg_selection/data/02_np_experiments.csv, 03_gene_loci.csv
         3_analysis_framing/data/01_background_genes.csv
Outputs: data/01_hit_rate_per_gene.csv    (one row per gene x group)
         data/01_hit_rate_summary.csv     (one row per group, aggregate)

Usage:
  uv run python analyses/2026-08-02-np_gene_expression_full_analysis/5_analyze/scripts/01_compute_hit_rates.py
"""
import sys
from pathlib import Path

import pandas as pd
from multiomics_explorer import GraphConnection

STEP_DIR = Path(__file__).resolve().parents[1]
ANALYSIS_ROOT = STEP_DIR.parent
sys.path.insert(0, str(ANALYSIS_ROOT / "4_methods"))
from hit_rate import hit_rate  # noqa: E402

KG_SELECTION_DATA = ANALYSIS_ROOT / "2_kg_selection" / "data"
FRAMING_DATA = ANALYSIS_ROOT / "3_analysis_framing" / "data"
DATA_DIR = STEP_DIR / "data"

POSITIVE_CONTROL_N = ["ntcA", "glnA", "amtB/amt1", "ureA"]
POSITIVE_CONTROL_P = ["pstS", "phoA", "phoB", "phoR"]


def pairs_from(df: pd.DataFrame) -> list[tuple[str, str]]:
    return list(df[["organism_name", "locus_tag"]].itertuples(index=False, name=None))


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
        "H2_Ngenes_Pexp": (pairs_from(n_target), p_exp_ids),
        "H2_Pgenes_Nexp": (pairs_from(p_target), n_exp_ids),
        "positive_control_N": (pairs_from(pos_n), n_exp_ids),
        "positive_control_P": (pairs_from(pos_p), p_exp_ids),
        "background_N": (pairs_from(background), n_exp_ids),
    }

    per_gene_rows = []
    summary_rows = []
    with GraphConnection() as conn:
        for group, (pairs, exp_ids) in groups.items():
            result = hit_rate(pairs, exp_ids, conn=conn)
            pg = result["per_gene"].copy()
            pg.insert(0, "group", group)
            per_gene_rows.append(pg)

            agg = result["aggregate"]
            summary_rows.append({"group": group, "n_genes_in_set": len(set(pairs)), **agg})
            print(f"{group}: {agg}")

    per_gene = pd.concat(per_gene_rows, ignore_index=True)
    summary = pd.DataFrame(summary_rows)

    per_gene.to_csv(DATA_DIR / "01_hit_rate_per_gene.csv", index=False)
    summary.to_csv(DATA_DIR / "01_hit_rate_summary.csv", index=False)
    print(f"\nWrote {len(per_gene)} per-gene rows and {len(summary)} group summaries to {DATA_DIR}")


if __name__ == "__main__":
    main()
