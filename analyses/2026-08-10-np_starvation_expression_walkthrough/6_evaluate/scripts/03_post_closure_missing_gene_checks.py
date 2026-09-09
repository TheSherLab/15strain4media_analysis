"""
Post-closure QC (2026-09-06) -- the researcher asked, after this analysis
was already closed, why specific target genes show "no data" for certain
experiments. Re-verified two distinct mechanisms directly against the live
KG rather than trusting this analysis's own prior summary text:

1. Read et al. 2017 MED4 RNA-seq (nitrogen): its table_scope_detail says
   "Top 50% of genes by expression level" -- confirmed this filter applies
   genome-wide (a gene missing at one timepoint is missing at all of them),
   using cynD/pipX/urtD as examples (all resolved and evidenced elsewhere
   in MED4).
2. phoE / locus PMT_2631 in MIT9313: confirmed this locus has zero DE rows
   in ANY MIT9313 experiment in the entire KG -- a platform/probe coverage
   gap, not a declared significance filter like (1).

See gaps_and_friction.md (2026-09-06 entry) for the write-up.

Usage: uv run analyses/2026-08-10-np_starvation_expression_walkthrough/6_evaluate/scripts/03_post_closure_missing_gene_checks.py
"""

from pathlib import Path

import pandas as pd
from multiomics_explorer import differential_expression_by_gene, GraphConnection, resolve_gene

BASE = Path(__file__).resolve().parent.parent
DATA_DIR = BASE / "data"

READ_EXP = "10.1038/ismej.2017.88_nitrogen_stress_ndepleted_pro99_medium_med4_rnaseq"
TARGET_LOCI = {"PMM0372": "cynD", "PMM0393": "pipX", "PMM0973": "urtD"}

TOLONEN_MIT9313 = "10.1038/msb4100087_nitrogen_nitrogen_deprivation_mit9313_mit9313_microarray"
PHOE_MIT9313_LOCUS = "PMT_2631"


def main() -> None:
    lines = []

    with GraphConnection() as conn:
        lines.append("=" * 70)
        lines.append("1. Read et al. 2017 MED4 RNA-seq -- ALL timepoints for cynD/pipX/urtD")
        lines.append("=" * 70)
        result = differential_expression_by_gene(organism="MED4", experiment_ids=[READ_EXP], limit=None, conn=conn)
        de = pd.DataFrame(result["results"])
        lines.append(f"Total rows: {len(de)}, distinct genes: {de['locus_tag'].nunique()}, "
                      f"timepoints: {sorted(de['timepoint'].dropna().unique().tolist())}")
        rows_out = []
        for locus, name in TARGET_LOCI.items():
            rows = de[de["locus_tag"] == locus]
            present = len(rows) > 0
            lines.append(f"  {name} ({locus}): {'present at ' + str(sorted(rows['timepoint'].tolist())) if present else 'NOT PRESENT at any timepoint'}")
            rows_out.append({"gene_name": name, "locus_tag": locus, "n_rows_in_read_experiment": len(rows)})

        lines.append("")
        lines.append("=" * 70)
        lines.append("2. phoE (PMT_2631, MIT9313) -- Tolonen all timepoints + all MIT9313 experiments")
        lines.append("=" * 70)
        try:
            g = resolve_gene(PHOE_MIT9313_LOCUS, organism="MIT9313", conn=conn)
            lines.append(f"Gene record: {g['results']}")
        except Exception as e:
            lines.append(f"resolve_gene failed: {e}")

        result2 = differential_expression_by_gene(organism="MIT9313", experiment_ids=[TOLONEN_MIT9313], limit=None, conn=conn)
        de2 = pd.DataFrame(result2["results"])
        lines.append(f"Tolonen MIT9313: {len(de2)} rows, {de2['locus_tag'].nunique()} distinct genes, "
                      f"timepoints: {sorted(de2['timepoint'].dropna().unique().tolist())}")
        rows2 = de2[de2["locus_tag"] == PHOE_MIT9313_LOCUS]
        lines.append(f"Rows for {PHOE_MIT9313_LOCUS} in Tolonen MIT9313: {len(rows2)}")

        result3 = differential_expression_by_gene(organism="MIT9313", locus_tags=[PHOE_MIT9313_LOCUS], limit=None, conn=conn)
        de3 = pd.DataFrame(result3["results"])
        lines.append(f"Rows for {PHOE_MIT9313_LOCUS} across ALL MIT9313 experiments in the KG: {len(de3)}")

    out_txt = DATA_DIR / "03_post_closure_missing_gene_checks.txt"
    out_txt.write_text("\n".join(lines), encoding="utf-8")
    pd.DataFrame(rows_out).to_csv(DATA_DIR / "03_read_target_gene_presence.csv", index=False)
    print("\n".join(lines))
    print(f"\nWrote {out_txt} and data/03_read_target_gene_presence.csv")


if __name__ == "__main__":
    main()
