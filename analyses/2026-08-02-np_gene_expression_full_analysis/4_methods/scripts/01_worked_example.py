"""Worked example / toy verification for hit_rate() (step 4 QC), before
applying it to the full 53-gene target set in step 5.

Traces ntcA and glnA (MED4's canonical nitrogen-regulation and
nitrogen-assimilation marker genes, from step 3's positive-control list)
through hit_rate() against MED4's 4 nitrogen experiments, and
cross-checks the function's aggregate counts against a hand-computed
tally from the raw differential_expression_by_gene rows -- confirming
the counting logic (per-row status tally, not e.g. double-counting
timepoints or dropping rows) before trusting it on 53+303 genes.

Usage:
  uv run python analyses/2026-08-02-np_gene_expression_full_analysis/4_methods/scripts/01_worked_example.py
"""
import sys
from pathlib import Path

import pandas as pd
from multiomics_explorer import GraphConnection, differential_expression_by_gene, resolve_gene

STEP_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(STEP_DIR))
from hit_rate import hit_rate  # noqa: E402

KG_SELECTION_DATA = STEP_DIR.parent / "2_kg_selection" / "data"
MED4 = "Prochlorococcus MED4"


def main() -> None:
    experiments = pd.read_csv(KG_SELECTION_DATA / "02_np_experiments.csv")
    n_exp_ids = experiments.loc[
        (experiments["n_or_p"] == "N") & (experiments["organism_name"] == MED4), "experiment_id"
    ].tolist()
    print(f"MED4 nitrogen experiments ({len(n_exp_ids)}): {n_exp_ids}")

    with GraphConnection() as conn:
        ntca = resolve_gene(identifier="ntcA", organism=MED4, limit=10, conn=conn)["results"]
        glna = resolve_gene(identifier="glnA", organism=MED4, limit=10, conn=conn)["results"]
        print(f"\nntcA -> {[(r['organism_name'], r['locus_tag']) for r in ntca]}")
        print(f"glnA -> {[(r['organism_name'], r['locus_tag']) for r in glna]}")

        pairs = [(MED4, r["locus_tag"]) for r in ntca] + [(MED4, r["locus_tag"]) for r in glna]
        result = hit_rate(pairs, n_exp_ids, conn=conn)

        print("\nhit_rate() per-gene output:")
        print(result["per_gene"].to_string(index=False))
        print("\nhit_rate() aggregate:", result["aggregate"])

        # Hand-check: re-pull raw DE rows directly and tally by hand.
        print("\n--- Hand-check (direct query, manual tally) ---")
        for organism, locus_tag in pairs:
            de = differential_expression_by_gene(
                organism=organism,
                locus_tags=[locus_tag],
                experiment_ids=n_exp_ids,
                significant_only=False,
                verbose=True,
                limit=None,
                conn=conn,
            )
            rows = de["results"]
            n_up = sum(1 for r in rows if r["expression_status"] == "significant_up")
            n_down = sum(1 for r in rows if r["expression_status"] == "significant_down")
            n_ns = sum(1 for r in rows if r["expression_status"] == "not_significant")
            print(
                f"{locus_tag}: {len(rows)} rows -> up={n_up} down={n_down} "
                f"not_sig={n_ns} (sum check: {n_up + n_down + n_ns == len(rows)})"
            )
            for r in rows:
                print(
                    f"    {r['experiment_id'][:50]:<50} tp={r['timepoint']:<8} "
                    f"status={r['expression_status']:<17} log2fc={r.get('log2fc')}"
                )


if __name__ == "__main__":
    main()
