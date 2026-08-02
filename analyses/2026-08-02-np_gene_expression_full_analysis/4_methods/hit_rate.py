"""Shared method for this analysis: the significant-hit rate.

For a set of (strain, locus_tag) genes and a set of experiment IDs, what
fraction of (gene x experiment x timepoint) tests come back significant,
split by direction (up/down)?

This single function is applied six ways across steps 5-6 to test the
three hypotheses locked in step 3:
  H1 (matched, directional):  N genes x N experiments, P genes x P experiments
  H2 (cross, non-directional): N genes x P experiments, P genes x N experiments
  H3 (noise check):            background genes x N experiments, x P experiments
It is also used once as a sanity check on the 8 positive-control genes.

Definitions:
- A "test" is one (locus_tag, experiment_id, timepoint) row returned by
  differential_expression_by_gene(significant_only=False). Genes with no
  row for a given experiment were not measured there and are excluded
  from that experiment's denominator (not counted as "not significant").
- experiment_ids may span multiple organisms (e.g. "all nitrogen
  experiments" includes both MED4 and MIT9313 experiment IDs);
  differential_expression_by_gene requires a single-organism experiment
  set per call, so this function first narrows the passed experiment_ids
  to each organism's own subset via list_experiments before querying DE.
- expression_status is one of {significant_up, significant_down,
  not_significant} in every experiment this analysis uses (verified in
  step 4's worked example below -- no ambiguous "not_known" cells here).
- rate_significant = (n_up + n_down) / n_tests
- rate_up, rate_down = n_up / n_tests, n_down / n_tests

Usage (see 4_methods/notebook.md for the worked example):
    from hit_rate import hit_rate
    result = hit_rate([("Prochlorococcus MED4", "PMM1082")], experiment_ids=[...])
"""
from collections import defaultdict

import pandas as pd
from multiomics_explorer import GraphConnection, differential_expression_by_gene, list_experiments


def hit_rate(
    organism_locus_pairs: list[tuple[str, str]],
    experiment_ids: list[str],
    conn: GraphConnection | None = None,
) -> dict:
    by_organism: dict[str, list[str]] = defaultdict(list)
    for organism, locus_tag in organism_locus_pairs:
        by_organism[organism].append(locus_tag)

    rows = []
    for organism, loci in by_organism.items():
        loci = sorted(set(loci))

        # experiment_ids may span multiple organisms -- narrow to this
        # organism's own subset before calling differential_expression_by_gene,
        # which requires a single-organism experiment set.
        exp_lookup = list_experiments(
            experiment_ids=experiment_ids, organism=organism, verbose=False, limit=None, conn=conn
        )
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
                    "organism_name": organism,
                    "locus_tag": r["locus_tag"],
                    "experiment_id": r["experiment_id"],
                    "timepoint": r.get("timepoint"),
                    "expression_status": r["expression_status"],
                }
            )

    de = pd.DataFrame(rows, columns=["organism_name", "locus_tag", "experiment_id", "timepoint", "expression_status"])

    if de.empty:
        per_gene = pd.DataFrame(
            columns=["organism_name", "locus_tag", "n_tests", "n_up", "n_down", "n_not_significant", "rate_significant"]
        )
    else:
        counts = (
            de.groupby(["organism_name", "locus_tag", "expression_status"])
            .size()
            .unstack(fill_value=0)
        )
        for col in ["significant_up", "significant_down", "not_significant"]:
            if col not in counts.columns:
                counts[col] = 0
        counts["n_tests"] = counts["significant_up"] + counts["significant_down"] + counts["not_significant"]
        counts["rate_significant"] = (counts["significant_up"] + counts["significant_down"]) / counts["n_tests"]
        per_gene = counts.reset_index().rename(
            columns={"significant_up": "n_up", "significant_down": "n_down", "not_significant": "n_not_significant"}
        )[["organism_name", "locus_tag", "n_tests", "n_up", "n_down", "n_not_significant", "rate_significant"]]

    n_tests = int(per_gene["n_tests"].sum())
    n_up = int(per_gene["n_up"].sum())
    n_down = int(per_gene["n_down"].sum())
    n_not_significant = int(per_gene["n_not_significant"].sum())

    aggregate = {
        "n_genes_tested": per_gene["locus_tag"].nunique() if not per_gene.empty else 0,
        "n_tests": n_tests,
        "n_up": n_up,
        "n_down": n_down,
        "n_not_significant": n_not_significant,
        "rate_significant": (n_up + n_down) / n_tests if n_tests else float("nan"),
        "rate_up": n_up / n_tests if n_tests else float("nan"),
        "rate_down": n_down / n_tests if n_tests else float("nan"),
    }

    return {"per_gene": per_gene, "aggregate": aggregate, "raw": de}
