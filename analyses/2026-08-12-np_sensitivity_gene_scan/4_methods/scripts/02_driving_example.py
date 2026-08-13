"""
Step 4 -- driving example: run the methods module end-to-end on one real
ortholog group, both statistical tests, both groupings (sensitivity,
ecotype), using the real gene/strain data step 2 and step 3 already pulled
and QC'd. Proves the module works against real data before step 5 wires it
up into the full ~5,000+-group loop.

Group selection rule (mechanical, not outcome-based -- chosen before looking
at any p-value): the first ortholog group, in sorted final_group_id order,
that has (a) real presence/absence variation across the 15 strains (carried
by more than 1 but fewer than all 15) and (b) real copy-number variation
(at least one strain with more than 1 copy) -- so the demo exercises both
test functions rather than landing on a trivial universal-single-copy gene
(which is a real, valid case, just not one that shows both code paths).

FDR correction and candidate classification are NOT run here: both require
the full batch of p-values from a given (test type x grouping) combination,
which only exists once step 5 runs the loop across all ortholog groups. This
example shows raw p-values and direction only.

Inputs:
  ../../2_kg_selection/data/03_gene_group_assignment.csv
  ../../3_analysis_framing/data/01_strain_group_membership.csv
Output:
  data/02_driving_example_readout.csv

Usage:
  uv run python analyses/2026-08-12-np_sensitivity_gene_scan/4_methods/scripts/02_driving_example.py
"""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from np_sensitivity_scan import copy_number_test, presence_absence_test  # noqa: E402

STEP_DIR = Path(__file__).resolve().parents[1]
KG_DATA = STEP_DIR.parent / "2_kg_selection" / "data"
FRAMING_DATA = STEP_DIR.parent / "3_analysis_framing" / "data"
OUT_DIR = STEP_DIR / "data"


def main() -> None:
    genes = pd.read_csv(KG_DATA / "03_gene_group_assignment.csv")
    strains = pd.read_csv(FRAMING_DATA / "01_strain_group_membership.csv")

    strain_gene_counts = dict(zip(strains["organism_name"], strains["gene_count"]))
    all_strains = list(strains["organism_name"])
    sensitive = list(strains.loc[strains["sensitivity_group"] == "N_sensitive", "organism_name"])
    mixed = list(strains.loc[strains["sensitivity_group"] == "mixed", "organism_name"])
    hl = list(strains.loc[strains["ecotype_group"] == "HL", "organism_name"])
    ll = list(strains.loc[strains["ecotype_group"] == "LL", "organism_name"])
    assert len(sensitive) == 9 and len(mixed) == 6, "sensitivity split must match the locked 9-vs-6 grouping"
    assert len(hl) == 10 and len(ll) == 5, "ecotype split must match step 3's 10-vs-5 grouping"

    tested = genes[genes["group_source"] != "true_orphan"]
    copy_matrix = (
        tested.groupby(["final_group_id", "organism_name"]).size().unstack(fill_value=0)
    )
    copy_matrix = copy_matrix.reindex(columns=all_strains, fill_value=0)

    def carries_count(row):
        return (row > 0).sum()

    def has_multi_copy(row):
        return (row > 1).any()

    candidates = copy_matrix[
        (copy_matrix.apply(carries_count, axis=1) > 1)
        & (copy_matrix.apply(carries_count, axis=1) < 15)
        & (copy_matrix.apply(has_multi_copy, axis=1))
    ]
    chosen_group_id = sorted(candidates.index)[0]
    copy_counts = copy_matrix.loc[chosen_group_id].to_dict()

    example_rows = tested[tested["final_group_id"] == chosen_group_id]
    gene_names = sorted(set(example_rows["gene_name"].dropna())) or ["(no gene_name recorded)"]
    product = example_rows["product"].mode().iat[0] if not example_rows["product"].mode().empty else ""

    print(f"Chosen ortholog group: {chosen_group_id}")
    print(f"  gene name(s) across strains: {gene_names}")
    print(f"  product: {product}")
    print(f"  copy counts by strain: {copy_counts}\n")

    sens_pa = presence_absence_test(copy_counts, sensitive, mixed, "N_sensitive", "mixed")
    eco_pa = presence_absence_test(copy_counts, hl, ll, "HL", "LL")
    sens_cn = copy_number_test(copy_counts, strain_gene_counts, sensitive, mixed, "N_sensitive", "mixed")
    eco_cn = copy_number_test(copy_counts, strain_gene_counts, hl, ll, "HL", "LL")

    print("Presence/absence vs sensitivity:", sens_pa)
    print("Presence/absence vs ecotype:    ", eco_pa)
    print("Copy number vs sensitivity:     ", sens_cn)
    print("Copy number vs ecotype:         ", eco_cn)

    rows = []
    for grouping, result in [("sensitivity", sens_pa), ("ecotype", eco_pa)]:
        rows.append({"final_group_id": chosen_group_id, "test": "presence_absence", "grouping": grouping, **result})
    for grouping, result in [("sensitivity", sens_cn), ("ecotype", eco_cn)]:
        if result is not None:
            rows.append({"final_group_id": chosen_group_id, "test": "copy_number", "grouping": grouping, **result})
        else:
            rows.append({"final_group_id": chosen_group_id, "test": "copy_number", "grouping": grouping, "p_value": None, "direction": "skipped (not variable)"})

    out = pd.DataFrame(rows)
    OUT_DIR.mkdir(exist_ok=True)
    out_path = OUT_DIR / "02_driving_example_readout.csv"
    out.to_csv(out_path, index=False)
    print(f"\nWrote {out_path}")


if __name__ == "__main__":
    main()
