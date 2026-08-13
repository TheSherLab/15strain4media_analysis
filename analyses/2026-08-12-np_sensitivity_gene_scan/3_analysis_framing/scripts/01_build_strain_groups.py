"""
Step 3 -- materialize the two strain groupings every gene will be tested
against: the researcher-supplied sensitivity phenotype (N-sensitive vs.
mixed), and the KG-verified ecotype clade collapsed to HL/LL (co-defined
2026-08-13, after the clade-confound discussion).

Inputs: KG (via multiomics_explorer.list_organisms, for `clade` --
        verified here, not assumed from memory)
Outputs: data/01_strain_group_membership.csv

Usage: uv run analyses/2026-08-12-np_sensitivity_gene_scan/3_analysis_framing/scripts/01_build_strain_groups.py
"""

from pathlib import Path

import pandas as pd
from multiomics_explorer import list_organisms

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
HL_CLADES = {"HLI", "HLII"}
LL_CLADES = {"LLI", "LLII", "LLIII", "LLIV"}


def main() -> None:
    result = list_organisms(organism_names=STRAINS, limit=15)
    rows = []
    for r in result["results"]:
        clade = r["clade"]
        ecotype = "HL" if clade in HL_CLADES else ("LL" if clade in LL_CLADES else None)
        assert ecotype is not None, f"unrecognized clade {clade!r} for {r['organism_name']}"
        rows.append({
            "organism_name": r["organism_name"],
            "clade": clade,
            "ecotype_group": ecotype,
            "sensitivity_group": "N_sensitive" if r["organism_name"] in N_SENSITIVE else "mixed",
            "gene_count": r["gene_count"],
        })
    df = pd.DataFrame(rows).sort_values("organism_name")
    assert len(df) == 15
    df.to_csv(DATA_DIR / "01_strain_group_membership.csv", index=False)
    print(f"Wrote {len(df)} rows to data/01_strain_group_membership.csv")
    print(df.to_string(index=False))
    print("\nCross-tab (sensitivity x ecotype):")
    print(pd.crosstab(df["sensitivity_group"], df["ecotype_group"]).to_string())


if __name__ == "__main__":
    main()
