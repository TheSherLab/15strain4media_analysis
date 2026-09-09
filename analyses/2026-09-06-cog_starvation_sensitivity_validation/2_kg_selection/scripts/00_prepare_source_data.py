"""
Step 2 (reopened 2026-09-09) -- prepare the source inputs for this analysis.

APPROACH CHANGE (2026-09-09, researcher-directed): the earlier approach used
one researcher-supplied Cyanorak ID per COG. That picks a single
"representative" gene per COG and misses the rest of the COG's membership
(e.g. COG0443 landed on a DUF3181 hypothetical instead of dnaK1/2/3). The
researcher's comparative-genomics significance test was run on COG *copy
number* -- i.e. every gene carrying that COG number in a strain -- so the
in-vivo expression validation must use the same gene set. We now take the
full per-strain COG membership straight from the researcher's own
whole-genome annotation (`WhloeGenome_AllStrains_Concated.xlsx`).

Inputs (external, not in this repo):
  .../Revisions/paper 29042025/Dataset 3.xlsx  (sheet "Significant NorMixed COGs ")
      -- the 41 COGs, Direction tag, protein names, functional ordering
  .../Analysis- home/Whole genome/WhloeGenome_AllStrains_Concated.xlsx
      -- every CDS in all 13 annotated study-strain genomes, with its
         cog_numbers and cluster_number (Cyanorak group)
  .../Analysis- home/Whole genome/cogs_with_sig_p_values.xlsx
      -- the 41 significant COGs with p-values + median copy numbers
         (used here only as a Direction cross-check)

Outputs:
  data/00_source_normixed_cogs.csv       -- 41 COGs, Direction, protein, ordering
  data/00_wholegenome_cog_members.csv    -- every gene carrying one of the 41
                                            COGs, per strain (the new gene universe)
  data/00_cog_to_ck_id_mapping.csv       -- SUPERSEDED single-CK-per-COG map,
                                            kept for the record / step-1 trace

Usage: uv run analyses/2026-09-06-cog_starvation_sensitivity_validation/2_kg_selection/scripts/00_prepare_source_data.py
"""

from pathlib import Path

import pandas as pd

BASE = Path(__file__).resolve().parent.parent
DATA_DIR = BASE / "data"

DESKTOP = Path(r"C:\Users\yaras\OneDrive\Desktop\Analysis for paper")
SOURCE_XLSX = DESKTOP / "Revisions" / "paper 29042025" / "Dataset 3.xlsx"
SHEET_NAME = "Significant NorMixed COGs "
WHOLEGENOME_XLSX = DESKTOP / "Analysis- home" / "Whole genome" / "WhloeGenome_AllStrains_Concated.xlsx"
SIG_PVAL_XLSX = DESKTOP / "Analysis- home" / "Whole genome" / "cogs_with_sig_p_values.xlsx"

# SUPERSEDED (2026-09-09): kept only so the step-1 -> step-2 trail stays in the
# repo. The analysis no longer resolves COGs through these single IDs.
COG_TO_CK_ID = {
    "COG2159": "CK_00002647", "COG1596": "CK_00048203", "COG3206": "CK_00001962",
    "COG0477": "CK_00002121", "COG2192": "CK_00002335", "COG2931": "CK_00042429",
    "COG0188": "CK_00008107", "COG1454": "CK_00049429", "COG0412": "CK_00001881",
    "COG1403": "CK_00000339", "COG0673": "CK_00000253", "COG2370": "CK_00000499",
    "COG3727": "CK_00011865", "COG0640": "CK_00001532", "COG0697": "CK_00000195",
    "COG0367": "CK_00057301", "COG3210": "CK_00001608", "COG0443": "CK_00000457",
    "COG5651": "CK_00006528", "COG1212": "CK_00048376", "COG1462": "CK_00002347",
    "COG4252": "CK_00002384", "COG3754": "CK_00007567", "COG2114": "CK_00002384",
    "COG0454": "CK_00001405", "COG1187": "CK_00000762", "COG2214": "CK_00000214",
    "COG1472": "CK_00000205", "COG0810": "CK_00001270", "COG0082": "CK_00000482",
    "COG3839": "CK_00001436", "COG0330": "CK_00041773", "COG0845": "CK_00000467",
    "COG2274": "CK_00056818", "COG1989": "CK_00001720", "COG0488": "CK_00000335",
    "COG0484": "CK_00000350", "COG1793": "CK_00001728", "COG1787": "CK_00003490",
    "COG2010": "CK_00001007", "COG0399": "CK_00009007",
}


def main() -> None:
    # --- 1. the 41-COG list: Direction, protein name, functional ordering ---
    df = pd.read_excel(SOURCE_XLSX, sheet_name=SHEET_NAME)
    df = df.dropna(subset=["COGs"]).reset_index(drop=True)
    df.columns = [c.strip() for c in df.columns]
    # COG0640 (smtB/ArsR) has a blank Direction cell; the researcher confirmed
    # (2026-09-09, on the ordered screenshot) it is N. Two tags: N (13), mixed (28).
    df["Direction"] = df["Direction"].fillna("N")
    assert set(df["Direction"]) == {"N", "mixed"}, sorted(set(df["Direction"]))
    for col in ("Protien", "Annotation", "Wider annotation", "General annotation"):
        if col in df.columns:
            df[col] = df[col].astype("string").str.replace("\xa0", " ", regex=False).str.strip()
    df["order"] = range(len(df))  # sheet order == the researcher's functional ordering

    out1 = DATA_DIR / "00_source_normixed_cogs.csv"
    df.to_csv(out1, index=False)
    print(f"Wrote {len(df)} rows to {out1}")
    print(df["Direction"].value_counts().to_string())
    cogs41 = set(df["COGs"])
    assert len(cogs41) == 41

    # Direction cross-check against cogs_with_sig_p_values.xlsx (uses N / "Other").
    sig = pd.read_excel(SIG_PVAL_XLSX)
    sig_dir = dict(zip(sig["cog_numbers"], sig["Direction"]))
    mism = []
    for _, r in df.iterrows():
        other = sig_dir.get(r["COGs"])
        want_n = r["Direction"] == "N"
        has_n = other == "N"
        if other is not None and want_n != has_n:
            mism.append((r["COGs"], r["Direction"], other))
    print(f"Direction cross-check vs cogs_with_sig_p_values: {len(mism)} mismatch(es): {mism}")
    # COG0640 is the one expected mismatch (blank in sig file -> 'Other'; researcher says N).

    # --- 2. full per-strain COG membership from the whole-genome annotation ---
    wg = pd.read_excel(WHOLEGENOME_XLSX)
    wg = wg[wg["cog_numbers"].isin(cogs41)].copy()
    wg["locus_name"] = wg["attributes"].str.extract(r"Name=([^;]+)")
    wg["ckpro_id"] = wg["attributes"].str.extract(r"ID=(CK_Pro_[^;]+)")
    wg["cluster_number"] = wg["attributes"].str.extract(r"cluster_number=(CK_\d+)")
    wg["product"] = wg["attributes"].str.extract(r"product=([^;]+)")
    wg["product"] = wg["product"].str.replace("%2C", ",", regex=False)
    dir_map = dict(zip(df["COGs"], df["Direction"]))
    wg["direction"] = wg["cog_numbers"].map(dir_map)

    members = wg[[
        "strain", "cog_numbers", "direction", "locus_name", "ckpro_id",
        "cluster_number", "product", "seqid", "start", "end", "strand",
    ]].sort_values(["cog_numbers", "strain", "start"]).reset_index(drop=True)
    out2 = DATA_DIR / "00_wholegenome_cog_members.csv"
    members.to_csv(out2, index=False)
    print(f"\nWrote {len(members)} gene rows to {out2}")
    print(f"  strains: {sorted(members['strain'].unique())}")
    print(f"  COGs covered: {members['cog_numbers'].nunique()} / 41"
          f"  (missing: {sorted(cogs41 - set(members['cog_numbers']))})")
    print(f"  distinct Cyanorak clusters: {members['cluster_number'].nunique()}")
    ev = ["MED4", "MIT9313", "NATL2A", "MIT9312"]
    print("\n  genes per COG in the 4 evidence strains (top 12 by size):")
    piv = (members[members["strain"].isin(ev)]
           .groupby(["cog_numbers", "strain"]).size().unstack(fill_value=0).reindex(columns=ev, fill_value=0))
    piv["max"] = piv.max(axis=1)
    print(piv.sort_values("max", ascending=False).head(12).to_string())

    # --- 3. superseded single-CK map (kept for the record) ---
    mapping_df = pd.DataFrame([{"cog_numbers": k, "CK_ID": v, "status": "SUPERSEDED_2026-09-09"}
                               for k, v in COG_TO_CK_ID.items()])
    out3 = DATA_DIR / "00_cog_to_ck_id_mapping.csv"
    mapping_df.to_csv(out3, index=False)
    print(f"\nWrote {len(mapping_df)} rows to {out3} (superseded; kept for the step-1 trace)")
    assert not (cogs41 - set(mapping_df["cog_numbers"]))


if __name__ == "__main__":
    main()
