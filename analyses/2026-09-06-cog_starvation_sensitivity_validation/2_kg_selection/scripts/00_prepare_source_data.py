"""
Step 2 -- copy and clean the two source inputs for this analysis:
1. The researcher's "Significant NorMixed COGs" sheet (41 COGs, Direction
   tag, p-values, protein name, annotation).
2. The researcher's own COG-to-Cyanorak-ID mapping (verified in step 1
   against the live KG: 41/41 COGs covered, 364 gene rows resolved across
   the 15 study strains, zero unmatched groups/strains).

Inputs (external, not in this repo):
  C:\\Users\\yaras\\OneDrive\\Desktop\\Analysis for paper\\Revisions\\paper 29042025\\Dataset 3.xlsx
    (sheet "Significant NorMixed COGs ")
  researcher-supplied COG->CK_ID table (transcribed from the researcher's
    screenshot in this analysis's step-1 dialogue)

Outputs: data/00_source_normixed_cogs.csv, data/00_cog_to_ck_id_mapping.csv

Usage: uv run analyses/2026-09-06-cog_starvation_sensitivity_validation/2_kg_selection/scripts/00_prepare_source_data.py
"""

from pathlib import Path

import pandas as pd

BASE = Path(__file__).resolve().parent.parent
DATA_DIR = BASE / "data"

SOURCE_XLSX = (
    r"C:\Users\yaras\OneDrive\Desktop\Analysis for paper\Revisions\paper 29042025\Dataset 3.xlsx"
)
SHEET_NAME = "Significant NorMixed COGs "

# Transcribed from the researcher's own COG->Cyanorak-ID table (step 1,
# verified complete and clean against the live KG before adoption).
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
    df = pd.read_excel(SOURCE_XLSX, sheet_name=SHEET_NAME)
    df = df.dropna(subset=["COGs"]).reset_index(drop=True)
    df.columns = [c.strip() for c in df.columns]
    df["Direction"] = df["Direction"].fillna("unlabeled")

    out1 = DATA_DIR / "00_source_normixed_cogs.csv"
    df.to_csv(out1, index=False)
    print(f"Wrote {len(df)} rows to {out1}")
    print(df["Direction"].value_counts())

    mapping_df = pd.DataFrame(
        [{"cog_numbers": k, "CK_ID": v} for k, v in COG_TO_CK_ID.items()]
    )
    out2 = DATA_DIR / "00_cog_to_ck_id_mapping.csv"
    mapping_df.to_csv(out2, index=False)
    print(f"Wrote {len(mapping_df)} rows to {out2}")

    # sanity: every COG in the source sheet has a mapping entry, and vice versa
    missing_mapping = set(df["COGs"]) - set(mapping_df["cog_numbers"])
    extra_mapping = set(mapping_df["cog_numbers"]) - set(df["COGs"])
    assert not missing_mapping, f"COGs missing a CK_ID mapping: {missing_mapping}"
    assert not extra_mapping, f"Mapping has COGs not in the source sheet: {extra_mapping}"
    print("Sanity check passed: source sheet and CK_ID mapping cover the same 41 COGs.")


if __name__ == "__main__":
    main()
