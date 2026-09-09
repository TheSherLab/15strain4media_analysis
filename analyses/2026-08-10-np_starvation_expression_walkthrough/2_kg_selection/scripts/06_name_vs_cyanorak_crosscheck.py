"""
Step 2 (redo, 2026-09-08) -- triangulate every gene's Cyanorak-ID
resolution against an independent name-based lookup.

02_resolve_genes.py keys ONLY on the researcher's Cyanorak ortholog-group
ID. The phoE bug (multi-member group) showed one failure mode; this script
checks a different one -- whether the Cyanorak ID in the spreadsheet points
at the gene its name/ORF says it should. For each (gene x in-scope strain):

  - resolve_gene(<gene_name>) and resolve_gene(<normalised ORF>) -> the
    locus tag(s) the KG gives for that identifier, independent of Cyanorak.
  - compare to the locus tag the Cyanorak-ID route produced
    (02_gene_locus_resolution.csv).

Verdicts per (gene, strain):
  agree            -- name route and Cyanorak route give the same locus
  disagree         -- both resolved, to DIFFERENT loci  (investigate!)
  name_ambiguous   -- name route gave >1 locus in the strain (can't confirm)
  name_only        -- name resolved, Cyanorak route did not
  cyanorak_only    -- Cyanorak route resolved, name route did not
  both_absent      -- neither resolved (expected: gene not in that strain)

Inputs: ../data/00_source_gene_list.csv, ../data/02_gene_locus_resolution.csv
Outputs: ../data/06_name_vs_cyanorak_crosscheck.csv

Usage: uv run analyses/2026-08-10-np_starvation_expression_walkthrough/2_kg_selection/scripts/06_name_vs_cyanorak_crosscheck.py
"""

import re
from pathlib import Path

import pandas as pd
from multiomics_explorer import resolve_gene, GraphConnection

DATA_DIR = Path(__file__).resolve().parent.parent / "data"
STRAINS = [
    "Prochlorococcus MED4",
    "Prochlorococcus MIT9312",
    "Prochlorococcus MIT9313",
    "Prochlorococcus NATL2A",
]


def candidate_identifiers(gene_name: str) -> list[str]:
    """The researcher's names include real gene names (phoE, pstS), slash
    aliases (amtB/amt1), and MED4 ORF stubs (PMM707 -> PMM0707)."""
    ids = set()
    for part in re.split(r"[/,]", gene_name):
        part = part.strip()
        if not part:
            continue
        ids.add(part)
        m = re.fullmatch(r"([A-Za-z_]+?)(\d+)", part)
        if m:  # PMM707 -> PMM0707 ; PMT9312_25 -> PMT9312_0025 (best-effort zero-pad to 4)
            stem, num = m.group(1), m.group(2)
            ids.add(f"{stem}{int(num):04d}")
    return sorted(ids)


def main() -> None:
    genes = pd.read_csv(DATA_DIR / "00_source_gene_list.csv")
    genes.columns = [c.strip() for c in genes.columns]
    res = pd.read_csv(DATA_DIR / "02_gene_locus_resolution.csv")
    ck_locus = {
        (r["gene_name"], r["organism_name"]): r["locus_tag"]
        for _, r in res.iterrows()
        if r["status"] == "resolved"
    }

    rows = []
    with GraphConnection() as conn:
        for _, g in genes.iterrows():
            gene_name = g["Gene_name"] if pd.notna(g["Gene_name"]) else "(unnamed)"
            n_or_p = g["N or P gene?"]
            ck_id = g["Cyanorak Ids"] if pd.notna(g["Cyanorak Ids"]) else None
            if gene_name == "(unnamed)":
                continue

            # name-route hits: {organism_name: set(locus_tags)}
            name_hits: dict[str, set] = {s: set() for s in STRAINS}
            for ident in candidate_identifiers(gene_name):
                try:
                    r = resolve_gene(identifier=ident, conn=conn)
                except Exception as e:  # noqa: BLE001
                    print(f"  resolve_gene({ident!r}) raised {e}")
                    continue
                for hit in r.get("results", []):
                    org = hit.get("organism_name")
                    if org in name_hits:
                        name_hits[org].add(hit["locus_tag"])

            for strain in STRAINS:
                ck = ck_locus.get((gene_name, strain))
                names = sorted(name_hits[strain])
                if ck and names:
                    if ck in names and len(names) == 1:
                        verdict = "agree"
                    elif ck in names:
                        verdict = "agree_name_also_has_others"
                    else:
                        verdict = "DISAGREE"
                elif ck and not names:
                    verdict = "cyanorak_only"
                elif names and not ck:
                    verdict = "name_only" if len(names) == 1 else "name_only_ambiguous"
                else:
                    verdict = "both_absent"
                rows.append({
                    "n_or_p": n_or_p, "gene_name": gene_name, "cyanorak_id": ck_id,
                    "organism_name": strain, "cyanorak_route_locus": ck,
                    "name_route_loci": "; ".join(names), "verdict": verdict,
                })

    df = pd.DataFrame(rows)
    out = DATA_DIR / "06_name_vs_cyanorak_crosscheck.csv"
    df.to_csv(out, index=False)
    print(f"Wrote {len(df)} rows to {out}\n")
    print(df["verdict"].value_counts().to_string())
    print("\n--- DISAGREE rows (name route and Cyanorak route give different loci) ---")
    dis = df[df["verdict"] == "DISAGREE"]
    print(dis.to_string(index=False) if len(dis) else "  (none)")
    print("\n--- name_only / name_only_ambiguous (Cyanorak route missed a locus the name found) ---")
    no = df[df["verdict"].str.startswith("name_only")]
    print(no.to_string(index=False) if len(no) else "  (none)")


if __name__ == "__main__":
    main()
