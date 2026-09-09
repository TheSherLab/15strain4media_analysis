"""
Step 2 (reopened 2026-09-09) -- triangulate every COG's Cyanorak-ID
resolution against two independent routes, mirroring the prior walkthrough
analysis's `2_kg_selection/scripts/06_name_vs_cyanorak_crosscheck.py`.

The researcher supplied ONE Cyanorak group ID per COG, applied to all 15
strains. Two failure modes this checks, for the 4 evidence-bearing strains
(MED4 / MIT9313 / NATL2A / MIT9312):

  1. name route:  resolve_gene(<protein-name token>) -> locus tag(s) in the
     strain, independent of the supplied Cyanorak ID. Compare to the
     CK-ID route locus.
  2. reverse homolog route:  for each locus the CK-ID route produced,
     gene_homologs(locus, source="cyanorak") -> which Cyanorak group(s)
     that locus actually belongs to. Confirm the supplied CK_ID is one.

Verdicts per (COG, strain): agree / agree_name_multi / DISAGREE /
name_only / name_only_multi / cyanorak_only / both_absent.
`reverse_group_ok` flags whether the CK-route locus actually carries the
supplied Cyanorak ID.

Inputs:  data/00_source_normixed_cogs.csv, data/02_gene_locus_resolution.csv
Outputs: data/06_name_vs_cyanorak_crosscheck.csv

Usage: uv run analyses/2026-09-06-cog_starvation_sensitivity_validation/2_kg_selection/scripts/05_name_vs_cyanorak_crosscheck.py
"""

import re
from pathlib import Path

import pandas as pd
from multiomics_explorer import resolve_gene, gene_homologs, GraphConnection

DATA_DIR = Path(__file__).resolve().parent.parent / "data"
EVIDENCE_STRAINS = [
    "Prochlorococcus MED4",
    "Prochlorococcus MIT9313",
    "Prochlorococcus NATL2A",
    "Prochlorococcus MIT9312",
]

STOPWORDS = {
    "and", "or", "the", "including", "family", "protein", "domain", "like",
    "type", "activity", "repair", "very", "short", "path", "member", "site",
    "specific", "a", "of", "with", "poorly",
}


def name_tokens(protein: str, annotation: str) -> list[str]:
    ids: set[str] = set()
    text = protein if isinstance(protein, str) and protein.strip() else ""
    for chunk in re.split(r"[/,()]", text):
        for tok in chunk.split():
            tok = tok.strip(" .-")
            if not tok or tok.lower() in STOPWORDS:
                continue
            if re.fullmatch(r"[A-Za-z][A-Za-z0-9]{1,7}", tok):
                ids.add(tok)
    if not ids and isinstance(annotation, str) and annotation.strip():
        ids.add(annotation.strip())
    return sorted(ids)


def main() -> None:
    source = pd.read_csv(DATA_DIR / "00_source_normixed_cogs.csv")
    source.columns = [c.strip() for c in source.columns]
    res = pd.read_csv(DATA_DIR / "02_gene_locus_resolution.csv")

    ck_route: dict[tuple, list] = {}
    for _, r in res.iterrows():
        if r["status"] == "resolved":
            ck_route.setdefault((r["cog_number"], r["organism_name"]), []).append(r["locus_tag"])
    cog_ck = {r["cog_number"]: r["CK_ID"] for _, r in res.iterrows()}

    # --- token -> {organism: set(locus)} : one resolve_gene per unique token ---
    per_cog_tokens = {
        g["COGs"]: name_tokens(g["Protien"], g["Annotation"]) for _, g in source.iterrows()
    }
    all_tokens = sorted({t for toks in per_cog_tokens.values() for t in toks})
    token_hits: dict[str, dict[str, set]] = {}
    # --- locus -> set(cyanorak groups) : batched gene_homologs ---
    all_ck_loci = sorted({l for loci in ck_route.values() for l in loci})
    locus_groups: dict[str, set] = {}

    with GraphConnection() as conn:
        for tok in all_tokens:
            hits: dict[str, set] = {s: set() for s in EVIDENCE_STRAINS}
            try:
                r = resolve_gene(identifier=tok, limit=80, conn=conn)
                for hit in r.get("results", []):
                    org = hit.get("organism_name")
                    if org in hits:
                        hits[org].add(hit["locus_tag"])
            except Exception as e:  # noqa: BLE001
                print(f"  resolve_gene({tok!r}) raised {e}")
            token_hits[tok] = hits

        for i in range(0, len(all_ck_loci), 40):
            batch = all_ck_loci[i:i + 40]
            gh = gene_homologs(locus_tags=batch, source="cyanorak", limit=200, conn=conn)
            for row in gh.get("results", []):
                locus_groups.setdefault(row["locus_tag"], set()).add(
                    row["group_id"].replace("cyanorak:", "")
                )

    rows = []
    for _, g in source.iterrows():
        cog = g["COGs"]
        toks = per_cog_tokens[cog]
        ck_id = cog_ck[cog]
        for strain in EVIDENCE_STRAINS:
            ck_loci = sorted(set(ck_route.get((cog, strain), [])))
            names = sorted({l for t in toks for l in token_hits.get(t, {}).get(strain, set())})

            if ck_loci and names:
                if set(ck_loci) == set(names):
                    verdict = "agree"
                elif set(ck_loci) & set(names):
                    verdict = "agree_name_multi"
                else:
                    verdict = "DISAGREE"
            elif ck_loci and not names:
                verdict = "cyanorak_only"
            elif names and not ck_loci:
                verdict = "name_only" if len(names) == 1 else "name_only_multi"
            else:
                verdict = "both_absent"

            reverse_ok = None
            reverse_groups = ""
            if ck_loci:
                gset = set().union(*(locus_groups.get(l, set()) for l in ck_loci))
                reverse_ok = ck_id in gset
                reverse_groups = "; ".join(sorted(gset))

            rows.append({
                "cog_number": cog, "protein_name": g["Protien"], "direction": g["Direction"],
                "cyanorak_id": ck_id, "organism_name": strain,
                "name_tokens": "; ".join(toks),
                "ck_route_loci": "; ".join(ck_loci),
                "name_route_loci": "; ".join(names),
                "verdict": verdict,
                "reverse_group_ok": reverse_ok,
                "reverse_groups_for_ck_loci": reverse_groups,
            })

    df = pd.DataFrame(rows)
    out = DATA_DIR / "06_name_vs_cyanorak_crosscheck.csv"
    df.to_csv(out, index=False)
    print(f"Wrote {len(df)} rows to {out}\n")
    print(df["verdict"].value_counts().to_string())
    for flag, mask in [
        ("DISAGREE", df["verdict"] == "DISAGREE"),
        ("name_only / name_only_multi", df["verdict"].str.startswith("name_only")),
        ("reverse_group_ok == False", df["reverse_group_ok"] == False),  # noqa: E712
    ]:
        sub = df[mask]
        print(f"\n--- {flag} ({len(sub)}) ---")
        if len(sub):
            print(sub[["cog_number", "protein_name", "organism_name", "cyanorak_id",
                       "ck_route_loci", "name_route_loci", "reverse_groups_for_ck_loci"]].to_string(index=False))


if __name__ == "__main__":
    main()
