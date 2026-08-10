"""
Step 5 -- strip plot: % upregulated per gene, N-target vs. P-target vs.
N-background, to visually back up the H3 bootstrap result (p<0.0001,
02_h3_summary.txt) with the actual per-gene distribution rather than just
a single p-value.

Background genes are grouped by locus_tag (not gene_name): unlike target
genes, which share the same gene_name across strains and so get one
pooled percentage per name (04_pct_upregulated.csv), background genes are
a different random sample per strain -- there is no cross-strain identity
to group by. A MED4 background gene can appear in up to 3 of the 4
nitrogen background experiments (proteomics/RNA-seq/microarray all being
MED4); a MIT9313 background gene appears in only 1 (Tolonen MIT9313 is
the only MIT9313 experiment in the background pool). This is coarser
than the target genes' cross-strain grouping -- noted in notebook.md, not
hidden.

Inputs: data/04_pct_upregulated.csv,
        ../../3_analysis_framing/data/02_negative_background_pool.csv
Outputs: data/05_pct_scatter_background.csv (background per-locus rates)
         figures/03_pct_scatter.png

Usage: uv run analyses/2026-08-10-np_starvation_expression_walkthrough/5_analyze/scripts/05_pct_scatter.py
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

BASE = Path(__file__).resolve().parent.parent
DATA_DIR = BASE / "data"
FIG_DIR = BASE / "figures"
FRAMING_DATA = BASE.parent / "3_analysis_framing" / "data"

COLOR_N = "#2a78d6"
COLOR_P = "#e34948"
COLOR_BG = "#898781"
INK = "#0b0b0b"
SECONDARY_INK = "#52514e"
MUTED = "#898781"
SURFACE = "#fcfcfb"
GRID = "#e1e0d9"


def main() -> None:
    target = pd.read_csv(DATA_DIR / "04_pct_upregulated.csv")
    target = target[target["n_tested"] > 0]

    bg = pd.read_csv(FRAMING_DATA / "02_negative_background_pool.csv")
    bg_summary = bg.groupby("locus_tag").agg(
        n_tested=("expression_status", "size"),
        n_up=("expression_status", lambda s: (s == "significant_up").sum()),
    ).reset_index()
    bg_summary["pct_up"] = 100 * bg_summary["n_up"] / bg_summary["n_tested"]
    bg_summary.to_csv(DATA_DIR / "05_pct_scatter_background.csv", index=False)
    print(f"Background: {len(bg_summary)} distinct genes (by locus_tag), "
          f"mean pct_up={bg_summary['pct_up'].mean():.1f}%, median={bg_summary['pct_up'].median():.1f}%")

    groups = [
        ("N target genes", target[target["n_or_p"] == "N"]["pct_up"].values, COLOR_N),
        ("N background pool", bg_summary["pct_up"].values, COLOR_BG),
        ("P target genes", target[target["n_or_p"] == "P"]["pct_up"].values, COLOR_P),
    ]
    for label, vals, _ in groups:
        print(f"{label}: n={len(vals)}, mean={np.mean(vals):.1f}%, median={np.median(vals):.1f}%")

    fig, ax = plt.subplots(figsize=(7, 6.5), dpi=300)
    fig.patch.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)

    rng = np.random.default_rng(7)
    for i, (label, vals, color) in enumerate(groups):
        jitter = rng.uniform(-0.16, 0.16, size=len(vals))
        alpha = 0.85 if len(vals) < 60 else 0.35
        size = 34 if len(vals) < 60 else 14
        ax.scatter(np.full(len(vals), i) + jitter, vals, color=color, alpha=alpha, s=size,
                   edgecolors="none", zorder=3)
        mean_val = np.mean(vals)
        ax.hlines(mean_val, i - 0.28, i + 0.28, color=INK, linewidth=2.2, zorder=4)
        ax.text(i, 104, f"mean {mean_val:.1f}%", ha="center", va="bottom", fontsize=9,
                 fontweight="bold", color=INK)

    ax.set_xlim(-0.6, 2.6)
    ax.set_xticks(range(len(groups)))
    ax.set_xticklabels([g[0] for g in groups], fontsize=10.5)
    ax.set_ylim(-3, 112)
    ax.set_ylabel("% of tested experiments where the gene was significantly upregulated", fontsize=9.5, color=MUTED)
    ax.tick_params(axis="y", labelsize=9, colors=MUTED)
    ax.tick_params(axis="x", labelsize=10.5, colors=SECONDARY_INK)
    ax.grid(axis="y", color=GRID, linewidth=0.8, zorder=0)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.spines["left"].set_color("#c3c2b7")
    ax.spines["bottom"].set_color("#c3c2b7")

    fig.suptitle("Per-gene upregulation rate: target genes vs. background",
                  x=0.05, y=0.98, ha="left", fontsize=13, fontweight="bold", color=INK)
    fig.text(0.05, 0.935, "each dot = one gene; black bar = group mean; background = nitrogen only (step 3)",
              ha="left", fontsize=9, color=SECONDARY_INK)

    fig.tight_layout(rect=[0, 0, 1, 0.91])
    out_path = FIG_DIR / "03_pct_scatter.png"
    fig.savefig(out_path, dpi=300, facecolor=SURFACE)
    print(f"\nWrote {out_path}")


if __name__ == "__main__":
    main()
