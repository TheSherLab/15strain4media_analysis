"""
Step 5 -- percentage of experiments each gene was found upregulated in.

For every gene: among the experiments where it was actually tested
(excludes "no data" cells), what fraction came back significant_up?
Computed against ALL 10 experiments per gene (not split by matched/cross),
so a gene's percentage reflects its total observed upregulation rate
across every experiment it had data in, whichever nutrient that
experiment was.

Inputs: data/01_target_gene_experiment_matrix.csv
Outputs: data/04_pct_upregulated.csv (all 61 genes)
         figures/02_pct_upregulated.png (genes with >=1 tested experiment)

Usage: uv run analyses/2026-08-10-np_starvation_expression_walkthrough/5_analyze/scripts/04_pct_upregulated.py
"""

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

BASE = Path(__file__).resolve().parent.parent
DATA_DIR = BASE / "data"
FIG_DIR = BASE / "figures"

TESTED_STATUSES = ["significant_up", "significant_down", "not_significant"]
COLOR_N = "#2a78d6"
COLOR_P = "#e34948"
INK = "#0b0b0b"
SECONDARY_INK = "#52514e"
MUTED = "#898781"
SURFACE = "#fcfcfb"
GRID = "#e1e0d9"


def main() -> None:
    matrix = pd.read_csv(DATA_DIR / "01_target_gene_experiment_matrix.csv")
    tested = matrix[matrix["status"].isin(TESTED_STATUSES)]

    summary = tested.groupby(["gene_name", "n_or_p"]).agg(
        n_tested=("status", "size"),
        n_up=("status", lambda s: (s == "significant_up").sum()),
        n_down=("status", lambda s: (s == "significant_down").sum()),
    ).reset_index()
    summary["pct_up"] = 100 * summary["n_up"] / summary["n_tested"]
    summary["pct_down"] = 100 * summary["n_down"] / summary["n_tested"]
    summary = summary.sort_values(["n_or_p", "pct_up"], ascending=[False, False])

    out_path = DATA_DIR / "04_pct_upregulated.csv"
    summary.to_csv(out_path, index=False)
    print(f"Wrote {len(summary)} rows to {out_path}")
    print(summary.to_string(index=False))

    # Figure: only genes tested in >=1 experiment (all 61 named genes qualify
    # except the 1 with zero locus/data everywhere).
    plot_df = summary[summary["n_tested"] > 0].copy()
    plot_df = plot_df.sort_values(["n_or_p", "pct_up"], ascending=[True, True])

    fig_height = 0.185 * len(plot_df) + 1.7
    fig, ax = plt.subplots(figsize=(8.8, fig_height), dpi=300)
    fig.patch.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)

    colors = [COLOR_N if r.n_or_p == "N" else COLOR_P for r in plot_df.itertuples()]
    y_pos = range(len(plot_df))
    ax.barh(list(y_pos), plot_df["pct_up"], color=colors, height=0.7, zorder=3)

    labels = [f"{r.gene_name} (n={r.n_tested})" for r in plot_df.itertuples()]
    ax.set_yticks(list(y_pos))
    ax.set_yticklabels(labels, fontsize=7, fontfamily="monospace")
    for tick, r in zip(ax.get_yticklabels(), plot_df.itertuples()):
        tick.set_color(COLOR_N if r.n_or_p == "N" else COLOR_P)

    ax.set_xlim(0, 105)
    ax.set_xlabel("% of tested experiments where the gene was significantly upregulated", fontsize=9, color=MUTED)
    ax.tick_params(axis="x", labelsize=8, colors=MUTED)
    ax.grid(axis="x", color=GRID, linewidth=0.8, zorder=0)
    for spine in ("top", "right", "left"):
        ax.spines[spine].set_visible(False)
    ax.spines["bottom"].set_color("#c3c2b7")

    fig.suptitle("How often each gene was upregulated, of the experiments it was tested in",
                  x=0.01, y=0.995, ha="left", va="top", fontsize=11.5, fontweight="bold", color=INK)
    fig.text(0.01, 0.995 - 0.28 / fig_height, "blue label = N-annotated gene, red label = P-annotated gene; "
                            "n = number of experiments the gene was actually tested in",
              ha="left", va="top", fontsize=8, color=SECONDARY_INK)

    fig.tight_layout(rect=[0, 0, 1, 1 - 0.55 / fig_height])
    out_fig = FIG_DIR / "02_pct_upregulated.png"
    fig.savefig(out_fig, dpi=300, facecolor=SURFACE)
    print(f"\nWrote {out_fig}")


if __name__ == "__main__":
    main()
