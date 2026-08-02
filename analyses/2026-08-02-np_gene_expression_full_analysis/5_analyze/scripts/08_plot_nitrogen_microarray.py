"""Microarray-specific figure for the same 53-gene x nitrogen-starvation
data as 07_plot_nitrogen_volcano.py, but a different chart type: a
volcano plot's y-axis needs a continuous, graded significance measure,
and this KG build's microarray padj is not that -- it's effectively
binary (0.01 = "significant" for 90 of 546 rows, 1.0 = "not
significant" for the other 456; no values in between). Plotting that on
a continuous -log10(padj) axis produces two flat, uninformative bands
(see 5_analyze/notebook.md).

Instead: a strip plot, x = log2FC, y = the two real categories the data
actually distinguishes ("Significant" / "Not significant"), colored by
each gene's original N/P annotation, points jittered within each row.
This shows exactly what the microarray data supports -- direction and
magnitude of fold-change, and a binary significance call -- without
implying false precision.

Input:  data/06_nitrogen_volcano_data.csv
Output: figures/05_nitrogen_microarray.png

Usage:
  uv run python analyses/2026-08-02-np_gene_expression_full_analysis/5_analyze/scripts/08_plot_nitrogen_microarray.py
"""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

STEP_DIR = Path(__file__).resolve().parents[1]
DATA_DIR = STEP_DIR / "data"
FIG_DIR = STEP_DIR / "figures"

N_COLOR = "#2a78d6"   # blue -- categorical slot 1, same as volcano plot
P_COLOR = "#eb6834"   # orange -- categorical slot 8, same as volcano plot
INK = "#0b0b0b"
MUTED = "#898781"
GRID = "#e1e0d9"

ROW_ORDER = ["Not significant", "Significant"]
RNG = np.random.default_rng(42)


def main() -> None:
    df = pd.read_csv(DATA_DIR / "06_nitrogen_volcano_data.csv")
    ma = df[df["omics_type"] == "MICROARRAY"].copy()
    ma["sig_label"] = np.where(ma["padj"] < 0.05, "Significant", "Not significant")
    ma["y"] = ma["sig_label"].map({"Not significant": 0, "Significant": 1})
    ma["y_jitter"] = ma["y"] + RNG.uniform(-0.18, 0.18, size=len(ma))

    fig, ax = plt.subplots(figsize=(10, 5), dpi=300)

    for label, color in [("N", N_COLOR), ("P", P_COLOR)]:
        g = ma[ma["n_or_p"] == label]
        ax.scatter(
            g["log2fc"], g["y_jitter"], s=18, alpha=0.6, color=color,
            edgecolors="none", label=f"{label}-annotated (n={len(g)})",
        )

    ax.axvline(0, color=MUTED, linewidth=0.8, linestyle=":")
    ax.set_yticks([0, 1])
    ax.set_yticklabels(ROW_ORDER, fontsize=10, color=INK)
    ax.set_ylim(-0.6, 1.6)
    ax.set_xlabel("log2 fold change", color=INK)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.grid(axis="x", color=GRID, linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)
    ax.tick_params(colors=MUTED)
    ax.legend(loc="upper left", frameon=False, fontsize=9)

    n_sig = int((ma["padj"] < 0.05).sum())
    ax.set_title(
        "Nitrogen-starvation DE, microarray only, all 53 genes -- "
        f"{n_sig}/{len(ma)} rows significant\n"
        "padj is binary in this KG build (0.01 vs 1.0, not continuous) -- shown as two rows, not a volcano y-axis.",
        fontsize=9.5, color=INK, loc="left",
    )
    fig.tight_layout()

    out_path = FIG_DIR / "05_nitrogen_microarray.png"
    fig.savefig(out_path, dpi=300)
    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
