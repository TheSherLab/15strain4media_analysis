"""One figure per nutrient (nitrogen, phosphorus), each with 3 panels
(RNASEQ, PROTEOMICS, MICROARRAY) side by side -- combines
07_plot_nitrogen_volcano.py and 08_plot_nitrogen_microarray.py's
approaches into a single figure per nutrient instead of two separate
files, per the researcher's request.

RNASEQ/PROTEOMICS panels are true volcano plots (log2FC vs
-log10(padj)) -- both platforms have continuous padj. The MICROARRAY
panel is a strip plot (log2FC vs "Significant"/"Not significant") --
this KG build's microarray padj is effectively binary (see
5_analyze/notebook.md), so a continuous y-axis there would imply
confidence gradation the data doesn't have. All panels colored by each
gene's original N/P annotation, with a legend.

Nitrogen figure uses data/06_nitrogen_volcano_data.csv (already
pulled). Phosphorus figure uses data/11_phosphorus_de_data.csv --
its N-annotated side is very sparse (2 genes, 6 rows total) because
phosphorus experiment tables are narrow and pre-filtered to
already-significant genes (see 2_kg_selection/notebook.md); shown as-is,
not padded or hidden.

Inputs:  data/06_nitrogen_volcano_data.csv, data/11_phosphorus_de_data.csv
Outputs: figures/07_nitrogen_three_panel.png, figures/08_phosphorus_three_panel.png

Usage:
  uv run python analyses/2026-08-02-np_gene_expression_full_analysis/5_analyze/scripts/12_plot_nutrient_three_panel.py
"""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

STEP_DIR = Path(__file__).resolve().parents[1]
DATA_DIR = STEP_DIR / "data"
FIG_DIR = STEP_DIR / "figures"

N_COLOR = "#2a78d6"
P_COLOR = "#eb6834"
INK = "#0b0b0b"
MUTED = "#898781"
GRID = "#e1e0d9"
RNG = np.random.default_rng(42)


def plot_volcano_panel(ax, sub: pd.DataFrame, platform: str, padj_floor: float) -> None:
    sub = sub.copy()
    sub["padj_floored"] = sub["padj"].clip(lower=padj_floor)
    sub["neg_log10_padj"] = -np.log10(sub["padj_floored"])

    for label, color in [("N", N_COLOR), ("P", P_COLOR)]:
        g = sub[sub["n_or_p"] == label]
        ax.scatter(
            g["log2fc"], g["neg_log10_padj"], s=18, alpha=0.65,
            color=color, edgecolors="none", label=f"{label}-annotated (n={len(g)})",
        )

    ax.axhline(-np.log10(0.05), color=MUTED, linewidth=1, linestyle="--")
    ax.axvline(0, color=MUTED, linewidth=0.8, linestyle=":")

    sig = sub[sub["padj"] < 0.05]
    top = sig.reindex(sig["log2fc"].abs().sort_values(ascending=False).index).head(3)
    for i, (_, r) in enumerate(top.iterrows()):
        ax.annotate(
            r["gene_name_raw"], xy=(r["log2fc"], r["neg_log10_padj"]),
            xytext=(4, 3 + 11 * i), textcoords="offset points", fontsize=7.5, color=INK,
        )

    ax.set_ylim(bottom=0)  # always anchor at 0 -- auto-scaling on tiny n (e.g. phosphorus
    # proteomics, 5 points) zooms into a razor-thin range and misrepresents how close
    # to "not significant" those points actually are.
    ax.set_title(platform, fontsize=10, color=INK)
    ax.set_xlabel("log2 fold change", color=INK)
    ax.set_ylabel("-log10(padj)", color=INK)
    ax.legend(loc="lower left", frameon=False, fontsize=7.5)


def plot_microarray_panel(ax, sub: pd.DataFrame) -> None:
    sub = sub.copy()
    sub["sig_label"] = np.where(sub["padj"] < 0.05, "Significant", "Not significant")
    sub["y"] = sub["sig_label"].map({"Not significant": 0, "Significant": 1})
    sub["y_jitter"] = sub["y"] + RNG.uniform(-0.18, 0.18, size=len(sub))

    # Clip x-axis against known extreme outliers (e.g. PMM0708/phoA, log2FC ~162,
    # MED4 phosphate microarray -- see 5_analyze/notebook.md) so a single point
    # doesn't compress the rest of the panel into an unreadable sliver.
    x_lo, x_hi = -10, 20
    n_clipped = int(((sub["log2fc"] < x_lo) | (sub["log2fc"] > x_hi)).sum())

    for label, color in [("N", N_COLOR), ("P", P_COLOR)]:
        g = sub[sub["n_or_p"] == label]
        ax.scatter(
            g["log2fc"], g["y_jitter"], s=18, alpha=0.6, color=color,
            edgecolors="none", label=f"{label}-annotated (n={len(g)})",
        )

    ax.axvline(0, color=MUTED, linewidth=0.8, linestyle=":")
    ax.set_xlim(x_lo, x_hi)
    ax.set_yticks([0, 1])
    ax.set_yticklabels(["Not\nsignificant", "Significant"], fontsize=9, color=INK)
    ax.set_ylim(-0.6, 1.6)
    ax.set_xlabel("log2 fold change", color=INK)
    n_sig = int((sub["padj"] < 0.05).sum())
    title = f"MICROARRAY (binary padj, {n_sig}/{len(sub)} significant)"
    if n_clipped:
        title += f"\n{n_clipped} known outlier(s) off-scale, see notebook"
    ax.set_title(title, fontsize=10, color=INK)
    ax.legend(loc="upper left", frameon=False, fontsize=7.5)


def make_figure(df: pd.DataFrame, nutrient_label: str, out_path: Path, extra_caveat: str = "") -> None:
    continuous = df[df["omics_type"].isin(["RNASEQ", "PROTEOMICS"])]
    n_zero = int((continuous["padj"] == 0).sum())
    nonzero = continuous.loc[continuous["padj"] > 0, "padj"]
    padj_floor = nonzero.min() / 10 if not nonzero.empty else 1e-12

    fig, axes = plt.subplots(1, 3, figsize=(15.5, 5.5), dpi=300)

    rnaseq = df[df["omics_type"] == "RNASEQ"]
    proteomics = df[df["omics_type"] == "PROTEOMICS"]
    microarray = df[df["omics_type"] == "MICROARRAY"]

    if not rnaseq.empty:
        plot_volcano_panel(axes[0], rnaseq, "RNASEQ", padj_floor)
    else:
        axes[0].set_title("RNASEQ (no data)", fontsize=10, color=INK)
    if not proteomics.empty:
        plot_volcano_panel(axes[1], proteomics, "PROTEOMICS", padj_floor)
    else:
        axes[1].set_title("PROTEOMICS (no data)", fontsize=10, color=INK)
    if not microarray.empty:
        plot_microarray_panel(axes[2], microarray)
    else:
        axes[2].set_title("MICROARRAY (no data)", fontsize=10, color=INK)

    for ax in axes:
        for spine in ("top", "right"):
            ax.spines[spine].set_visible(False)
        ax.grid(color=GRID, linewidth=0.8, zorder=0)
        ax.set_axisbelow(True)
        ax.tick_params(colors=MUTED)

    subtitle = (
        f"{nutrient_label}-starvation DE, all 53 genes, colored by original N/P annotation. "
        f"RNASEQ/PROTEOMICS: volcano ({n_zero} padj=0 points floored to {padj_floor:.1e}). "
        "MICROARRAY: strip plot (binary padj here, not continuous -- see notebook)."
    )
    if extra_caveat:
        subtitle += f"\n{extra_caveat}"
    fig.suptitle(subtitle, fontsize=9, color=INK, x=0.01, y=0.99, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.87 if extra_caveat else 0.90))

    fig.savefig(out_path, dpi=300)
    print(f"Wrote {out_path}")


def main() -> None:
    nitrogen = pd.read_csv(DATA_DIR / "06_nitrogen_volcano_data.csv")
    phosphorus = pd.read_csv(DATA_DIR / "11_phosphorus_de_data.csv")

    make_figure(nitrogen, "Nitrogen", FIG_DIR / "07_nitrogen_three_panel.png")
    make_figure(
        phosphorus, "Phosphorus", FIG_DIR / "08_phosphorus_three_panel.png",
        extra_caveat=(
            "N-annotated side is sparse here (2 genes, 6 rows total) -- phosphorus experiment "
            "tables are narrow and pre-filtered to already-significant genes."
        ),
    )


if __name__ == "__main__":
    main()
