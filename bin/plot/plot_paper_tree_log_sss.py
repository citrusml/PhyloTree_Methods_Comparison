#!/usr/bin/env python3
"""
Plot Log-Scale SSS vs nRF and Accuracy for results_paper_tree
Replicating Matsui & Iwasaki (2020) Fig. 3a with logarithmic horizontal axis.
"""

import os
import shutil
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from statsmodels.nonparametric.smoothers_lowess import lowess

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
PAPER_TREE_DIR = PROJECT_ROOT / "results" / "results_paper_tree"
CSV_PATH = PAPER_TREE_DIR / "benchmark_summary.csv"

# Consistent styling matching Matsui & Iwasaki (2020) and project guidelines
plt.rcParams.update({
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "font.family": "sans-serif",
    "axes.edgecolor": "#333333",
    "axes.linewidth": 1.0,
    "grid.color": "#e0e0e0",
    "grid.linestyle": "--",
    "grid.alpha": 0.6,
    "savefig.dpi": 300,
    "figure.dpi": 300,
})

PIPELINE_COLORS = {
    "PWA+NJ": "#1f77b4",       # Blue
    "MSA+NJ": "#ff7f0e",       # Orange
    "MSA+ML": "#2ca02c",       # Green (IQ-TREE 2)
    "MSA+RAXML": "#006400",    # Dark Green (RAxML -f d)
    "GS": "#9467bd",           # Purple (Matsui & Iwasaki 2020 Graph Splitting)
    "TRUE_PWA+NJ": "#17becf",  # Cyan
    "TRUE_MSA+NJ": "#8c564b",  # Brown
    "TRUE_MSA+ML": "#e377c2",  # Pink (True IQ-TREE 2)
    "TRUE_MSA+RAXML": "#c49c94", # Rosy brown (True RAxML -f d)
    "TRUE_DIST+NJ": "#7f7f7f"  # Grey
}

preferred_order = [
    "PWA+NJ", "MSA+NJ", "MSA+ML", "MSA+RAXML", "GS",
    "TRUE_PWA+NJ", "TRUE_MSA+NJ", "TRUE_MSA+ML", "TRUE_MSA+RAXML", "TRUE_DIST+NJ"
]


def generate_log_sss_plots():
    if not CSV_PATH.exists():
        raise FileNotFoundError(f"Data file {CSV_PATH} not found.")

    df = pd.read_csv(CSV_PATH)
    pipes = [p for p in preferred_order if p in df["pipeline"].unique()] + [
        p for p in df["pipeline"].unique() if p not in preferred_order
    ]

    fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(16, 7))

    for pipe in pipes:
        sub = df[(df["pipeline"] == pipe) & df["sss_mean"].notna() & df["nrf_distance"].notna()].copy()
        if len(sub) < 3:
            continue
        sub = sub.sort_values("sss_mean")
        x = sub["sss_mean"].values
        y_err = sub["nrf_distance"].values
        y_acc = 1.0 - y_err

        color = PIPELINE_COLORS.get(pipe, "#333333")
        is_ctrl = "TRUE" in pipe
        lw_width = 2.0 if is_ctrl else 2.5
        scatter_alpha = 0.06 if is_ctrl else 0.12

        # LOWESS smoothing across SSS
        sm_err = lowess(y_err, x, frac=0.35, it=3)
        sm_acc = lowess(y_acc, x, frac=0.35, it=3)

        z_line = 3 if is_ctrl else 5
        z_scat = 2 if is_ctrl else 4

        ax1.plot(sm_err[:, 0], sm_err[:, 1], color=color, linewidth=lw_width, label=pipe, zorder=z_line)
        ax2.plot(sm_acc[:, 0], sm_acc[:, 1], color=color, linewidth=lw_width, label=pipe, zorder=z_line)

        ax1.scatter(x, y_err, color=color, alpha=scatter_alpha, s=12, edgecolors="none", zorder=z_scat)
        ax2.scatter(x, y_acc, color=color, alpha=scatter_alpha, s=12, edgecolors="none", zorder=z_scat)

    # X-axis ticks for log scale
    xticks_major = [0.015, 0.02, 0.03, 0.04, 0.06, 0.10, 0.15, 0.20, 0.30, 0.50]
    xtick_labels = ["0.015", "0.02", "0.03", "0.04", "0.06", "0.10", "0.15", "0.20", "0.30", "0.50"]

    for ax, title, ylabel, ylim, loc in [
        (
            ax1,
            "Topological Error (nRF) vs Sequence Similarity Score (SSS)\n[Log Scale]",
            "Normalized RF Distance (Lower is Better)",
            (-0.02, 0.72),
            "upper right",
        ),
        (
            ax2,
            "Topological Accuracy (1 - nRF) vs Sequence Similarity Score (SSS)\n[Matsui & Iwasaki (2020) Fig. 3a Replication, Log Scale]",
            "Correct Topology Ratio (1 - nRF, Higher is Better)",
            (0.28, 1.02),
            "lower right",
        ),
    ]:
        # Vertical reference lines for biological regimes
        ax.axvline(
            0.06,
            color="#d62728",
            linestyle="--",
            linewidth=1.5,
            alpha=0.85,
            label=r"Crossover Threshold ($SSS = 0.06$)",
            zorder=2,
        )
        ax.axvline(
            0.03,
            color="#8b0000",
            linestyle=":",
            linewidth=1.8,
            alpha=0.85,
            label=r"Extreme Breakdown ($SSS \leq 0.03$)",
            zorder=2,
        )

        ax.set_xscale("log")
        ax.set_xlim(0.012, 0.55)
        ax.set_xticks(xticks_major)
        ax.set_xticklabels(xtick_labels, fontsize=9)
        ax.set_ylim(ylim)

        ax.set_title(title, fontsize=12, fontweight="bold", pad=12)
        ax.set_xlabel(r"Average Sequence Similarity Score ($\bar{w}$, log scale)", fontsize=11, fontweight="bold")
        ax.set_ylabel(ylabel, fontsize=11, fontweight="bold")
        ax.grid(True, which="both", linestyle="--", alpha=0.45)
        ax.legend(fontsize=9, loc=loc, framealpha=0.92, facecolor="white", edgecolor="#cccccc")

    plt.tight_layout()

    # Backup existing linear image if not backed up yet
    orig_path = PAPER_TREE_DIR / "nrf_vs_sss_curves.png"
    linear_path = PAPER_TREE_DIR / "nrf_vs_sss_linear.png"
    if orig_path.exists() and not linear_path.exists():
        shutil.copy(orig_path, linear_path)
        print(f"Preserved linear plot to: {linear_path}")

    # Save log scale plots
    log_path = PAPER_TREE_DIR / "nrf_vs_sss_curves_log.png"
    plt.savefig(log_path, dpi=300)
    print(f"Saved: {log_path}")

    # Also save to nrf_vs_sss_curves.png as requested by user
    plt.savefig(orig_path, dpi=300)
    print(f"Updated: {orig_path}")
    plt.close()


if __name__ == "__main__":
    generate_log_sss_plots()
