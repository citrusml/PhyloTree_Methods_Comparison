#!/usr/bin/env python3
"""
True Pairwise Alignment (TRUE_PWA+NJ) Benchmark Plotting Script
Location: bin/plot_true_pwa_benchmark.py

Generates:
1. Heatmap of normalized RF distance (nRF) across Distance D and Sequence Length L.
2. Line plots showing nRF vs Distance D for each Sequence Length L.
3. Summary statistics table of nRF (mean ± std) across all conditions.
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

plt.rcParams.update({
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "font.family": "sans-serif",
    "axes.edgecolor": "#333333",
    "axes.linewidth": 1.0,
    "grid.color": "#e0e0e0",
    "grid.linestyle": "--",
    "grid.alpha": 0.7,
})

def plot_true_pwa_heatmap(df, outdir):
    """Plots heatmap of mean nRF across Distance D and Length L."""
    nrf_col = "nrf_distance" if "nrf_distance" in df.columns else "nrf"
    agg = df.groupby(["distance", "length"])[nrf_col].mean().reset_index()
    pivot = agg.pivot(index="distance", columns="length", values=nrf_col)
    pivot = pivot.sort_index(ascending=False)

    plt.figure(figsize=(8, 6), dpi=300)
    ax = sns.heatmap(
        pivot,
        annot=True,
        fmt=".4f",
        cmap="YlGnBu_r",
        cbar_kws={"label": "Mean Normalized RF Distance (nRF)"},
        linewidths=0.5,
        vmin=0.0,
        vmax=max(0.2, pivot.values.max())
    )
    plt.title("True Pairwise Alignment (TRUE_PWA+NJ) Accuracy Map\n(N=32 Taxa, LG+G4, alpha=1.0)", fontsize=13, fontweight="bold", pad=12)
    plt.xlabel("Sequence Length L (aa)", fontsize=11, fontweight="bold")
    plt.ylabel("Evolutionary Distance D (subst/site)", fontsize=11, fontweight="bold")
    plt.tight_layout()

    out_png = os.path.join(outdir, "true_pwa_heatmap.png")
    out_pdf = os.path.join(outdir, "true_pwa_heatmap.pdf")
    plt.savefig(out_png, dpi=300)
    plt.savefig(out_pdf)
    plt.close()
    print(f"Heatmap saved to {out_png}")

def plot_true_pwa_curves(df, outdir):
    """Plots line curves of nRF vs Distance D for each Length L."""
    nrf_col = "nrf_distance" if "nrf_distance" in df.columns else "nrf"
    agg = df.groupby(["distance", "length"])[nrf_col].agg(["mean", "std", "count"]).reset_index()
    agg["se"] = agg["std"] / np.sqrt(agg["count"])

    lengths = sorted(df["length"].unique())
    palette = sns.color_palette("tab10", n_colors=len(lengths))

    plt.figure(figsize=(8, 6), dpi=300)
    for idx, L in enumerate(lengths):
        sub = agg[agg["length"] == L].sort_values("distance")
        plt.errorbar(
            sub["distance"],
            sub["mean"],
            yerr=sub["se"],
            label=f"L = {L} aa",
            marker="o",
            capsize=4,
            linewidth=2,
            color=palette[idx]
        )

    plt.title("TRUE_PWA+NJ: Error Scaling with Evolutionary Distance\n(N=32 Taxa, LG+G4, alpha=1.0)", fontsize=13, fontweight="bold", pad=12)
    plt.xlabel("Evolutionary Distance D (subst/site)", fontsize=11, fontweight="bold")
    plt.ylabel("Mean Normalized RF Distance (nRF)", fontsize=11, fontweight="bold")
    plt.grid(True)
    plt.ylim(bottom=-0.01)
    plt.legend(title="Sequence Length", frameon=True, facecolor="white", framealpha=0.9)
    plt.tight_layout()

    out_png = os.path.join(outdir, "true_pwa_scaling_curves.png")
    out_pdf = os.path.join(outdir, "true_pwa_scaling_curves.pdf")
    plt.savefig(out_png, dpi=300)
    plt.savefig(out_pdf)
    plt.close()
    print(f"Scaling curves saved to {out_png}")

def export_summary_table(df, outdir):
    """Exports summary pivot table with mean ± std of nRF."""
    nrf_col = "nrf_distance" if "nrf_distance" in df.columns else "nrf"
    agg = df.groupby(["distance", "length"])[nrf_col].agg(
        mean="mean",
        std="std",
        median="median",
        count="count"
    ).reset_index()

    out_csv = os.path.join(outdir, "benchmark_true_pwa_statistics.csv")
    agg.to_csv(out_csv, index=False)
    print(f"Statistics table saved to {out_csv}")

def main():
    parser = argparse.ArgumentParser(description="Plot TRUE_PWA+NJ Benchmark Results")
    parser.add_argument("--csv", required=True, help="Input summary CSV file (benchmark_true_pwa_summary.csv)")
    parser.add_argument("--outdir", default="results/results_true_pwa", help="Output directory for plots")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    df = pd.read_csv(args.csv)

    # Ensure required columns
    required_cols = {"distance", "length"}
    if not required_cols.issubset(df.columns) or ("nrf" not in df.columns and "nrf_distance" not in df.columns):
        raise ValueError(f"CSV missing required columns: {required_cols - set(df.columns)}")

    plot_true_pwa_heatmap(df, args.outdir)
    plot_true_pwa_curves(df, args.outdir)
    export_summary_table(df, args.outdir)

if __name__ == "__main__":
    main()
