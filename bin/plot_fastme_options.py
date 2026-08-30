#!/usr/bin/env python3
"""
FastME Options Benchmark Plotting Script (Experiment 10)
Location: bin/plot_fastme_options.py

Generates:
1. Multi-panel line plots of nRF vs Evolutionary Distance across Sequence Lengths (L=100, 300, 500, 1000, 1500).
2. 2D Regime Map comparing PWA+FastME_SPR, MSA+FastME_LG_G, and MSA+ML accuracy.
3. Summary statistics pivot tables.
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

PIPELINE_COLORS = {
    "PWA+FastME_SPR":   "#17becf",  # Cyan
    "MSA+FastME_LG_G":  "#bcbd22",  # Olive / Yellow-green
    "MSA+ML":           "#2ca02c",  # Green
    "PWA+NJ":           "#1f77b4",  # Blue (if present)
    "MSA+NJ":           "#ff7f0e",  # Orange (if present)
}

PIPELINE_MARKERS = {
    "PWA+FastME_SPR":   "o",
    "MSA+FastME_LG_G":  "s",
    "MSA+ML":           "^",
    "PWA+NJ":           "x",
    "MSA+NJ":           "d",
}

PIPELINE_ORDER = ["PWA+FastME_SPR", "MSA+FastME_LG_G", "MSA+ML"]

def plot_scaling_curves(df, outdir):
    """Plots Normalized RF Distance vs Evolutionary Distance D for each Sequence Length L."""
    nrf_col = "nrf_distance" if "nrf_distance" in df.columns else "nrf"
    lengths = sorted(df["length"].unique())
    distances = sorted(df["distance"].unique())

    ncols = len(lengths)
    fig, axes = plt.subplots(1, ncols, figsize=(4.5 * ncols, 4.5), sharey=True)
    if ncols == 1:
        axes = [axes]

    for idx, l in enumerate(lengths):
        ax = axes[idx]
        sub_df = df[df["length"] == l]

        for pipe in PIPELINE_ORDER:
            pipe_df = sub_df[sub_df["pipeline"] == pipe]
            if pipe_df.empty:
                continue

            agg = pipe_df.groupby("distance")[nrf_col].agg(["mean", "std", "count"]).reset_index()
            agg["se"] = agg["std"] / np.sqrt(agg["count"])

            ax.errorbar(
                agg["distance"],
                agg["mean"],
                yerr=agg["se"],
                label=pipe,
                color=PIPELINE_COLORS.get(pipe, "#333"),
                marker=PIPELINE_MARKERS.get(pipe, "o"),
                capsize=3,
                linewidth=2.0,
                markersize=6
            )

        ax.set_title(f"Length L = {l} aa", fontsize=12, fontweight="bold", pad=8)
        ax.set_xlabel("Evolutionary Distance (D)", fontsize=11, fontweight="bold")
        ax.set_ylim(-0.02, 1.02)
        ax.set_xticks(distances)
        ax.grid(True)

        if idx == 0:
            ax.set_ylabel("Normalized RF Distance", fontsize=11, fontweight="bold")

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(
        handles, labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 1.05),
        ncol=len(PIPELINE_ORDER),
        fontsize=11,
        frameon=True,
        facecolor="white",
        framealpha=0.95
    )

    plt.suptitle("FastME Options Benchmark: Topological Error (nRF) vs Distance D (N=32 Taxa)", fontsize=13, fontweight="bold", y=1.12)
    plt.tight_layout()

    out_png = os.path.join(outdir, "scaling_curves_fastme_options.png")
    out_pdf = os.path.join(outdir, "scaling_curves_fastme_options.pdf")
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.savefig(out_pdf, bbox_inches="tight")
    plt.close()
    print(f"Scaling curves saved to {out_png}")

def plot_regime_map(df, outdir):
    """Plots 2D Regime Map (Best Performing Method per D x L condition)."""
    nrf_col = "nrf_distance" if "nrf_distance" in df.columns else "nrf"
    agg = df.groupby(["distance", "length", "pipeline"])[nrf_col].mean().reset_index()

    idx_min = agg.groupby(["distance", "length"])[nrf_col].idxmin()
    best_df = agg.loc[idx_min].copy()

    pipes_present = [p for p in PIPELINE_ORDER if p in best_df["pipeline"].unique()]
    pipe_to_int = {pipe: i for i, pipe in enumerate(pipes_present)}
    best_df["winner_val"] = best_df["pipeline"].map(pipe_to_int)

    pivot = best_df.pivot(index="distance", columns="length", values="winner_val")

    fig, ax = plt.subplots(figsize=(7, 5))
    cmap = sns.color_palette([PIPELINE_COLORS[p] for p in pipes_present])

    sns.heatmap(
        pivot,
        annot=True,
        cmap=cmap,
        cbar=False,
        ax=ax,
        linewidths=1.0,
        linecolor="white",
        fmt=".0f"
    )

    ax.set_title("FastME Options Regime Map (Lowest Mean nRF)", fontsize=13, fontweight="bold", pad=12)
    ax.set_xlabel("Sequence Length (L)", fontsize=11, fontweight="bold")
    ax.set_ylabel("Evolutionary Distance (D)", fontsize=11, fontweight="bold")
    ax.invert_yaxis()

    out_png = os.path.join(outdir, "regime_map_fastme_options.png")
    out_pdf = os.path.join(outdir, "regime_map_fastme_options.pdf")
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.savefig(out_pdf, bbox_inches="tight")
    plt.close()
    print(f"Regime map saved to {out_png}")

def export_statistics(df, outdir):
    """Exports summary statistics table."""
    nrf_col = "nrf_distance" if "nrf_distance" in df.columns else "nrf"
    group_cols = ["pipeline", "distance", "length"]
    agg = df.groupby(group_cols)[nrf_col].agg(
        mean="mean",
        std="std",
        median="median",
        count="count"
    ).reset_index()

    out_csv = os.path.join(outdir, "benchmark_fastme_options_statistics.csv")
    agg.to_csv(out_csv, index=False)
    print(f"Summary statistics table saved to {out_csv}")

def main():
    parser = argparse.ArgumentParser(description="Plot FastME Options Benchmark Results (Experiment 10)")
    parser.add_argument("--csv", required=True, help="Input summary CSV file")
    parser.add_argument("--outdir", default="results/results_fastme_options", help="Output directory for plots")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    df = pd.read_csv(args.csv)

    plot_scaling_curves(df, args.outdir)
    export_statistics(df, args.outdir)
    try:
        plot_regime_map(df, args.outdir)
    except Exception as e:
        print(f"Note: Could not generate regime map: {e}")

if __name__ == "__main__":
    main()
