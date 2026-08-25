#!/usr/bin/env python3
"""
Gamma Distance Benchmark Plotting Script
Location: bin/plot_gamma_benchmark.py

Generates:
1. Multi-panel line plots of nRF vs Simulation Alpha across Distance D and Sequence Length L
   using Gamma-corrected evolutionary distance (dist_model=gamma_poisson, alpha=1.0).
2. Direct comparison between PWA+NJ, MSA+NJ, and MSA+ML under gamma distance correction.
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
    "PWA+NJ": "#1f77b4",  # Blue
    "MSA+NJ": "#ff7f0e",  # Orange
    "MSA+ML": "#2ca02c",  # Green
}

PIPELINE_MARKERS = {
    "PWA+NJ": "o",
    "MSA+NJ": "s",
    "MSA+ML": "^",
}

PIPELINE_ORDER = ["PWA+NJ", "MSA+NJ", "MSA+ML"]

def plot_gamma_alpha_scaling(df, outdir):
    """Plots nRF vs Simulation Alpha across Distance (rows) and Length (cols)."""
    if "alpha" not in df.columns:
        print("Warning: 'alpha' column not found in dataframe.")
        return

    nrf_col = "nrf_distance" if "nrf_distance" in df.columns else "nrf"
    distances = sorted(df["distance"].unique())
    lengths = sorted(df["length"].unique())
    alphas = sorted(df["alpha"].unique())

    fig, axes = plt.subplots(
        nrows=len(distances),
        ncols=len(lengths),
        figsize=(16, 14),
        sharex=True,
        sharey=True
    )

    for r, d in enumerate(distances):
        for c, l in enumerate(lengths):
            ax = axes[r, c] if len(distances) > 1 and len(lengths) > 1 else (axes[c] if len(distances) == 1 else axes[r])
            sub_df = df[(df["distance"] == d) & (df["length"] == l)]

            for pipe in PIPELINE_ORDER:
                pipe_df = sub_df[sub_df["pipeline"] == pipe]
                if pipe_df.empty:
                    continue

                agg = pipe_df.groupby("alpha")[nrf_col].agg(["mean", "std", "count"]).reset_index()
                agg["se"] = agg["std"] / np.sqrt(agg["count"])

                ax.errorbar(
                    agg["alpha"],
                    agg["mean"],
                    yerr=agg["se"],
                    label=pipe,
                    color=PIPELINE_COLORS.get(pipe, "#333"),
                    marker=PIPELINE_MARKERS.get(pipe, "o"),
                    capsize=3,
                    linewidth=1.8,
                    markersize=6
                )

            ax.set_title(f"D = {d}, L = {l} aa", fontsize=11, fontweight="bold", pad=6)
            ax.grid(True)
            ax.set_ylim(-0.02, 1.02)
            ax.set_xscale("log", base=2)
            ax.set_xticks(alphas)
            ax.get_xaxis().set_major_formatter(plt.ScalarFormatter())

            if c == 0:
                ax.set_ylabel("Normalized RF Distance", fontsize=10, fontweight="bold")
            if r == len(distances) - 1:
                ax.set_xlabel("Simulation Alpha (Rate Heterogeneity)", fontsize=10, fontweight="bold")

    handles, labels = axes[0, 0].get_legend_handles_labels() if len(distances) > 1 and len(lengths) > 1 else axes[0].get_legend_handles_labels()
    fig.legend(
        handles, labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.99),
        ncol=len(PIPELINE_ORDER),
        fontsize=12,
        frameon=True,
        facecolor="white",
        framealpha=0.95
    )

    plt.suptitle("Gamma Distance Correction: Topological Accuracy vs Rate Heterogeneity\n(Distance Formula: Gamma-Poisson with alpha=1.0, N=32 Taxa)", fontsize=14, fontweight="bold", y=1.02)
    plt.tight_layout()

    out_png = os.path.join(outdir, "gamma_distance_alpha_scaling.png")
    out_pdf = os.path.join(outdir, "gamma_distance_alpha_scaling.pdf")
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.savefig(out_pdf, bbox_inches="tight")
    plt.close()
    print(f"Gamma scaling plot saved to {out_png}")

def export_gamma_statistics(df, outdir):
    """Exports summary statistics table for gamma benchmark."""
    nrf_col = "nrf_distance" if "nrf_distance" in df.columns else "nrf"
    group_cols = ["pipeline", "alpha", "distance", "length"]
    agg = df.groupby(group_cols)[nrf_col].agg(
        mean="mean",
        std="std",
        median="median",
        count="count"
    ).reset_index()

    out_csv = os.path.join(outdir, "benchmark_gamma_statistics.csv")
    agg.to_csv(out_csv, index=False)
    print(f"Summary statistics table saved to {out_csv}")

def main():
    parser = argparse.ArgumentParser(description="Plot Gamma Distance Benchmark Results")
    parser.add_argument("--csv", required=True, help="Input summary CSV file (benchmark_gamma_summary.csv)")
    parser.add_argument("--outdir", default="results/results_gamma", help="Output directory for plots")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    df = pd.read_csv(args.csv)

    plot_gamma_alpha_scaling(df, args.outdir)
    export_gamma_statistics(df, args.outdir)

if __name__ == "__main__":
    main()
