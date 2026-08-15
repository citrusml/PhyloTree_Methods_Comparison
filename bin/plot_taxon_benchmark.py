#!/usr/bin/env python3
"""
Taxon Scaling Benchmark Plotting Script
Location: bin/plot_taxon_benchmark.py

Generates:
1. Line plots showing nRF vs Taxon Count (N) for each (Distance, Length) condition.
2. Comprehensive multi-panel plots comparing PWA+NJ, MSA+NJ, and MSA+ML scaling behavior.
3. Summary pivot table of nRF (mean ± std) across Taxon Counts.
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

def plot_taxon_scaling_lines(df, outdir):
    """
    Plots nRF vs Taxon Count N (8, 32, 64, 128) across Distance (rows) and Length (cols).
    """
    distances = sorted(df["distance"].unique())
    lengths = sorted(df["length"].unique())
    taxa = sorted(df["num_taxa"].unique())

    fig, axes = plt.subplots(
        nrows=len(distances),
        ncols=len(lengths),
        figsize=(18, 16),
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

                stats = pipe_df.groupby("num_taxa")["nrf_distance"].agg(["mean", "std", "count"]).reset_index()
                ax.plot(
                    stats["num_taxa"],
                    stats["mean"],
                    marker=PIPELINE_MARKERS.get(pipe, "o"),
                    markersize=6,
                    linewidth=2.0,
                    color=PIPELINE_COLORS.get(pipe, "#888"),
                    label=pipe
                )
                # Error band (std)
                ax.fill_between(
                    stats["num_taxa"],
                    np.maximum(0, stats["mean"] - stats["std"]),
                    np.minimum(1, stats["mean"] + stats["std"]),
                    color=PIPELINE_COLORS.get(pipe, "#888"),
                    alpha=0.15
                )

            ax.grid(True, linestyle="--", alpha=0.6)
            ax.set_ylim(-0.02, 0.65)
            ax.set_xticks(taxa)

            if r == 0:
                ax.set_title(f"Length L = {l} aa", fontsize=12, fontweight="bold", pad=8)
            if c == 0:
                ax.set_ylabel(f"D = {d}\nNormalized RF", fontsize=11, fontweight="bold")
            if r == len(distances) - 1:
                ax.set_xlabel("Number of Taxa (N)", fontsize=11, fontweight="bold")

    handles, labels = axes[0, 0].get_legend_handles_labels() if len(distances) > 1 and len(lengths) > 1 else axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 0.995), ncol=3, fontsize=13, frameon=True)

    plt.suptitle("Phylogenetic Accuracy (Normalized RF) vs Taxon Count (N = 8, 32, 64, 128)\nacross Evolutionary Distance (D) and Sequence Length (L)", fontsize=15, fontweight="bold", y=1.02)
    plt.tight_layout()

    out_png = os.path.join(outdir, "taxon_scaling_curves.png")
    out_pdf = os.path.join(outdir, "taxon_scaling_curves.pdf")
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.savefig(out_pdf, bbox_inches="tight")
    plt.close()
    print(f"Saved taxon scaling curve plot: {out_png}")
    print(f"Saved taxon scaling curve plot (PDF): {out_pdf}")

def generate_taxon_summary_tables(df, outdir):
    """Generates detailed CSV summary and pivot tables across Taxon counts."""
    grouped = df.groupby(["num_taxa", "distance", "length", "pipeline"])["nrf_distance"]
    stats_df = grouped.agg(
        replicates="count",
        nrf_mean="mean",
        nrf_variance=lambda x: np.var(x, ddof=1) if len(x) > 1 else 0.0,
        nrf_std=lambda x: np.std(x, ddof=1) if len(x) > 1 else 0.0,
        nrf_median="median"
    ).reset_index()

    stats_df["nrf_mean"] = stats_df["nrf_mean"].round(5)
    stats_df["nrf_variance"] = stats_df["nrf_variance"].round(6)
    stats_df["nrf_std"] = stats_df["nrf_std"].round(5)

    out_csv = os.path.join(outdir, "benchmark_taxon_statistics.csv")
    stats_df.to_csv(out_csv, index=False)
    print(f"Saved taxon statistics summary: {out_csv}")

def main():
    parser = argparse.ArgumentParser(description="Plot Taxon Benchmark scaling results (N=8, 32, 64, 128)")
    parser.add_argument("--csv", required=True, help="Path to benchmark_taxon_summary.csv")
    parser.add_argument("--outdir", default=".", help="Output directory for plots and summary tables")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    df = pd.read_csv(args.csv)

    print(f"=== Taxon Benchmark Plotting ===")
    print(f"Loaded {len(df)} records across taxa: {sorted(df['num_taxa'].unique())}")

    plot_taxon_scaling_lines(df, args.outdir)
    generate_taxon_summary_tables(df, args.outdir)

    print("Taxon plotting completed successfully!")

if __name__ == "__main__":
    main()
