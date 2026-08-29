#!/usr/bin/env python3
"""
Alpha Scaling Benchmark Plotting Script
Location: bin/plot_alpha_benchmark.py

Generates:
1. Line plots showing nRF vs Gamma Shape Alpha (0.25, 0.5, 1.0, 2.0) across (Distance, Length) conditions.
2. Distance response curves stratified by Alpha for each pipeline.
3. Summary pivot tables of nRF (mean ± std) across Alpha values.
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

def plot_alpha_scaling_lines(df, outdir):
    """
    Plots nRF vs Gamma Shape Parameter Alpha across Distance (rows) and Length (cols).
    """
    if "alpha" not in df.columns:
        print("Warning: 'alpha' column not found in dataframe. Skipping alpha line plot.")
        return

    distances = sorted(df["distance"].unique())
    lengths = sorted(df["length"].unique())
    alphas = sorted(df["alpha"].unique())

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

                stats = pipe_df.groupby("alpha")["nrf_distance"].agg(["mean", "std", "count"]).reset_index()
                ax.plot(
                    stats["alpha"],
                    stats["mean"],
                    marker=PIPELINE_MARKERS.get(pipe, "o"),
                    markersize=6,
                    linewidth=2.0,
                    color=PIPELINE_COLORS.get(pipe, "#888"),
                    label=pipe
                )
                # Error band (std)
                ax.fill_between(
                    stats["alpha"],
                    np.maximum(0, stats["mean"] - stats["std"]),
                    np.minimum(1, stats["mean"] + stats["std"]),
                    color=PIPELINE_COLORS.get(pipe, "#888"),
                    alpha=0.15
                )

            ax.grid(True, linestyle="--", alpha=0.6)
            ax.set_ylim(-0.02, 0.65)
            ax.set_xscale("log", base=2)
            ax.set_xticks(alphas)
            ax.get_xaxis().set_major_formatter(plt.ScalarFormatter())

            if r == 0:
                ax.set_title(f"Length L = {l} aa", fontsize=12, fontweight="bold", pad=8)
            if c == 0:
                ax.set_ylabel(f"D = {d}\nNormalized RF", fontsize=11, fontweight="bold")
            if r == len(distances) - 1:
                ax.set_xlabel(r"Gamma Shape Parameter $\alpha$", fontsize=11, fontweight="bold")

    handles, labels = axes[0, 0].get_legend_handles_labels() if len(distances) > 1 and len(lengths) > 1 else axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 0.995), ncol=3, fontsize=13, frameon=True)

    plt.suptitle(r"Phylogenetic Accuracy (Normalized RF) vs Rate Heterogeneity ($\alpha \in \{0.25, 0.5, 1.0, 2.0\}$)" + "\n" + r"across Evolutionary Distance (D) and Sequence Length (L) [Taxa N = 32]", fontsize=15, fontweight="bold", y=1.02)
    plt.tight_layout()

    out_png = os.path.join(outdir, "alpha_scaling_curves.png")
    out_pdf = os.path.join(outdir, "alpha_scaling_curves.pdf")
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.savefig(out_pdf, bbox_inches="tight")
    plt.close()
    print(f"Saved alpha scaling curve plot: {out_png}")
    print(f"Saved alpha scaling curve plot (PDF): {out_pdf}")

def generate_alpha_summary_tables(df, outdir):
    """Generates detailed CSV summary and pivot tables across Alpha values."""
    groupby_cols = ["alpha", "distance", "length", "pipeline"] if "alpha" in df.columns else ["distance", "length", "pipeline"]
    grouped = df.groupby(groupby_cols)["nrf_distance"]
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

    out_csv = os.path.join(outdir, "benchmark_alpha_statistics.csv")
    stats_df.to_csv(out_csv, index=False)
    print(f"Saved alpha statistics summary: {out_csv}")

def main():
    parser = argparse.ArgumentParser(description="Plot Alpha Benchmark results (alpha=0.25, 0.5, 1.0, 2.0)")
    parser.add_argument("--csv", required=True, help="Path to benchmark_alpha_summary.csv")
    parser.add_argument("--outdir", default=".", help="Output directory for plots and summary tables")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    df = pd.read_csv(args.csv)

    print(f"=== Alpha Benchmark Plotting ===")
    if "alpha" in df.columns:
        print(f"Loaded {len(df)} records across alphas: {sorted(df['alpha'].unique())}")
    else:
        print(f"Loaded {len(df)} records (no explicit 'alpha' column)")

    plot_alpha_scaling_lines(df, args.outdir)
    generate_alpha_summary_tables(df, args.outdir)

    print("Alpha plotting completed successfully!")

if __name__ == "__main__":
    main()
