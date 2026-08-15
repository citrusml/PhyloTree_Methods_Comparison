#!/usr/bin/env python3
"""
Analysis Script: Normalized RF Distance Boxplots with Mean and Variance Annotations
Location: bin/plot_nrf_boxplots.py

Compares PWA+NJ, MSA+NJ, and MSA+ML across all 25 (Distance x Length) conditions.
Annotates each box with mean (μ) and variance (σ²) on labels / plots.
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

PIPELINE_ORDER = ["PWA+NJ", "MSA+NJ", "MSA+ML"]

def find_default_csv():
    candidates = [
        "results/benchmark_summary.csv",
        "results/results_ver2/benchmark_summary.csv",
        "benchmark_summary.csv"
    ]
    for c in candidates:
        if os.path.exists(c):
            return c
    return candidates[0]

def plot_5x5_grid_boxplots(df, outdir):
    distances = sorted(df["distance"].unique())
    lengths = sorted(df["length"].unique())

    fig, axes = plt.subplots(
        nrows=len(distances),
        ncols=len(lengths),
        figsize=(24, 20),
        sharey=True
    )

    for r, d in enumerate(distances):
        for c, l in enumerate(lengths):
            ax = axes[r, c]
            sub_df = df[(df["distance"] == d) & (df["length"] == l)]

            labels = []
            box_data = []
            colors = []

            for pipe in PIPELINE_ORDER:
                pipe_vals = sub_df[sub_df["pipeline"] == pipe]["nrf_distance"].values
                if len(pipe_vals) > 0:
                    mean_val = np.mean(pipe_vals)
                    var_val = np.var(pipe_vals, ddof=1) if len(pipe_vals) > 1 else 0.0
                    box_data.append(pipe_vals)
                    colors.append(PIPELINE_COLORS.get(pipe, "#888888"))
                    labels.append(f"{pipe}\nμ={mean_val:.3f}\nσ²={var_val:.4f}")
                else:
                    box_data.append([])
                    colors.append(PIPELINE_COLORS.get(pipe, "#888888"))
                    labels.append(f"{pipe}\nN/A")

            bp = ax.boxplot(
                box_data,
                tick_labels=labels,
                patch_artist=True,
                showmeans=True,
                meanline=True,
                widths=0.6,
                medianprops=dict(color="black", linewidth=1.5),
                meanprops=dict(color="red", linestyle="--", linewidth=1.5),
                whiskerprops=dict(color="#444444", linewidth=1.2),
                capprops=dict(color="#444444", linewidth=1.2),
                flierprops=dict(marker='o', markersize=3, markerfacecolor='gray', alpha=0.5)
            )

            for patch, col in zip(bp['boxes'], colors):
                patch.set_facecolor(col)
                patch.set_alpha(0.65)
                patch.set_edgecolor("#333333")

            ax.grid(True, axis="y", linestyle="--", alpha=0.6)
            ax.set_ylim(-0.05, 1.05)
            ax.tick_params(axis='x', labelsize=8.5)
            ax.tick_params(axis='y', labelsize=9)

            if r == 0:
                ax.set_title(f"Length L = {l} aa", fontsize=13, fontweight="bold", pad=8)
            if c == 0:
                ax.set_ylabel(f"D = {d}\nNormalized RF", fontsize=11, fontweight="bold")

    legend_elements = [
        plt.Line2D([0], [0], color='black', lw=1.5, label='Median'),
        plt.Line2D([0], [0], color='red', linestyle='--', lw=1.5, label='Mean (μ)'),
        plt.Rectangle((0, 0), 1, 1, facecolor=PIPELINE_COLORS["PWA+NJ"], alpha=0.65, edgecolor='#333', label='PWA+NJ'),
        plt.Rectangle((0, 0), 1, 1, facecolor=PIPELINE_COLORS["MSA+NJ"], alpha=0.65, edgecolor='#333', label='MSA+NJ'),
        plt.Rectangle((0, 0), 1, 1, facecolor=PIPELINE_COLORS["MSA+ML"], alpha=0.65, edgecolor='#333', label='MSA+ML')
    ]
    fig.legend(handles=legend_elements, loc='upper center', bbox_to_anchor=(0.5, 0.995), ncol=5, fontsize=12, frameon=True)

    plt.suptitle("Normalized RF Distance Boxplots across Evolutionary Distance (D) and Sequence Length (L)\nLabels: Mean (μ) and Sample Variance (σ²)", fontsize=16, fontweight="bold", y=1.015)
    plt.tight_layout()

    out_png = os.path.join(outdir, "nrf_boxplots_5x5_grid.png")
    out_pdf = os.path.join(outdir, "nrf_boxplots_5x5_grid.pdf")
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.savefig(out_pdf, bbox_inches="tight")
    plt.close()
    print(f"Saved 5x5 grid boxplot: {out_png}")
    print(f"Saved 5x5 grid boxplot (PDF): {out_pdf}")

def plot_distance_faceted_boxplots(df, outdir):
    distances = sorted(df["distance"].unique())
    lengths = sorted(df["length"].unique())
    sub_outdir = os.path.join(outdir, "by_distance")
    os.makedirs(sub_outdir, exist_ok=True)

    for d in distances:
        fig, axes = plt.subplots(1, len(lengths), figsize=(20, 5.5), sharey=True)
        sub_d = df[df["distance"] == d]

        for i, l in enumerate(lengths):
            ax = axes[i]
            sub_df = sub_d[sub_d["length"] == l]

            labels = []
            box_data = []
            colors = []

            for pipe in PIPELINE_ORDER:
                pipe_vals = sub_df[sub_df["pipeline"] == pipe]["nrf_distance"].values
                if len(pipe_vals) > 0:
                    mean_val = np.mean(pipe_vals)
                    var_val = np.var(pipe_vals, ddof=1) if len(pipe_vals) > 1 else 0.0
                    box_data.append(pipe_vals)
                    colors.append(PIPELINE_COLORS.get(pipe, "#888888"))
                    labels.append(f"{pipe}\nμ={mean_val:.3f}\nσ²={var_val:.4f}")
                else:
                    box_data.append([])
                    colors.append(PIPELINE_COLORS.get(pipe, "#888888"))
                    labels.append(f"{pipe}\nN/A")

            bp = ax.boxplot(
                box_data,
                tick_labels=labels,
                patch_artist=True,
                showmeans=True,
                meanline=True,
                widths=0.6,
                medianprops=dict(color="black", linewidth=1.5),
                meanprops=dict(color="red", linestyle="--", linewidth=1.5),
                whiskerprops=dict(color="#444444", linewidth=1.2),
                capprops=dict(color="#444444", linewidth=1.2),
                flierprops=dict(marker='o', markersize=3, markerfacecolor='gray', alpha=0.5)
            )

            for patch, col in zip(bp['boxes'], colors):
                patch.set_facecolor(col)
                patch.set_alpha(0.65)
                patch.set_edgecolor("#333333")

            ax.set_title(f"L = {l} aa", fontsize=12, fontweight="bold")
            ax.grid(True, axis="y", linestyle="--", alpha=0.6)
            ax.set_ylim(-0.05, 1.05)
            ax.tick_params(axis='x', labelsize=9)
            if i == 0:
                ax.set_ylabel("Normalized RF (nRF)", fontsize=11, fontweight="bold")

        fig.suptitle(f"Normalized RF Boxplot Comparison at Distance D = {d} substitutions/site", fontsize=14, fontweight="bold", y=1.02)
        plt.tight_layout()

        out_path = os.path.join(sub_outdir, f"nrf_boxplots_D{d}.png")
        plt.savefig(out_path, dpi=300, bbox_inches="tight")
        plt.close()

    print(f"Saved distance-faceted boxplots to {sub_outdir}/")

def generate_summary_stats_table(df, outdir):
    stats_list = []
    grouped = df.groupby(["distance", "length", "pipeline"])

    for (d, l, pipe), group in grouped:
        nrf = group["nrf_distance"].values
        stats_list.append({
            "distance": d,
            "length": l,
            "pipeline": pipe,
            "replicates": len(nrf),
            "nrf_mean": round(float(np.mean(nrf)), 5),
            "nrf_variance": round(float(np.var(nrf, ddof=1)), 6) if len(nrf) > 1 else 0.0,
            "nrf_std": round(float(np.std(nrf, ddof=1)), 5) if len(nrf) > 1 else 0.0,
            "nrf_median": round(float(np.median(nrf)), 5),
            "nrf_q25": round(float(np.percentile(nrf, 25)), 5),
            "nrf_q75": round(float(np.percentile(nrf, 75)), 5),
            "nrf_min": round(float(np.min(nrf)), 5),
            "nrf_max": round(float(np.max(nrf)), 5),
        })

    stats_df = pd.DataFrame(stats_list)
    stats_df = stats_df.sort_values(by=["distance", "length", "pipeline"])

    csv_path = os.path.join(outdir, "nrf_summary_statistics.csv")
    stats_df.to_csv(csv_path, index=False)
    print(f"Saved summary statistics table: {csv_path}")

    pivot_mean = stats_df.pivot(index=["distance", "length"], columns="pipeline", values="nrf_mean")
    pivot_var  = stats_df.pivot(index=["distance", "length"], columns="pipeline", values="nrf_variance")
    
    pivot_combined = pd.DataFrame()
    for p in PIPELINE_ORDER:
        if p in pivot_mean.columns and p in pivot_var.columns:
            pivot_combined[f"{p} (mean ± var)"] = pivot_mean[p].map(lambda m: f"{m:.3f}") + " ± " + pivot_var[p].map(lambda v: f"{v:.4f}")
    
    pivot_path = os.path.join(outdir, "nrf_comparison_pivot.csv")
    pivot_combined.to_csv(pivot_path)
    print(f"Saved comparison pivot table: {pivot_path}")

def main():
    parser = argparse.ArgumentParser(description="Generate nRF Boxplots with Mean and Variance for all conditions")
    parser.add_argument("--csv", default=find_default_csv(), help="Path to benchmark_summary.csv")
    parser.add_argument("--outdir", default="results/figures", help="Output directory for plots and tables (default: results/figures)")
    args = parser.parse_args()

    print(f"=== Normalized RF Boxplot Analysis ===")
    print(f"Input CSV       : {args.csv}")
    print(f"Output directory: {args.outdir}")

    if not os.path.exists(args.csv):
        print(f"Error: CSV file not found: {args.csv}")
        sys.exit(1)

    os.makedirs(args.outdir, exist_ok=True)
    df = pd.read_csv(args.csv)

    print(f"Loaded {len(df)} evaluation records.")

    plot_5x5_grid_boxplots(df, args.outdir)
    plot_distance_faceted_boxplots(df, args.outdir)
    generate_summary_stats_table(df, args.outdir)

    print(f"\nAll figures and tables saved to: {os.path.abspath(args.outdir)}")

if __name__ == "__main__":
    main()
