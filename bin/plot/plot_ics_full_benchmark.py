#!/usr/bin/env python3
"""
Full-length Invariant Category Sites (ICS Full) Benchmark Plotting Script
Location: bin/plot/plot_ics_full_benchmark.py

Generates:
1. Multi-panel line plots of NRF vs Evolutionary Distance (D) across Sequence Lengths (L).
2. Comparison between PWA+NJ (PSA+NJ), MSA+NJ, and MSA+ML under 100% full-length ICS with indels.
3. Summary statistics table.
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

PIPELINE_COLORS = {
    "PSA+NJ": "#1f77b4",  # Blue
    "PWA+NJ": "#1f77b4",  # Blue alias
    "MSA+NJ": "#ff7f0e",  # Orange
    "MSA+ML": "#2ca02c",  # Green
}

PIPELINE_ORDER = ["PSA+NJ", "MSA+NJ", "MSA+ML"]

def plot_ics_full_results(df: pd.DataFrame, outdir: str) -> None:
    """Plots NRF Distance vs Evolutionary Distance across sequence lengths."""
    sns.set_theme(style="whitegrid")
    
    # Standardize pipeline name
    df = df.copy()
    df["pipeline"] = df["pipeline"].str.replace("PWA", "PSA")
    
    nrf_col = "nrf_distance" if "nrf_distance" in df.columns else "nrf"
    lengths = sorted(df["length"].unique())

    fig, axes = plt.subplots(
        1, len(lengths),
        figsize=(4.0 * len(lengths), 4.2),
        sharey=True
    )
    if len(lengths) == 1:
        axes = [axes]

    for idx, length in enumerate(lengths):
        ax = axes[idx]
        sub = df[df["length"] == length]

        sns.lineplot(
            data=sub, x="distance", y=nrf_col, hue="pipeline", style="pipeline",
            hue_order=PIPELINE_ORDER, style_order=PIPELINE_ORDER,
            palette=PIPELINE_COLORS, markers=True, dashes=False, ax=ax,
            errorbar="sd", linewidth=1.8
        )
        ax.set_title(f"Length L = {length} aa", fontweight="bold", fontsize=11)
        ax.set_xlabel("Evolutionary Distance (D)", fontsize=10)
        ax.set_ylabel("NRF Distance" if idx == 0 else "", fontsize=10)
        ax.set_ylim(-0.02, 1.02)

        if idx < len(lengths) - 1 and ax.get_legend():
            ax.get_legend().remove()

    if axes[-1].get_legend():
        axes[-1].legend(bbox_to_anchor=(1.05, 1.0), loc="upper left", frameon=True, title="Pipeline")

    plt.subplots_adjust(wspace=0.08, top=0.85, right=0.85)
    plt.suptitle("Topology Accuracy: Full-length ICS Model with Indels (ICS = 100%, D = 0.1 to 3.0)",
                 fontweight="bold", fontsize=13, y=0.98)

    out_png = os.path.join(outdir, "ics_full_benchmark.png")
    out_pdf = os.path.join(outdir, "ics_full_benchmark.pdf")
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.savefig(out_pdf, bbox_inches="tight")
    plt.close()
    print(f"Plot successfully saved to {out_png}")

def export_statistics(df: pd.DataFrame, outdir: str) -> None:
    """Exports summary statistics table."""
    df = df.copy()
    df["pipeline"] = df["pipeline"].str.replace("PWA", "PSA")
    nrf_col = "nrf_distance" if "nrf_distance" in df.columns else "nrf"

    group_cols = ["pipeline", "distance", "length"]
    agg = df.groupby(group_cols)[nrf_col].agg(
        mean="mean",
        std="std",
        median="median",
        count="count"
    ).reset_index()

    out_csv = os.path.join(outdir, "benchmark_ics_full_statistics.csv")
    agg.to_csv(out_csv, index=False)
    print(f"Summary statistics table saved to {out_csv}")

def main() -> None:
    parser = argparse.ArgumentParser(description="Plot Full-length ICS Benchmark Results")
    parser.add_argument("--csv", required=True, help="Input summary CSV file (benchmark_ics_full_summary.csv)")
    parser.add_argument("--outdir", default="results/results_ics_full", help="Output directory for plots")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    df = pd.read_csv(args.csv)

    plot_ics_full_results(df, args.outdir)
    export_statistics(df, args.outdir)

if __name__ == "__main__":
    main()
