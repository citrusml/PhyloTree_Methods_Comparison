#!/usr/bin/env python3
"""
Plot Regime Map Script
Generates 2D phase diagrams (heatmaps) of evolutionary distance (D) vs sequence length (L)
showing relative performance difference ΔRF = RF(PWA+NJ) - RF(MSA+ML) and Regime Boundary (ΔRF = 0).
"""

import sys
import os
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def main():
    parser = argparse.ArgumentParser(description="Generate Regime Map Heatmaps")
    parser.add_argument("--csv", required=True, help="Aggregated results CSV file")
    parser.add_argument("--outdir", default="results", help="Output directory for plots")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    df = pd.read_csv(args.csv)
    if df.empty:
        print("CSV file is empty. Skipping plot generation.")
        return

    # Pivot table for mean nRF by distance, length, and pipeline
    pivot_df = df.groupby(["distance", "length", "pipeline"])["nrf_distance"].mean().reset_index()

    pwa = pivot_df[pivot_df["pipeline"] == "PWA+NJ"].set_index(["length", "distance"])["nrf_distance"]
    msa_ml = pivot_df[pivot_df["pipeline"] == "MSA+ML"].set_index(["length", "distance"])["nrf_distance"]
    msa_nj = pivot_df[pivot_df["pipeline"] == "MSA+NJ"].set_index(["length", "distance"])["nrf_distance"]

    delta_rf = pwa - msa_ml
    delta_rf_grid = delta_rf.unstack(level="distance")

    plt.figure(figsize=(9, 7))
    sns.set_theme(style="whitegrid", font_scale=1.1)

    # Diverging heatmap: Blue = PWA+NJ wins (ΔnRF < 0), Red = MSA+ML wins (ΔnRF > 0), White = Regime Boundary (ΔnRF = 0)
    cmap = sns.diverging_palette(240, 10, as_cmap=True)
    ax = sns.heatmap(
        delta_rf_grid,
        annot=True,
        fmt=".3f",
        cmap=cmap,
        center=0,
        cbar_kws={'label': r'$\Delta nRF = nRF_{\text{PWA+NJ}} - nRF_{\text{MSA+ML}}$'}
    )

    plt.title("Regime Map: PWA+NJ vs MSA+ML Relative Performance", fontsize=14, fontweight="bold", pad=15)
    plt.xlabel("Evolutionary Distance D (substitutions/site)", fontsize=12)
    plt.ylabel("Sequence Length L (amino acids)", fontsize=12)

    plot_path = os.path.join(args.outdir, "regime_map_delta_nrf.png")
    plt.tight_layout()
    plt.savefig(plot_path, dpi=300)
    plt.close()

    print(f"Regime Map heatmap successfully saved to {plot_path}")

if __name__ == "__main__":
    main()
