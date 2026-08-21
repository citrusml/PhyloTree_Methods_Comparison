#!/usr/bin/env python3
"""
Method Comparison Regime Maps Script:
Generates square-shaped 2D phase diagrams (heatmaps) exactly like `regime_map_delta_nrf.png` for:
1. MSA+ML vs MSA+NJ  (Reconstruction Engine: ML vs NJ)
2. MSA+NJ vs PWA+NJ  (Alignment Strategy: MSA vs PWA)
3. PWA+NJ vs MSA+ML  (Overall Pipeline Comparison)

Features:
- Exact square-cell aspect ratio matching regime_map_delta_nrf.png
- Diverging color palettes with zero-centering
- Clean 3-decimal annotation (.3f) or optional significance stars
- Output in both individual square PNGs and side-by-side comparison
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def compute_delta_grids(df, metric="nrf_distance"):
    """
    Computes mean delta grids across distance and length.
    """
    # Pivot table for mean metric by distance, length, and pipeline
    pivot_df = df.groupby(["distance", "length", "pipeline"])[metric].mean().reset_index()

    msa_ml = pivot_df[pivot_df["pipeline"] == "MSA+ML"].set_index(["length", "distance"])[metric]
    msa_nj = pivot_df[pivot_df["pipeline"] == "MSA+NJ"].set_index(["length", "distance"])[metric]
    pwa_nj = pivot_df[pivot_df["pipeline"] == "PWA+NJ"].set_index(["length", "distance"])[metric]

    # Delta 1: MSA+ML vs MSA+NJ -> (MSA+NJ - MSA+ML) : > 0 means MSA+ML wins (red), < 0 means MSA+NJ wins (blue)
    delta_ml_vs_nj = (msa_nj - msa_ml).unstack(level="distance")

    # Delta 2: MSA+NJ vs PWA+NJ -> (PWA+NJ - MSA+NJ) : > 0 means MSA+NJ wins (red), < 0 means PWA+NJ wins (blue)
    delta_nj_vs_pwa = (pwa_nj - msa_nj).unstack(level="distance")

    # Delta 3: PWA+NJ vs MSA+ML -> (PWA+NJ - MSA+ML) : > 0 means MSA+ML wins (red), < 0 means PWA+NJ wins (blue)
    delta_pwa_vs_ml = (pwa_nj - msa_ml).unstack(level="distance")

    # Absolute grids
    grid_ml = msa_ml.unstack(level="distance")
    grid_nj = msa_nj.unstack(level="distance")
    grid_pwa = pwa_nj.unstack(level="distance")

    return {
        "ml_vs_nj": delta_ml_vs_nj,
        "nj_vs_pwa": delta_nj_vs_pwa,
        "pwa_vs_ml": delta_pwa_vs_ml,
        "abs_ml": grid_ml,
        "abs_nj": grid_nj,
        "abs_pwa": grid_pwa
    }

def plot_square_regime_map(grid_df, title, cbar_label, out_path, vlim=None):
    """
    Plots a square heatmap matching regime_map_delta_nrf.png styling.
    """
    plt.figure(figsize=(9, 7))
    sns.set_theme(style="whitegrid", font_scale=1.1)

    cmap = sns.diverging_palette(240, 10, as_cmap=True)

    if vlim is None:
        max_abs = max(abs(grid_df.min().min()), abs(grid_df.max().max()))
        vlim = max(0.06, np.ceil(max_abs * 20) / 20)

    ax = sns.heatmap(
        grid_df,
        annot=True,
        fmt=".3f",
        cmap=cmap,
        center=0,
        square=True,
        cbar_kws={'label': cbar_label}
    )

    plt.title(title, fontsize=14, fontweight="bold", pad=15)
    plt.xlabel("Evolutionary Distance D (substitutions/site)", fontsize=12)
    plt.ylabel("Sequence Length L (amino acids)", fontsize=12)

    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    pdf_path = out_path.replace(".png", ".pdf")
    plt.savefig(pdf_path, dpi=300)
    plt.close()
    print(f"Saved square heatmap: {out_path} and {pdf_path}")

def plot_side_by_side_square_maps(grids, out_path):
    """
    Plots the two key regime maps side-by-side with square cells.
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 8))
    sns.set_theme(style="whitegrid", font_scale=1.1)

    cmap = sns.diverging_palette(240, 10, as_cmap=True)

    # 1. MSA+ML vs MSA+NJ
    sns.heatmap(
        grids["ml_vs_nj"],
        annot=True,
        fmt=".3f",
        cmap=cmap,
        center=0,
        square=True,
        ax=ax1,
        cbar_kws={'label': r'$\Delta nRF = nRF_{\text{MSA+NJ}} - nRF_{\text{MSA+ML}}$', 'shrink': 0.8}
    )
    ax1.set_title("Regime Map: MSA+ML vs MSA+NJ\n(Reconstruction Engine Comparison)", fontsize=13, fontweight="bold", pad=12)
    ax1.set_xlabel("Evolutionary Distance D (substitutions/site)", fontsize=11)
    ax1.set_ylabel("Sequence Length L (amino acids)", fontsize=11)

    # 2. MSA+NJ vs PWA+NJ
    sns.heatmap(
        grids["nj_vs_pwa"],
        annot=True,
        fmt=".3f",
        cmap=cmap,
        center=0,
        square=True,
        ax=ax2,
        cbar_kws={'label': r'$\Delta nRF = nRF_{\text{PWA+NJ}} - nRF_{\text{MSA+NJ}}$', 'shrink': 0.8}
    )
    ax2.set_title("Regime Map: MSA+NJ vs PWA+NJ\n(Alignment Strategy Comparison)", fontsize=13, fontweight="bold", pad=12)
    ax2.set_xlabel("Evolutionary Distance D (substitutions/site)", fontsize=11)
    ax2.set_ylabel("Sequence Length L (amino acids)", fontsize=11)

    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    plt.savefig(out_path.replace(".png", ".pdf"), dpi=300)
    plt.close()
    print(f"Saved side-by-side square maps: {out_path}")

def main():
    parser = argparse.ArgumentParser(
        description="Generate Square-shaped Regime Map Heatmaps (MSA+ML vs MSA+NJ and MSA+NJ vs PWA+NJ)"
    )
    parser.add_argument("--csv", required=True, help="Path to benchmark_summary.csv")
    parser.add_argument("--outdir", default=None, help="Output directory for heatmap figures")
    parser.add_argument("--taxa", type=int, default=None, help="Filter by specific num_taxa (optional)")
    parser.add_argument("--metric", default="nrf_distance", choices=["nrf_distance", "rf_distance"],
                        help="Distance metric to evaluate (default: nrf_distance)")
    args = parser.parse_args()

    if not os.path.exists(args.csv):
        print(f"Error: CSV file not found: {args.csv}", file=sys.stderr)
        sys.exit(1)

    if args.outdir is None:
        args.outdir = os.path.dirname(args.csv)
    os.makedirs(args.outdir, exist_ok=True)

    print(f"Loading data from: {args.csv}")
    df = pd.read_csv(args.csv)

    if args.taxa is not None and "num_taxa" in df.columns:
        df = df[df["num_taxa"] == args.taxa]
        print(f"Filtered for num_taxa = {args.taxa} (rows: {len(df)})")

    grids = compute_delta_grids(df, metric=args.metric)

    # 1. Square Heatmap: MSA+ML vs MSA+NJ
    plot_square_regime_map(
        grid_df=grids["ml_vs_nj"],
        title="Regime Map: MSA+ML vs MSA+NJ Relative Performance",
        cbar_label=r'$\Delta nRF = nRF_{\text{MSA+NJ}} - nRF_{\text{MSA+ML}}$',
        out_path=os.path.join(args.outdir, "regime_map_msa_ml_vs_msa_nj.png")
    )

    # 2. Square Heatmap: MSA+NJ vs PWA+NJ
    plot_square_regime_map(
        grid_df=grids["nj_vs_pwa"],
        title="Regime Map: MSA+NJ vs PWA+NJ Relative Performance",
        cbar_label=r'$\Delta nRF = nRF_{\text{PWA+NJ}} - nRF_{\text{MSA+NJ}}$',
        out_path=os.path.join(args.outdir, "regime_map_msa_nj_vs_pwa_nj.png")
    )

    # 3. Square Heatmap: PWA+NJ vs MSA+ML (Original Regime Map)
    plot_square_regime_map(
        grid_df=grids["pwa_vs_ml"],
        title="Regime Map: PWA+NJ vs MSA+ML Relative Performance",
        cbar_label=r'$\Delta nRF = nRF_{\text{PWA+NJ}} - nRF_{\text{MSA+ML}}$',
        out_path=os.path.join(args.outdir, "regime_map_pwa_nj_vs_msa_ml.png")
    )

    # 4. Side-by-side comparison
    plot_side_by_side_square_maps(
        grids=grids,
        out_path=os.path.join(args.outdir, "regime_map_comparisons_side_by_side.png")
    )

    print("\nAll square regime map heatmaps generated successfully!")

if __name__ == "__main__":
    main()
