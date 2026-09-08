#!/usr/bin/env python3
"""
plot_dataset_similarity.py
==========================
Generates publication-quality figures comparing Sequence Similarity Scores (SSS)
and Pairwise Identity across datasets in dataset/ against the Superfamily thresholds
from Matsui & Iwasaki (Syst. Biol., 2020).

Outputs saved in dataset/results/:
- dataset_similarity_multiplot.png
- sss_distribution_comparison.png
- identity_distribution_comparison.png
- sss_vs_identity_correlation.png
"""

import os
import sys
from pathlib import Path
from typing import Dict, Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent.parent
DATASET_DIR = SCRIPT_DIR.parent
DEFAULT_RESULTS_DIR = DATASET_DIR / "results"


def setup_matplotlib_style():
    """Sets publication-level styling for Matplotlib and Seaborn."""
    plt.style.use("seaborn-v0_8-whitegrid")
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.size": 12,
        "axes.labelsize": 14,
        "axes.titlesize": 15,
        "axes.titleweight": "bold",
        "xtick.labelsize": 12,
        "ytick.labelsize": 12,
        "legend.fontsize": 11,
        "legend.title_fontsize": 12,
        "figure.titlesize": 16,
        "figure.dpi": 300,
        "savefig.dpi": 300,
        "savefig.bbox": "tight"
    })


def plot_sss_distribution(df_hogs: pd.DataFrame, out_path: Path):
    """Plots SSS distribution violin/strip plot with Superfamily baseline markers."""
    fig, ax = plt.subplots(figsize=(9, 6))

    palette = {"amn": "#2b5c8f", "euk": "#2ca02c", "bac": "#d95f02", "act": "#7570b3"}
    labels = {
        "amn": "Amniota\n(10 taxa)",
        "euk": "Eukaryota\n(50 taxa)",
        "bac": "Bacteria\n(50 taxa)",
        "act": "Actinobacteria\n(50 taxa)"
    }

    df = df_hogs.copy()
    df["dataset_label"] = df["dataset"].map(lambda k: labels.get(k, k))

    sns.violinplot(
        data=df,
        x="dataset_label",
        y="mean_sss",
        hue="dataset_label",
        legend=False,
        palette=[palette.get(k, "#333333") for k in df["dataset"].unique()],
        inner="quartile",
        cut=0,
        ax=ax,
        alpha=0.6
    )
    sns.stripplot(
        data=df,
        x="dataset_label",
        y="mean_sss",
        color="black",
        alpha=0.4,
        size=4,
        jitter=0.2,
        ax=ax
    )

    # Reference lines for Superfamily / Twilight zone
    ax.axhline(0.10, color="#d62728", linestyle="--", linewidth=2.0, label="Superfamily Upper Bound (SSS=0.10, GS/NJ > ML)")
    ax.axhline(0.07, color="#ff7f0e", linestyle=":", linewidth=2.0, label="SCOPe 243 Superfamilies Median (SSS=0.07)")

    # Fill twilight zone region
    ax.axhspan(0.0, 0.10, color="#d62728", alpha=0.12, label="Superfamily / Twilight Zone (SSS ≤ 0.10)")

    ax.set_ylim(-0.02, 1.02)
    ax.set_ylabel("Mean Sequence Similarity Score (SSS, $w$)")
    ax.set_xlabel("Dataset / Taxonomic Group")
    ax.set_title("Distribution of Sequence Similarity Scores (SSS) across HOGs")
    ax.legend(loc="upper right", frameon=True, framealpha=0.9)

    plt.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)
    print(f"Saved SSS distribution plot to: {out_path}")


def plot_identity_distribution(df_hogs: pd.DataFrame, out_path: Path):
    """Plots Pairwise Identity (%) distribution violin plot with Twilight Zone baseline."""
    fig, ax = plt.subplots(figsize=(9, 6))

    palette = {"amn": "#2b5c8f", "euk": "#2ca02c", "bac": "#d95f02", "act": "#7570b3"}
    labels = {
        "amn": "Amniota\n(10 taxa)",
        "euk": "Eukaryota\n(50 taxa)",
        "bac": "Bacteria\n(50 taxa)",
        "act": "Actinobacteria\n(50 taxa)"
    }

    df = df_hogs.copy()
    df["dataset_label"] = df["dataset"].map(lambda k: labels.get(k, k))
    df["mean_identity_pct"] = df["mean_identity"] * 100.0

    sns.violinplot(
        data=df,
        x="dataset_label",
        y="mean_identity_pct",
        hue="dataset_label",
        legend=False,
        palette=[palette.get(k, "#333333") for k in df["dataset"].unique()],
        inner="quartile",
        cut=0,
        ax=ax,
        alpha=0.6
    )
    sns.stripplot(
        data=df,
        x="dataset_label",
        y="mean_identity_pct",
        color="black",
        alpha=0.4,
        size=4,
        jitter=0.2,
        ax=ax
    )

    # Reference lines for Twilight zone (10-20%)
    ax.axhline(20.0, color="#d62728", linestyle="--", linewidth=2.0, label="Twilight Zone Boundary (20% Identity)")
    ax.axhline(10.0, color="#9467bd", linestyle=":", linewidth=2.0, label="Midnight Zone Boundary (10% Identity)")
    ax.axhspan(0.0, 20.0, color="#d62728", alpha=0.12, label="Superfamily Twilight Zone (10% - 20%)")

    ax.set_ylim(-2.0, 102.0)
    ax.set_ylabel("Mean Pairwise Sequence Identity (%)")
    ax.set_xlabel("Dataset / Taxonomic Group")
    ax.set_title("Distribution of Sequence Identity across HOGs")
    ax.legend(loc="lower right", frameon=True, framealpha=0.9)

    plt.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)
    print(f"Saved Identity distribution plot to: {out_path}")


def plot_sss_vs_identity_correlation(df_hogs: pd.DataFrame, out_path: Path):
    """Scatter plot of SSS vs Identity across all individual HOGs."""
    fig, ax = plt.subplots(figsize=(8, 6))

    palette = {"amn": "#2b5c8f", "euk": "#2ca02c", "bac": "#d95f02", "act": "#7570b3"}
    labels = {
        "amn": "Amniota (10 taxa)",
        "euk": "Eukaryota (50 taxa)",
        "bac": "Bacteria (50 taxa)",
        "act": "Actinobacteria (50 taxa)"
    }

    df = df_hogs.copy()
    df["mean_identity_pct"] = df["mean_identity"] * 100.0

    for d_key, grp in df.groupby("dataset"):
        ax.scatter(
            grp["mean_identity_pct"],
            grp["mean_sss"],
            label=labels.get(d_key, d_key),
            color=palette.get(d_key, "#333333"),
            alpha=0.75,
            edgecolors="w",
            s=55
        )

    # Shaded superfamily region
    ax.axhspan(0.0, 0.10, color="#d62728", alpha=0.12, label="Superfamily Regime (SSS ≤ 0.10)")
    ax.axvline(20.0, color="#d62728", linestyle="--", linewidth=1.5, label="Twilight Zone (20% Id)")

    ax.set_xlim(0, 100)
    ax.set_ylim(0, 1.0)
    ax.set_xlabel("Mean Pairwise Sequence Identity (%)")
    ax.set_ylabel("Mean Sequence Similarity Score (SSS, $w$)")
    ax.set_title("Correlation between SSS and Pairwise Identity in HOGs")
    ax.legend(loc="upper left", frameon=True, framealpha=0.9)

    plt.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)
    print(f"Saved SSS vs Identity scatter plot to: {out_path}")


def plot_multiplot_summary(df_hogs: pd.DataFrame, out_path: Path):
    """Creates a 2-panel composite figure for overview."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))

    palette = {"amn": "#2b5c8f", "euk": "#2ca02c", "bac": "#d95f02", "act": "#7570b3"}
    labels = {"amn": "Amniota", "euk": "Eukaryota", "bac": "Bacteria", "act": "Actinobacteria"}

    df = df_hogs.copy()
    df["dataset_label"] = df["dataset"].map(lambda k: labels.get(k, k))
    df["mean_identity_pct"] = df["mean_identity"] * 100.0

    # Panel A: SSS
    sns.violinplot(
        data=df, x="dataset_label", y="mean_sss",
        hue="dataset_label", legend=False,
        palette=[palette.get(k, "#333333") for k in df["dataset"].unique()],
        inner="quartile", cut=0, ax=ax1, alpha=0.6
    )
    sns.stripplot(data=df, x="dataset_label", y="mean_sss", color="black", alpha=0.35, size=3.5, jitter=0.2, ax=ax1)
    ax1.axhline(0.10, color="#d62728", linestyle="--", linewidth=2.0, label="Superfamily Upper Bound (SSS ≤ 0.10)")
    ax1.axhspan(0.0, 0.10, color="#d62728", alpha=0.12)
    ax1.set_ylim(-0.02, 1.02)
    ax1.set_ylabel("Mean SSS ($w$)")
    ax1.set_xlabel("Dataset")
    ax1.set_title("(A) Sequence Similarity Score (SSS) Distribution")
    ax1.legend(loc="upper right", frameon=True)

    # Panel B: Identity
    sns.violinplot(
        data=df, x="dataset_label", y="mean_identity_pct",
        hue="dataset_label", legend=False,
        palette=[palette.get(k, "#333333") for k in df["dataset"].unique()],
        inner="quartile", cut=0, ax=ax2, alpha=0.6
    )
    sns.stripplot(data=df, x="dataset_label", y="mean_identity_pct", color="black", alpha=0.35, size=3.5, jitter=0.2, ax=ax2)
    ax2.axhline(20.0, color="#d62728", linestyle="--", linewidth=2.0, label="Twilight Zone (20% Identity)")
    ax2.axhspan(0.0, 20.0, color="#d62728", alpha=0.12)
    ax2.set_ylim(-2.0, 102.0)
    ax2.set_ylabel("Mean Identity (%)")
    ax2.set_xlabel("Dataset")
    ax2.set_title("(B) Pairwise Sequence Identity (%) Distribution")
    ax2.legend(loc="lower right", frameon=True)

    plt.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)
    print(f"Saved multiplot summary to: {out_path}")


def main():
    setup_matplotlib_style()
    csv_path = DEFAULT_RESULTS_DIR / "hog_similarity_summary.csv"
    if not csv_path.exists():
        print(f"Error: Summary CSV not found at {csv_path}. Run analyze_dataset_similarity.py first.")
        sys.exit(1)

    df_hogs = pd.read_csv(csv_path)
    print(f"Loaded {len(df_hogs)} HOG records from {csv_path}.")

    plot_sss_distribution(df_hogs, DEFAULT_RESULTS_DIR / "sss_distribution_comparison.png")
    plot_identity_distribution(df_hogs, DEFAULT_RESULTS_DIR / "identity_distribution_comparison.png")
    plot_sss_vs_identity_correlation(df_hogs, DEFAULT_RESULTS_DIR / "sss_vs_identity_correlation.png")
    plot_multiplot_summary(df_hogs, DEFAULT_RESULTS_DIR / "dataset_similarity_multiplot.png")
    print("All plots generated successfully.")


if __name__ == "__main__":
    main()
