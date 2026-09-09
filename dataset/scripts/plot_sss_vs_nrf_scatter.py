#!/usr/bin/env python3
"""
plot_sss_vs_nrf_scatter.py
==========================
Generates scatterplots of SSS (X-axis) vs nRF (Y-axis) colored by method (hue),
separated by biological group (amn, euk, acti, bac):
  1. sss_vs_nrf_by_group_ml_ref.png (2x2 grid: MSA+NJ, PSA+NJ, GS vs MSA+ML)
  2. sss_vs_nrf_by_group_all_pairs.png (2x2 grid: all 6 method pairs)
  3. Individual scatterplots for each biological group (amn, euk, acti, bac)
"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

SCRIPT_DIR = Path(__file__).resolve().parent
RESULTS_DIR = SCRIPT_DIR.parent / "results"

plt.rcParams.update({
    "font.sans-serif": ["Arial", "DejaVu Sans", "Helvetica", "sans-serif"],
    "font.size": 11,
    "axes.labelsize": 12,
    "axes.titlesize": 13,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 10,
    "figure.titlesize": 15,
    "lines.linewidth": 1.8,
    "axes.grid": True,
    "grid.alpha": 0.3,
    "grid.linestyle": "--"
})

GROUPS = [
    ("amn", "Amniota (10 taxa)"),
    ("euk", "Eukaryota (50 taxa)"),
    ("act", "Actinobacteria (50 taxa)"),
    ("bac", "Bacteria (50 taxa)")
]


def load_merged_data() -> pd.DataFrame:
    """Loads and merges concordance and SSS data."""
    conc_file = RESULTS_DIR / "method_concordance_summary.csv"
    sim_file = RESULTS_DIR / "hog_similarity_summary.csv"

    if not conc_file.exists() or not sim_file.exists():
        raise FileNotFoundError(f"Missing required CSV files: {conc_file} or {sim_file}")

    df_conc = pd.read_csv(conc_file)
    df_sim = pd.read_csv(sim_file)

    cols_sim = ["dataset", "hog_id", "mean_sss", "median_sss", "mean_identity"]
    df_merged = pd.merge(df_conc, df_sim[cols_sim], on=["dataset", "hog_id"], how="inner")
    return df_merged


def plot_ml_reference_scatter(df: pd.DataFrame, out_path: Path):
    """
    2x2 multi-panel scatter plot where methods (MSA+NJ, PSA+NJ, GS)
    are evaluated against MSA+ML (the standard ML benchmark).
    """
    method_comparisons = [
        ("nrf_msa_ml_vs_msa_nj", "MSA+NJ (vs MSA+ML)", "#e67e22", "o"),
        ("nrf_msa_ml_vs_psa_nj", "PSA+NJ (vs MSA+ML)", "#27ae60", "s"),
        ("nrf_msa_ml_vs_gs",     "GS (vs MSA+ML)",     "#c0392b", "^")
    ]

    fig, axes = plt.subplots(2, 2, figsize=(15, 11), sharey=True)
    axes = axes.flatten()

    for idx, (d_key, title_name) in enumerate(GROUPS):
        ax = axes[idx]
        sub = df[df["dataset"] == d_key]

        for col, label, color, marker in method_comparisons:
            valid = sub.dropna(subset=["mean_sss", col])
            ax.scatter(
                valid["mean_sss"],
                valid[col],
                label=label,
                color=color,
                marker=marker,
                alpha=0.65,
                s=45,
                edgecolors="white",
                linewidth=0.5
            )

            # Trend line
            if len(valid) > 2:
                sns.regplot(
                    data=valid,
                    x="mean_sss",
                    y=col,
                    scatter=False,
                    ax=ax,
                    color=color,
                    line_kws={"linewidth": 2, "alpha": 0.85},
                    ci=None
                )

        ax.set_title(title_name, fontweight="bold", fontsize=13)
        ax.set_xlabel("Sequence Similarity Score (SSS: $\\bar{w}$)")
        ax.set_ylabel("Normalized RF (nRF vs MSA+ML)" if idx % 2 == 0 else "")
        ax.set_ylim(-0.03, 1.03)
        ax.set_xlim(-0.02, 1.02)
        ax.legend(frameon=True, loc="upper right", framealpha=0.9)

    plt.suptitle("SSS vs. nRF (Tree Distance from MSA+ML) Across Biological Groups", fontweight="bold", y=0.98)
    plt.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out_path}")


def plot_all_pairs_scatter(df: pd.DataFrame, out_path: Path):
    """
    2x2 multi-panel scatter plot with all 6 pairwise method comparisons colored by method pair.
    """
    all_pairs = [
        ("nrf_msa_ml_vs_msa_nj", "MSA+ML vs MSA+NJ", "#e67e22", "o"),
        ("nrf_msa_ml_vs_psa_nj", "MSA+ML vs PSA+NJ", "#27ae60", "s"),
        ("nrf_msa_ml_vs_gs",     "MSA+ML vs GS",     "#c0392b", "^"),
        ("nrf_msa_nj_vs_psa_nj", "MSA+NJ vs PSA+NJ", "#2980b9", "D"),
        ("nrf_msa_nj_vs_gs",     "MSA+NJ vs GS",     "#8e44ad", "v"),
        ("nrf_psa_nj_vs_gs",     "PSA+NJ vs GS",     "#16a085", "P")
    ]

    fig, axes = plt.subplots(2, 2, figsize=(16, 12), sharey=True)
    axes = axes.flatten()

    for idx, (d_key, title_name) in enumerate(GROUPS):
        ax = axes[idx]
        sub = df[df["dataset"] == d_key]

        for col, label, color, marker in all_pairs:
            valid = sub.dropna(subset=["mean_sss", col])
            ax.scatter(
                valid["mean_sss"],
                valid[col],
                label=label,
                color=color,
                marker=marker,
                alpha=0.6,
                s=35,
                edgecolors="none"
            )

            # Trend line
            if len(valid) > 2:
                sns.regplot(
                    data=valid,
                    x="mean_sss",
                    y=col,
                    scatter=False,
                    ax=ax,
                    color=color,
                    line_kws={"linewidth": 1.5, "alpha": 0.8},
                    ci=None
                )

        ax.set_title(title_name, fontweight="bold", fontsize=13)
        ax.set_xlabel("Sequence Similarity Score (SSS: $\\bar{w}$)")
        ax.set_ylabel("Normalized RF (nRF)" if idx % 2 == 0 else "")
        ax.set_ylim(-0.03, 1.03)
        ax.set_xlim(-0.02, 1.02)
        ax.legend(frameon=True, loc="upper right", framealpha=0.85, fontsize=8.5)

    plt.suptitle("SSS vs. Topological Distance (nRF) for All Method Pairs", fontweight="bold", y=0.98)
    plt.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out_path}")


def plot_individual_group_scatters(df: pd.DataFrame):
    """Generates a separate dedicated high-resolution plot for each biological group."""
    method_comparisons = [
        ("nrf_msa_ml_vs_msa_nj", "MSA+NJ (vs MSA+ML)", "#e67e22", "o"),
        ("nrf_msa_ml_vs_psa_nj", "PSA+NJ (vs MSA+ML)", "#27ae60", "s"),
        ("nrf_msa_ml_vs_gs",     "GS (vs MSA+ML)",     "#c0392b", "^"),
        ("nrf_msa_nj_vs_psa_nj", "MSA+NJ vs PSA+NJ",   "#2980b9", "D"),
        ("nrf_psa_nj_vs_gs",     "PSA+NJ vs GS",       "#16a085", "P")
    ]

    for d_key, title_name in GROUPS:
        fig, ax = plt.subplots(figsize=(8.5, 6))
        sub = df[df["dataset"] == d_key]

        for col, label, color, marker in method_comparisons:
            valid = sub.dropna(subset=["mean_sss", col])
            ax.scatter(
                valid["mean_sss"],
                valid[col],
                label=label,
                color=color,
                marker=marker,
                alpha=0.7,
                s=50,
                edgecolors="black",
                linewidth=0.4
            )
            if len(valid) > 2:
                sns.regplot(
                    data=valid,
                    x="mean_sss",
                    y=col,
                    scatter=False,
                    ax=ax,
                    color=color,
                    line_kws={"linewidth": 2.2, "alpha": 0.85},
                    ci=None
                )

        ax.set_title(f"SSS vs. nRF: {title_name} (100 HOGs)", fontweight="bold", fontsize=13)
        ax.set_xlabel("Sequence Similarity Score (SSS: $\\bar{w}$)", fontsize=12)
        ax.set_ylabel("Normalized Robinson-Foulds Distance (nRF)", fontsize=12)
        ax.set_ylim(-0.03, 1.03)
        ax.set_xlim(-0.02, 1.02)
        ax.legend(frameon=True, loc="upper right", framealpha=0.9)
        plt.tight_layout()

        out_name = RESULTS_DIR / f"scatter_sss_vs_nrf_{d_key}.png"
        fig.savefig(out_name, dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(f"Saved: {out_name}")


def main():
    df = load_merged_data()
    print(f"Loaded {len(df)} merged records.")

    out_ml_ref = RESULTS_DIR / "scatter_sss_vs_nrf_by_group_ml_ref.png"
    out_all_pairs = RESULTS_DIR / "scatter_sss_vs_nrf_by_group_all_pairs.png"

    plot_ml_reference_scatter(df, out_ml_ref)
    plot_all_pairs_scatter(df, out_all_pairs)
    plot_individual_group_scatters(df)
    print("All scatterplots generated successfully!")


if __name__ == "__main__":
    main()
