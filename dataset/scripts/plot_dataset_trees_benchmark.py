#!/usr/bin/env python3
"""
plot_dataset_trees_benchmark.py
===============================
Generates publication-quality visualizations (300 DPI) for the real dataset
phylogenetic tree benchmark across Amniota, Eukaryota, Bacteria, and Actinobacteria:
  1. Method Concordance Violin/Box Plots (pairwise nRF across biological groups)
  2. Average nRF Heatmap among the 4 methods (MSA+ML, MSA+NJ, PSA+NJ, GS)
  3. SSS vs. Method Disagreement (nRF) Correlation Plots
"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

SCRIPT_DIR = Path(__file__).resolve().parent
RESULTS_DIR = SCRIPT_DIR.parent / "results"

# Set publication style
plt.rcParams.update({
    "font.sans-serif": ["Arial", "DejaVu Sans", "Helvetica", "sans-serif"],
    "font.size": 11,
    "axes.labelsize": 12,
    "axes.titlesize": 13,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 10,
    "figure.titlesize": 14,
    "lines.linewidth": 1.5,
    "axes.grid": True,
    "grid.alpha": 0.3,
    "grid.linestyle": "--"
})

GROUP_COLORS = {
    "Amniota": "#2b5c8f",
    "Eukaryota": "#3282b8",
    "Actinobacteria": "#e07a5f",
    "Bacteria": "#d9534f"
}

METHOD_LABELS = {
    "msa_ml": "MSA+ML (IQ-TREE 2)",
    "msa_nj": "MSA+NJ (RapidNJ)",
    "psa_nj": "PSA+NJ (RapidNJ)",
    "gs": "GS (Graph Splitting)"
}


def plot_method_concordance_violins(df_concordance: pd.DataFrame, out_path: Path):
    """Plots violin/box plots of normalized Robinson-Foulds (nRF) distance across groups."""
    pairs = [
        ("nrf_msa_ml_vs_msa_nj", "MSA+ML vs MSA+NJ"),
        ("nrf_msa_ml_vs_psa_nj", "MSA+ML vs PSA+NJ"),
        ("nrf_msa_ml_vs_gs", "MSA+ML vs GS"),
        ("nrf_msa_nj_vs_psa_nj", "MSA+NJ vs PSA+NJ"),
        ("nrf_msa_nj_vs_gs", "MSA+NJ vs GS"),
        ("nrf_psa_nj_vs_gs", "PSA+NJ vs GS")
    ]

    records = []
    for _, row in df_concordance.iterrows():
        g_name = row["dataset_name"]
        for col, pair_label in pairs:
            val = row.get(col)
            if pd.notna(val):
                records.append({
                    "Group": g_name,
                    "Comparison": pair_label,
                    "nRF": val
                })

    df_long = pd.DataFrame(records)
    if df_long.empty:
        return

    order_groups = ["Amniota", "Eukaryota", "Actinobacteria", "Bacteria"]
    order_groups = [g for g in order_groups if g in df_long["Group"].unique()]

    fig, axes = plt.subplots(2, 3, figsize=(16, 9), sharey=True)
    axes = axes.flatten()

    for idx, (_, pair_label) in enumerate(pairs):
        ax = axes[idx]
        sub_df = df_long[df_long["Comparison"] == pair_label]

        sns.boxplot(
            data=sub_df,
            x="Group",
            y="nRF",
            order=order_groups,
            palette=[GROUP_COLORS.get(g, "#555") for g in order_groups],
            ax=ax,
            width=0.45,
            boxprops=dict(alpha=0.7),
            showmeans=True,
            meanprops={"marker": "o", "markerfacecolor": "white", "markeredgecolor": "black"}
        )
        sns.stripplot(
            data=sub_df,
            x="Group",
            y="nRF",
            order=order_groups,
            color="black",
            alpha=0.3,
            size=4,
            jitter=0.2,
            ax=ax
        )

        ax.set_title(pair_label, fontweight="bold")
        ax.set_xlabel("")
        ax.set_ylabel("Normalized RF (nRF)" if idx % 3 == 0 else "")
        ax.set_ylim(-0.05, 1.05)

    plt.suptitle("Topological Disagreement (nRF) Across Biological Ortholog Groups (400 HOGs)", fontweight="bold", y=0.98)
    plt.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out_path}")


def plot_method_heatmaps(df_concordance: pd.DataFrame, out_path: Path):
    """Plots 2x2 grid of average pairwise nRF heatmaps for each biological group."""
    methods = ["msa_ml", "msa_nj", "psa_nj", "gs"]
    labels = ["MSA+ML", "MSA+NJ", "PSA+NJ", "GS"]

    groups = ["Amniota", "Eukaryota", "Actinobacteria", "Bacteria"]
    groups = [g for g in groups if g in df_concordance["dataset_name"].unique()]

    fig, axes = plt.subplots(1, len(groups), figsize=(4.5 * len(groups), 4), sharey=True)
    if len(groups) == 1:
        axes = [axes]

    for idx, g_name in enumerate(groups):
        sub = df_concordance[df_concordance["dataset_name"] == g_name]
        mat = np.zeros((4, 4))

        for i in range(4):
            for j in range(4):
                if i == j:
                    mat[i, j] = 0.0
                elif i < j:
                    col = f"nrf_{methods[i]}_vs_{methods[j]}"
                    val = sub[col].mean()
                    mat[i, j] = val
                    mat[j, i] = val

        ax = axes[idx]
        sns.heatmap(
            mat,
            annot=True,
            fmt=".3f",
            cmap="YlOrRd",
            vmin=0.0,
            vmax=0.8,
            xticklabels=labels,
            yticklabels=labels if idx == 0 else False,
            cbar=(idx == len(groups) - 1),
            cbar_kws={"label": "Mean nRF"} if idx == len(groups) - 1 else {},
            ax=ax
        )
        ax.set_title(f"{g_name} (N={len(sub)})", fontweight="bold")

    plt.suptitle("Mean Pairwise Topological Distance (nRF) Heatmaps", fontweight="bold", y=1.02)
    plt.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out_path}")


def plot_sss_vs_nrf_correlation(df_concordance: pd.DataFrame, out_path: Path):
    """Correlates Sequence Similarity Score (SSS) with method disagreement (nRF)."""
    hog_sim_file = RESULTS_DIR / "hog_similarity_summary.csv"
    if not hog_sim_file.exists():
        print(f"Skipping SSS correlation: {hog_sim_file} not found.")
        return

    df_sim = pd.read_csv(hog_sim_file)
    # Merge on dataset and hog_id
    df_merged = pd.merge(df_concordance, df_sim, on=["dataset", "hog_id"], how="inner")
    if df_merged.empty or "mean_sss" not in df_merged.columns:
        return

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    comparisons = [
        ("nrf_msa_ml_vs_gs", "MSA+ML vs GS"),
        ("nrf_msa_ml_vs_psa_nj", "MSA+ML vs PSA+NJ"),
        ("nrf_psa_nj_vs_gs", "PSA+NJ vs GS")
    ]

    for idx, (col, title) in enumerate(comparisons):
        ax = axes[idx]
        for g_name, grp in df_merged.groupby("dataset_name"):
            ax.scatter(
                grp["mean_sss"],
                grp[col],
                label=g_name,
                color=GROUP_COLORS.get(g_name, "#555"),
                alpha=0.6,
                s=40,
                edgecolors="none"
            )

        # Overall trendline
        valid = df_merged.dropna(subset=["mean_sss", col])
        if len(valid) > 2:
            sns.regplot(
                data=valid,
                x="mean_sss",
                y=col,
                scatter=False,
                ax=ax,
                color="black",
                line_kws={"linestyle": "--", "linewidth": 2, "label": "Trendline"}
            )

        ax.set_title(title, fontweight="bold")
        ax.set_xlabel("Mean Sequence Similarity Score (SSS: $\\bar{w}$)")
        ax.set_ylabel("Normalized RF (nRF)")
        ax.set_ylim(-0.05, 1.05)
        ax.legend(frameon=True, loc="upper right")

    plt.suptitle("Impact of Sequence Similarity on Phylogenetic Method Divergence", fontweight="bold", y=1.02)
    plt.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out_path}")


def main():
    concordance_file = RESULTS_DIR / "method_concordance_summary.csv"
    if not concordance_file.exists():
        print(f"Error: {concordance_file} not found.")
        sys.exit(1)

    df_concordance = pd.read_csv(concordance_file)
    print(f"Loaded {len(df_concordance)} records from {concordance_file}")

    violins_png = RESULTS_DIR / "method_concordance_violins.png"
    heatmaps_png = RESULTS_DIR / "method_concordance_heatmaps.png"
    corr_png = RESULTS_DIR / "sss_vs_method_divergence.png"

    plot_method_concordance_violins(df_concordance, violins_png)
    plot_method_heatmaps(df_concordance, heatmaps_png)
    plot_sss_vs_nrf_correlation(df_concordance, corr_png)
    print("All benchmark plots generated successfully!")


if __name__ == "__main__":
    main()
