#!/usr/bin/env python3
"""
plot_cluster_tcs_comparison.py
==============================
Compares Taxonomy Congruence Scores (TCS), win rates, and ranks
among phylogenetic methods (MSA+ML, MSA+NJ, PSA+NJ, GS) across the 4 branch length clusters:
  Cluster 1: Ultra-diverged / Broad Long-branch (Bacteria dominated, mean branch 0.54)
  Cluster 2: Moderate Diverged / Progressive (Bacteria/Actino, mean branch 0.36)
  Cluster 3: Conservative / Short-branch (Amniota/Eukaryota, mean branch 0.10)
  Cluster 4: Long-branch Heterogeneous / LBA Risk (Eukaryota/Actino, max branch 3.14, skew 4.23)

Outputs:
1. cluster_tcs_summary.csv
2. plot_cluster_tcs_comparison.png (TCS distributions, win rates, and average ranks)
"""

from pathlib import Path
import numpy as np
import pandas as pd
import scipy.stats as stats
import matplotlib.pyplot as plt
import seaborn as sns

SCRIPT_DIR = Path(__file__).resolve().parent
BASE_DIR = SCRIPT_DIR.parent
RESULTS_DIR = BASE_DIR / "results"

plt.rcParams.update({
    "font.sans-serif": ["Arial", "DejaVu Sans", "Helvetica", "sans-serif"],
    "font.size": 11,
    "axes.labelsize": 12,
    "axes.titlesize": 13,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 10,
    "figure.titlesize": 15,
    "lines.linewidth": 1.5,
    "axes.grid": True,
    "grid.alpha": 0.3,
    "grid.linestyle": "--"
})

CLUSTER_NAMES = {
    1: "C1: Ultra-diverged\n(n=31, All-Long)",
    2: "C2: Mod. Diverged\n(n=113, Balanced)",
    3: "C3: Conservative\n(n=167, Short)",
    4: "C4: LBA-Risk\n(n=89, Extreme Long)"
}

METHOD_COLORS = {
    "msa_ml": "#2b5c8f",  # Navy
    "msa_nj": "#e07a5f",  # Coral
    "psa_nj": "#81b29a",  # Sage
    "gs": "#f2cc8f"       # Gold
}

METHOD_LABELS = {
    "msa_ml": "MSA+ML (IQ-TREE 2)",
    "msa_nj": "MSA+NJ (RapidNJ)",
    "psa_nj": "PSA+NJ (RapidNJ)",
    "gs": "GS (Graph Splitting)"
}


def main():
    df = pd.read_csv(RESULTS_DIR / "branch_length_clusters.csv")
    methods = ["msa_ml", "msa_nj", "psa_nj", "gs"]
    
    # 1. Compute per-HOG win shares and ranks
    def evaluate_hog(row):
        scores = [row[f"tcs_{m}"] for m in methods]
        max_s = max(scores)
        winners = [m for m in methods if row[f"tcs_{m}"] == max_s]
        w = 1.0 / len(winners)
        
        ranks = stats.rankdata([-s for s in scores])
        res = {}
        for idx, m in enumerate(methods):
            res[f"win_{m}"] = w if m in winners else 0.0
            res[f"rank_{m}"] = ranks[idx]
        return pd.Series(res)

    eval_df = df.apply(evaluate_hog, axis=1)
    df_full = pd.concat([df, eval_df], axis=1)

    # 2. Build summary table
    summary_rows = []
    for cl in sorted(df["cluster"].unique()):
        sub = df_full[df_full["cluster"] == cl]
        n_hogs = len(sub)
        
        row = {
            "cluster": cl,
            "cluster_name": CLUSTER_NAMES[cl].replace("\n", " "),
            "num_hogs": n_hogs,
            "mean_sss": sub["mean_sss"].mean(),
            "mean_branch_ml": sub["msa_ml_mean"].mean(),
            "max_branch_ml": sub["msa_ml_max"].mean(),
            "skew_ml": sub["msa_ml_skew"].mean()
        }
        
        for m in methods:
            row[f"mean_tcs_{m}"] = sub[f"tcs_{m}"].mean()
            row[f"median_tcs_{m}"] = sub[f"tcs_{m}"].median()
            row[f"win_rate_{m}"] = sub[f"win_{m}"].mean() * 100
            row[f"avg_rank_{m}"] = sub[f"rank_{m}"].mean()

        # Head to head: ML vs NJ, ML vs GS, PSA vs MSA
        row["h2h_ml_over_msa_nj%"] = (sub["tcs_msa_ml"] > sub["tcs_msa_nj"]).mean() * 100
        row["h2h_ml_over_gs%"] = (sub["tcs_msa_ml"] > sub["tcs_gs"]).mean() * 100
        row["h2h_psa_over_msa_nj%"] = (sub["tcs_psa_nj"] > sub["tcs_msa_nj"]).mean() * 100
        
        summary_rows.append(row)

    summary_df = pd.DataFrame(summary_rows)
    out_csv = RESULTS_DIR / "cluster_tcs_summary.csv"
    summary_df.to_csv(out_csv, index=False)
    print(f"Saved cluster TCS summary to {out_csv}")

    # 3. Create Multi-panel Visualization
    fig, axes = plt.subplots(1, 3, figsize=(20, 6.5))

    # (A) Win Rate per Cluster (Grouped Bar Chart)
    ax0 = axes[0]
    bar_width = 0.2
    x = np.arange(len(summary_df))

    for idx, m in enumerate(methods):
        win_rates = summary_df[f"win_rate_{m}"]
        bars = ax0.bar(
            x + (idx - 1.5) * bar_width,
            win_rates,
            width=bar_width,
            label=METHOD_LABELS[m],
            color=METHOD_COLORS[m],
            alpha=0.9,
            edgecolor="black",
            linewidth=0.8
        )
        for bar in bars:
            h = bar.get_height()
            ax0.text(
                bar.get_x() + bar.get_width() / 2.,
                h + 0.8,
                f"{h:.1f}%",
                ha="center",
                va="bottom",
                fontsize=8.5,
                fontweight="bold"
            )

    ax0.set_title("A. Win Rate (% of Top TCS) by Cluster", fontsize=13, fontweight="bold")
    ax0.set_xticks(x)
    ax0.set_xticklabels([CLUSTER_NAMES[c] for c in summary_df["cluster"]], fontsize=10)
    ax0.set_ylabel("Win Rate within HOG (%)")
    ax0.set_ylim(0, 50)
    ax0.legend(loc="upper left", frameon=True, fontsize=9.5)

    # (B) Average Rank per Cluster (Inverted, 1 is best)
    ax1 = axes[1]
    for idx, m in enumerate(methods):
        ranks = summary_df[f"avg_rank_{m}"]
        ax1.plot(
            x,
            ranks,
            marker="o",
            markersize=8,
            linewidth=2.2,
            label=METHOD_LABELS[m],
            color=METHOD_COLORS[m]
        )
        for xi, yi in zip(x, ranks):
            ax1.text(xi, yi - 0.06, f"{yi:.2f}", ha="center", va="top", fontsize=9, fontweight="bold", color=METHOD_COLORS[m])

    ax1.set_title("B. Average Rank by Cluster (Lower is Better)", fontsize=13, fontweight="bold")
    ax1.set_xticks(x)
    ax1.set_xticklabels([CLUSTER_NAMES[c] for c in summary_df["cluster"]], fontsize=10)
    ax1.set_ylabel("Average Rank (1 = Top, 4 = Worst)")
    ax1.set_ylim(3.3, 1.9)  # Invert y-axis so 1 is at top
    ax1.legend(loc="lower left", frameon=True, fontsize=9.5)

    # (C) Mean TCS by Cluster
    ax2 = axes[2]
    for idx, m in enumerate(methods):
        means = summary_df[f"mean_tcs_{m}"]
        bars = ax2.bar(
            x + (idx - 1.5) * bar_width,
            means,
            width=bar_width,
            label=METHOD_LABELS[m],
            color=METHOD_COLORS[m],
            alpha=0.9,
            edgecolor="black",
            linewidth=0.8
        )
        for bar in bars:
            h = bar.get_height()
            ax2.text(
                bar.get_x() + bar.get_width() / 2.,
                h + 0.4,
                f"{h:.1f}",
                ha="center",
                va="bottom",
                fontsize=8.5
            )

    ax2.set_title("C. Mean TCS Score by Cluster", fontsize=13, fontweight="bold")
    ax2.set_xticks(x)
    ax2.set_xticklabels([CLUSTER_NAMES[c] for c in summary_df["cluster"]], fontsize=10)
    ax2.set_ylabel("Mean TCS (Taxonomy Congruence Score)")
    ax2.set_ylim(0, 32)
    ax2.legend(loc="upper left", frameon=True, fontsize=9.5)

    plt.suptitle("Phylogenetic Method Performance (TCS) Across Branch Length Clusters", fontsize=15, y=1.02)
    plt.tight_layout()
    out_png = RESULTS_DIR / "plot_cluster_tcs_comparison.png"
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved visualization plot to {out_png}")

    # Print nicely formatted summary table
    print("\n--- CLUSTER TCS SUMMARY TABLE ---")
    cols_display = [
        "cluster", "num_hogs", "mean_sss", "max_branch_ml",
        "win_rate_msa_ml", "win_rate_msa_nj", "win_rate_psa_nj", "win_rate_gs",
        "avg_rank_msa_ml", "avg_rank_psa_nj", "avg_rank_gs"
    ]
    print(summary_df[cols_display].to_string(index=False))


if __name__ == "__main__":
    main()
