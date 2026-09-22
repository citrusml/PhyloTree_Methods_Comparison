#!/usr/bin/env python3
"""
plot_group_win_rates.py
=======================
Plots the within-HOG win rates (%) and mean ranks across biological groups
(Amniota, Eukaryota, Actinobacteria, Bacteria) and Overall (Total 400 HOGs).
Corresponds to Table 6.1 in memo/20260911_dataset_results.md.

Output:
  dataset/results/plot_group_win_rates.png
"""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent.parent
DATASET_DIR = SCRIPT_DIR.parent
RESULTS_DIR = DATASET_DIR / "results"

plt.rcParams.update({
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
    "font.size": 11,
    "axes.labelsize": 12,
    "axes.titlesize": 13,
    "xtick.labelsize": 10.5,
    "ytick.labelsize": 10.5,
    "legend.fontsize": 10,
    "figure.titlesize": 14.5,
    "axes.grid": True,
    "grid.alpha": 0.3,
    "grid.linestyle": "--"
})

METHODS = [
    ("msa_ml", "MSA+ML (IQ-TREE 2)", "#2980b9"),
    ("msa_nj", "MSA+NJ (RapidNJ)",   "#e67e22"),
    ("psa_nj", "PSA+NJ (RapidNJ)",   "#27ae60"),
    ("gs",     "GS (Graph Splitting)", "#c0392b")
]

GROUPS = [
    ("amn", "Amniota\n(N=100)"),
    ("euk", "Eukaryota\n(N=100)"),
    ("act", "Actinobacteria\n(N=100)"),
    ("bac", "Bacteria\n(N=100)"),
    ("all", "Total (All Groups)\n(N=400)")
]


def main():
    csv_path = RESULTS_DIR / "within_hog_win_rates.csv"
    if not csv_path.exists():
        raise FileNotFoundError(f"{csv_path} not found.")
    df = pd.read_csv(csv_path)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))

    x = np.arange(len(GROUPS))
    width = 0.20

    # (A) Win Rate %
    for i, (m_key, m_name, color) in enumerate(METHODS):
        sub = df[df["Method_Key"] == m_key].set_index("Group").reindex([g[0] for g in GROUPS])
        vals = sub["Win_Rate_%"].values
        bars = ax1.bar(x + (i - 1.5) * width, vals, width, label=m_name, color=color, alpha=0.9, edgecolor="black", linewidth=0.6)
        for bar in bars:
            h = bar.get_height()
            ax1.annotate(f"{h:.1f}%",
                         xy=(bar.get_x() + bar.get_width() / 2, h),
                         xytext=(0, 3), textcoords="offset points",
                         ha="center", va="bottom", fontsize=8.5, fontweight="bold")

    ax1.set_title("(A) Win Rate (% HOGs where Method is Best in TCS)", fontweight="bold", fontsize=12.5)
    ax1.set_ylabel("Win Rate (%) within Same HOG", fontsize=11.5)
    ax1.set_xticks(x)
    ax1.set_xticklabels([g[1] for g in GROUPS], fontsize=10.5)
    ax1.set_ylim(0, 52)
    ax1.axvline(3.5, color="gray", linestyle=":", linewidth=1.5, alpha=0.8)
    ax1.text(3.6, 48, "Summary", fontsize=10, style="italic", color="#555555")
    ax1.legend(frameon=True, loc="upper right", fontsize=10)

    # (B) Mean Rank
    for i, (m_key, m_name, color) in enumerate(METHODS):
        sub = df[df["Method_Key"] == m_key].set_index("Group").reindex([g[0] for g in GROUPS])
        vals = sub["Mean_Rank"].values
        bars = ax2.bar(x + (i - 1.5) * width, vals, width, label=m_name, color=color, alpha=0.9, edgecolor="black", linewidth=0.6)
        for bar in bars:
            h = bar.get_height()
            ax2.annotate(f"{h:.2f}",
                         xy=(bar.get_x() + bar.get_width() / 2, h),
                         xytext=(0, 3), textcoords="offset points",
                         ha="center", va="bottom", fontsize=8.5, fontweight="bold")

    ax2.set_title("(B) Mean Rank within HOG (1.0 = Best, 4.0 = Worst)", fontweight="bold", fontsize=12.5)
    ax2.set_ylabel("Mean Rank (Lower is Better)", fontsize=11.5)
    ax2.set_xticks(x)
    ax2.set_xticklabels([g[1] for g in GROUPS], fontsize=10.5)
    ax2.set_ylim(1.0, 4.0)
    ax2.axvline(3.5, color="gray", linestyle=":", linewidth=1.5, alpha=0.8)
    ax2.text(3.6, 3.85, "Summary", fontsize=10, style="italic", color="#555555")
    ax2.legend(frameon=True, loc="upper right", fontsize=10)

    plt.suptitle("Within-HOG Performance by Biological Group (TCS Evaluation Across 400 Real HOGs)", fontweight="bold", fontsize=14, y=0.98)
    plt.tight_layout()

    out_file = RESULTS_DIR / "plot_group_win_rates.png"
    fig.savefig(out_file, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out_file}")


if __name__ == "__main__":
    main()
