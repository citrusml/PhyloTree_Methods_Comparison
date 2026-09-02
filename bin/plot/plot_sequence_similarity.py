#!/usr/bin/env python3
"""
Plot Sequence Similarity Analysis
Visualizes:
1. Sequence Identity vs Evolutionary Distance D (across lengths L) with clear legend
2. Comparison between True MSA Identity vs MAFFT MSA Identity vs PSA Identity (Alignment Bias / Over-alignment)
3. Gap Proportion vs Evolutionary Distance D (True MSA vs MAFFT MSA across lengths L) with clear legend
"""

import os
import sys
from pathlib import Path

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

PROJECT_ROOT = Path(__file__).resolve().parent.parent
ANALYSIS_DIR = PROJECT_ROOT / "analysis"
SIMILARITY_DIR = ANALYSIS_DIR / "similarity"
CSV_PATH = SIMILARITY_DIR / "similarity_summary.csv"
BENCHMARK_PATH = ANALYSIS_DIR / "all_experiments_summary.csv"

# Global style settings
plt.style.use('seaborn-v0_8-whitegrid' if 'seaborn-v0_8-whitegrid' in plt.style.available else 'default')
plt.rcParams.update({
    'font.size': 11,
    'figure.titlesize': 14,
    'axes.titlesize': 12,
    'axes.labelsize': 11,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 10,
    'legend.title_fontsize': 11,
    'savefig.dpi': 300
})

def plot_similarity_analysis():
    if not CSV_PATH.exists():
        print(f"Error: {CSV_PATH} does not exist.")
        return

    df_sim = pd.read_csv(CSV_PATH)
    print(f"Loaded {len(df_sim)} similarity records.")

    # =========================================================================
    # 【図1】 配列同一性（% Identity）の減衰スペクトラム (D = 0.1 ~ 6.0)
    # =========================================================================
    lengths = [300, 500, 1000, 1500]
    fig, axes = plt.subplots(1, len(lengths), figsize=(16, 4.2), sharey=True)

    for idx, length in enumerate(lengths):
        ax = axes[idx]
        sub = df_sim[df_sim['length'] == length]

        # 3つのアライメント基準での同一性を描画
        sns.lineplot(data=sub, x='distance', y='true_identity_mean', ax=ax,
                     color='#9467bd', marker='o', label='True MSA Identity', errorbar=('ci', 95), linewidth=1.8)
        sns.lineplot(data=sub, x='distance', y='mafft_identity_mean', ax=ax,
                     color='#ff7f0e', marker='s', label='MAFFT MSA Identity', errorbar=('ci', 95), linewidth=1.8)
        sns.lineplot(data=sub, x='distance', y='psa_identity_mean', ax=ax,
                     color='#1f77b4', marker='^', label='PSA Identity (Needleman-Wunsch)', errorbar=('ci', 95), linewidth=1.8)

        # ランダム期待値 5% (1/20) の水平線
        #ax.axhline(0.05, color='gray', linestyle=':', linewidth=1.2, label='Random Expectation (5%)')
        # Twilight zone (20-35%) の帯
        #ax.axhspan(0.20, 0.35, color='orange', alpha=0.15, label='Twilight Zone (20-35%)')

        ax.set_title(f"Length L = {length} aa", fontweight='bold')
        ax.set_xlabel("Evolutionary Distance (D)")
        ax.set_ylabel("Pairwise Sequence Identity" if idx == 0 else "")
        ax.set_ylim(-0.02, 1.02)

        # 各サブプロット内部の凡例はいったん削除
        if ax.get_legend():
            ax.get_legend().remove()

    # 1番目のサブプロットからハンドルとラベルを抽出し、1番右のサブプロットの外側に配置
    handles, labels = axes[0].get_legend_handles_labels()
    axes[-1].legend(handles, labels, bbox_to_anchor=(1.05, 1.0), loc='upper left', frameon=True, title="Identity Metric")

    plt.subplots_adjust(wspace=0.08, top=0.85, right=0.82)
    plt.suptitle("Sequence Similarity Decay Spectrum across Evolutionary Distance (D = 0.1 to 6.0)", 
                 fontweight='bold', y=0.98)
    
    out_fig1 = SIMILARITY_DIR / "sequence_identity_spectrum.png"
    plt.savefig(out_fig1, bbox_inches='tight')
    plt.close()
    print(f"Saved: {out_fig1}")

    # =========================================================================
    # 【図2】 アライメントバイアス（MAFFT / PSA による見かけ上の一致度過大評価）
    # =========================================================================
    df_sim['mafft_overalignment'] = df_sim['mafft_identity_mean'] - df_sim['true_identity_mean']
    df_sim['psa_overalignment'] = df_sim['psa_identity_mean'] - df_sim['true_identity_mean']

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))

    pivot_mafft = df_sim.pivot_table(index='distance', columns='length', values='mafft_overalignment', aggfunc='mean')
    sns.heatmap(pivot_mafft, ax=ax1, annot=True, fmt='+.3f', cmap='coolwarm', center=0,
                cbar_kws={'label': 'Mean Identity Difference (MAFFT - True)'})
    ax1.set_title("MAFFT MSA Over-alignment Bias\\n(MAFFT Identity - True Identity)", fontweight='bold')
    ax1.set_xlabel("Sequence Length (L)")
    ax1.set_ylabel("Evolutionary Distance (D)")
    ax1.invert_yaxis()

    pivot_psa = df_sim.pivot_table(index='distance', columns='length', values='psa_overalignment', aggfunc='mean')
    sns.heatmap(pivot_psa, ax=ax2, annot=True, fmt='+.3f', cmap='coolwarm', center=0,
                cbar_kws={'label': 'Mean Identity Difference (PSA - True)'})
    ax2.set_title("Needleman-Wunsch PSA Over-alignment Bias\\n(PSA Identity - True Identity)", fontweight='bold')
    ax2.set_xlabel("Sequence Length (L)")
    ax2.set_ylabel("Evolutionary Distance (D)")
    ax2.invert_yaxis()

    plt.suptitle("Over-alignment Bias in High Evolutionary Distance Regime (D = 0.1 to 6.0)", 
                 fontweight='bold', y=1.02)
    plt.tight_layout()
    
    out_fig2 = SIMILARITY_DIR / "overalignment_bias_heatmap.png"
    plt.savefig(out_fig2, bbox_inches='tight')
    plt.close()
    print(f"Saved: {out_fig2}")

    # =========================================================================
    # 【図3】 ギャップ割合（Gap Proportion）の推移 (True MSA vs MAFFT MSA)
    # =========================================================================
    fig, axes = plt.subplots(1, len(lengths), figsize=(16, 4.2), sharey=True)

    for idx, length in enumerate(lengths):
        ax = axes[idx]
        sub = df_sim[df_sim['length'] == length]

        # 1. True MSA のギャップ率
        sns.lineplot(data=sub, x='distance', y='true_gap_mean', ax=ax,
                     color='#9467bd', marker='o', label='True MSA Gap Proportion', 
                     errorbar=('ci', 95), linewidth=1.8)
        # 2. MAFFT MSA のギャップ率
        sns.lineplot(data=sub, x='distance', y='mafft_gap_mean', ax=ax,
                     color='#ff7f0e', marker='s', linestyle='--', label='MAFFT MSA Gap Proportion', 
                     errorbar=('ci', 95), linewidth=1.8)

        ax.set_title(f"Length L = {length} aa", fontweight='bold')
        ax.set_xlabel("Evolutionary Distance (D)")
        ax.set_ylabel("Mean Gap Proportion per Taxon" if idx == 0 else "")
        ax.set_ylim(-0.02, 1.02)

        if ax.get_legend():
            ax.get_legend().remove()

    handles, labels = axes[0].get_legend_handles_labels()
    axes[-1].legend(handles, labels, bbox_to_anchor=(1.05, 1.0), loc='upper left', frameon=True, title="Gap Metric")

    plt.subplots_adjust(wspace=0.08, top=0.85, right=0.82)
    plt.suptitle("Mean Gap Proportion across Evolutionary Distance (True MSA vs MAFFT MSA, D = 0.1 to 6.0)", 
                 fontweight='bold', y=0.98)

    out_fig3 = SIMILARITY_DIR / "gap_proportion_vs_distance.png"
    plt.savefig(out_fig3, bbox_inches='tight')
    plt.close()
    print(f"Saved: {out_fig3}")

    print("\n🎉 All sequence similarity plots generated successfully!")

if __name__ == "__main__":
    plot_similarity_analysis()
