#!/usr/bin/env python3
"""
Comprehensive Post-Processing, Reporting, and Artifact Organization Script: generate_all_reports.py
Integrates all analysis, visualization, sequence length tracking, and replication file organization:

1. Phylogenetic Tree Benchmark Analysis:
   - 2D Regime Map Heatmap (ΔnRF: PWA+NJ vs MSA+ML)
   - Pairwise Method Comparisons (ML vs NJ, NJ vs PWA, PWA vs ML, 3-panel comparison)
   - Condition-wise nRF Boxplots (5x5 grid with mean μ and variance σ²)
   - Summary Statistics CSV

2. Sequence Length & Residue Retention Analysis:
   - Remaining sequence length distributions under indels (violin plots & condition boxplots)
   - Sequence length summary statistics table (CSV)

3. Replication Artifact Harvesting & Organization:
   - Collects true trees, unaligned fasta, true MSA, inferred trees, and matrices from work/ into results/replications/
"""

import os
import re
import sys
import glob
import json
import shutil
import argparse
from concurrent.futures import ThreadPoolExecutor
from collections import defaultdict
from typing import Dict, List, Tuple, Optional, Any, Union
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
    "PWA+NJ": "#1f77b4",       # Blue
    "MSA+NJ": "#ff7f0e",       # Orange
    "MSA+ML": "#2ca02c",       # Green (IQ-TREE 2)
    "MSA+RAXML": "#006400",    # Dark Green (RAxML -f d)
    "MSA+BI": "#d62728",       # Red (MrBayes Bayesian Inference)
    "GS": "#9467bd",           # Purple (Matsui & Iwasaki 2020 Graph Splitting)
    "TRUE_PWA+NJ": "#17becf",  # Cyan
    "TRUE_MSA+NJ": "#8c564b",  # Brown
    "TRUE_MSA+ML": "#e377c2",  # Pink (True IQ-TREE 2)
    "TRUE_MSA+RAXML": "#c49c94", # Rosy brown (True RAxML -f d)
    "TRUE_MSA+BI": "#bcbd22",  # Olive
    "TRUE_DIST+NJ": "#7f7f7f", # Grey
    "PWA+FastME": "#e377c2",
    "MSA+FastME": "#bcbd22",
    "PWA+FastME_SPR": "#e377c2",
    "MSA+FastME_LG_G": "#bcbd22",
}

# ==========================================
# 1. Phylogenetic Benchmark Visualizations
# ==========================================

def generate_regime_map(df, outdir):
    """Generates 2D phase diagram of Distance vs Length showing ΔnRF = PWA+NJ - MSA+ML."""
    pivot_df = df.groupby(["distance", "length", "pipeline"])["nrf_distance"].mean().reset_index()

    pwa = pivot_df[pivot_df["pipeline"] == "PWA+NJ"].set_index(["length", "distance"])["nrf_distance"]
    msa_ml = pivot_df[pivot_df["pipeline"] == "MSA+ML"].set_index(["length", "distance"])["nrf_distance"]

    if pwa.empty or msa_ml.empty:
        print("Note: PWA+NJ or MSA+ML missing in dataset. Skipping primary regime map.")
        return

    delta_rf = (pwa - msa_ml).unstack(level="distance")

    plt.figure(figsize=(9, 7))
    sns.set_theme(style="whitegrid", font_scale=1.1)
    cmap = sns.diverging_palette(240, 10, as_cmap=True)

    ax = sns.heatmap(
        delta_rf,
        annot=True,
        fmt=".3f",
        cmap=cmap,
        center=0,
        cbar_kws={'label': r'$\Delta nRF = nRF_{\text{PWA+NJ}} - nRF_{\text{MSA+ML}}$'}
    )

    plt.title("Regime Map: PWA+NJ vs MSA+ML Relative Performance", fontsize=14, fontweight="bold", pad=15)
    plt.xlabel("Evolutionary Distance D (substitutions/site)", fontsize=12)
    plt.ylabel("Sequence Length L (amino acids)", fontsize=12)

    plot_path = os.path.join(outdir, "regime_map_delta_nrf.png")
    plt.tight_layout()
    plt.savefig(plot_path, dpi=300)
    plt.close()
    print(f"Generated: {plot_path}")

def generate_method_comparisons(df, outdir):
    """Generates comprehensive pairwise comparison heatmaps and absolute performance heatmaps across all methods."""
    comp_dir = os.path.join(outdir, "method_comparisons")
    os.makedirs(comp_dir, exist_ok=True)

    pivot_df = df.groupby(["distance", "length", "pipeline"])["nrf_distance"].mean().reset_index()
    pipes = df["pipeline"].unique()

    def get_grid(pname):
        sub = pivot_df[pivot_df["pipeline"] == pname]
        if sub.empty:
            return None
        return sub.set_index(["length", "distance"])["nrf_distance"].unstack(level="distance")

    # -------------------------------------------------------------
    # 1. Standard Inferred Pipelines Comparison
    # -------------------------------------------------------------
    # ML vs NJ
    if "MSA+ML" in pipes and "MSA+NJ" in pipes:
        grid_ml = get_grid("MSA+ML")
        grid_nj = get_grid("MSA+NJ")
        delta = grid_nj - grid_ml  # >0: ML wins (red), <0: NJ wins (blue)
        _plot_heatmap(delta, "Reconstruction Engine: MSA+NJ vs MSA+ML\n(ΔnRF = MSA+NJ - MSA+ML)",
                      r"$\Delta nRF$ (>0: ML wins, <0: NJ wins)", os.path.join(comp_dir, "comparison_ml_vs_nj.png"))

    # NJ vs PWA
    if "MSA+NJ" in pipes and "PWA+NJ" in pipes:
        grid_pwa = get_grid("PWA+NJ")
        grid_nj = get_grid("MSA+NJ")
        delta = grid_pwa - grid_nj  # >0: MSA wins (red), <0: PWA wins (blue)
        _plot_heatmap(delta, "Alignment Strategy: PWA+NJ vs MSA+NJ\n(ΔnRF = PWA+NJ - MSA+NJ)",
                      r"$\Delta nRF$ (>0: MSA wins, <0: PWA wins)", os.path.join(comp_dir, "comparison_nj_vs_pwa.png"))

    # PWA vs ML
    if "PWA+NJ" in pipes and "MSA+ML" in pipes:
        grid_pwa = get_grid("PWA+NJ")
        grid_ml = get_grid("MSA+ML")
        delta = grid_pwa - grid_ml
        _plot_heatmap(delta, "Overall Pipeline: PWA+NJ vs MSA+ML\n(ΔnRF = PWA+NJ - MSA+ML)",
                      r"$\Delta nRF$ (>0: ML wins, <0: PWA wins)", os.path.join(comp_dir, "comparison_pwa_vs_ml.png"))

    # -------------------------------------------------------------
    # 2. True MSA Benchmark Comparisons (Experiment 1 / Alignment Error Analysis)
    # -------------------------------------------------------------
    # Alignment Error Penalty on NJ (MSA+NJ vs TRUE_MSA+NJ)
    if "MSA+NJ" in pipes and "TRUE_MSA+NJ" in pipes:
        grid_msa_nj = get_grid("MSA+NJ")
        grid_true_nj = get_grid("TRUE_MSA+NJ")
        delta = grid_msa_nj - grid_true_nj  # >0: Alignment error degraded NJ accuracy
        _plot_heatmap(delta, "Alignment Error Penalty (NJ): MSA+NJ vs TRUE_MSA+NJ\n(ΔnRF = MSA+NJ - TRUE_MSA+NJ)",
                      r"$\Delta nRF$ (>0: Alignment Error Degraded Tree)", os.path.join(comp_dir, "comparison_alignment_loss_nj.png"))

    # Alignment Error Penalty on ML (MSA+ML vs TRUE_MSA+ML)
    if "MSA+ML" in pipes and "TRUE_MSA+ML" in pipes:
        grid_msa_ml = get_grid("MSA+ML")
        grid_true_ml = get_grid("TRUE_MSA+ML")
        delta = grid_msa_ml - grid_true_ml  # >0: Alignment error degraded ML accuracy
        _plot_heatmap(delta, "Alignment Error Penalty (ML): MSA+ML vs TRUE_MSA+ML\n(ΔnRF = MSA+ML - TRUE_MSA+ML)",
                      r"$\Delta nRF$ (>0: Alignment Error Degraded Tree)", os.path.join(comp_dir, "comparison_alignment_loss_ml.png"))

    # PWA+NJ vs TRUE_PWA+NJ
    if "PWA+NJ" in pipes and "TRUE_PWA+NJ" in pipes:
        grid_pwa = get_grid("PWA+NJ")
        grid_true_pwa = get_grid("TRUE_PWA+NJ")
        delta = grid_pwa - grid_true_pwa
        _plot_heatmap(delta, "PWA vs True PWA Alignment: PWA+NJ vs TRUE_PWA+NJ\n(ΔnRF = PWA+NJ - TRUE_PWA+NJ)",
                      r"$\Delta nRF$ (>0: True PWA wins, <0: PWA wins)", os.path.join(comp_dir, "comparison_pwa_vs_true_pwa_nj.png"))

    # PWA+NJ vs TRUE_MSA+NJ (How close PWA gets to ideal True MSA NJ)
    if "PWA+NJ" in pipes and "TRUE_MSA+NJ" in pipes:
        grid_pwa = get_grid("PWA+NJ")
        grid_true_nj = get_grid("TRUE_MSA+NJ")
        delta = grid_pwa - grid_true_nj
        _plot_heatmap(delta, "PWA vs Ideal Alignment: PWA+NJ vs TRUE_MSA+NJ\n(ΔnRF = PWA+NJ - TRUE_MSA+NJ)",
                      r"$\Delta nRF$ (>0: True MSA wins, <0: PWA wins)", os.path.join(comp_dir, "comparison_pwa_vs_true_msa_nj.png"))

    # PWA+NJ vs TRUE_MSA+ML
    if "PWA+NJ" in pipes and "TRUE_MSA+ML" in pipes:
        grid_pwa = get_grid("PWA+NJ")
        grid_true_ml = get_grid("TRUE_MSA+ML")
        delta = grid_pwa - grid_true_ml
        _plot_heatmap(delta, "PWA vs Ideal ML: PWA+NJ vs TRUE_MSA+ML\n(ΔnRF = PWA+NJ - TRUE_MSA+ML)",
                      r"$\Delta nRF$ (>0: True MSA+ML wins, <0: PWA wins)", os.path.join(comp_dir, "comparison_pwa_vs_true_msa_ml.png"))

    # TRUE_MSA+NJ vs TRUE_MSA+ML (Intrinsic algorithm difference without alignment error)
    if "TRUE_MSA+NJ" in pipes and "TRUE_MSA+ML" in pipes:
        grid_true_nj = get_grid("TRUE_MSA+NJ")
        grid_true_ml = get_grid("TRUE_MSA+ML")
        delta = grid_true_nj - grid_true_ml
        _plot_heatmap(delta, "Intrinsic Algorithm Performance: TRUE_MSA+NJ vs TRUE_MSA+ML\n(ΔnRF = TRUE_MSA+NJ - TRUE_MSA+ML)",
                      r"$\Delta nRF$ (>0: ML wins, <0: NJ wins)", os.path.join(comp_dir, "comparison_true_ml_vs_true_nj.png"))

    # GS (Graph Splitting) Comparisons
    if "GS" in pipes:
        grid_gs = get_grid("GS")
        if "PWA+NJ" in pipes:
            delta = grid_pwa - grid_gs
            _plot_heatmap(delta, "Pairwise Benchmark: PWA+NJ vs GS (Graph Splitting)\n(ΔnRF = PWA+NJ - GS)",
                          r"$\Delta nRF$ (>0: GS wins, <0: PWA wins)", os.path.join(comp_dir, "comparison_pwa_vs_gs.png"))
        if "MSA+ML" in pipes:
            delta = grid_ml - grid_gs
            _plot_heatmap(delta, "Regime Map: MSA+ML vs GS (Graph Splitting)\n(ΔnRF = MSA+ML - GS)",
                          r"$\Delta nRF$ (>0: GS wins, <0: ML wins)", os.path.join(comp_dir, "comparison_ml_vs_gs.png"))
        if "MSA+NJ" in pipes:
            delta = grid_nj - grid_gs
            _plot_heatmap(delta, "Regime Map: MSA+NJ vs GS (Graph Splitting)\n(ΔnRF = MSA+NJ - GS)",
                          r"$\Delta nRF$ (>0: GS wins, <0: NJ wins)", os.path.join(comp_dir, "comparison_nj_vs_gs.png"))

    # MSA+BI (MrBayes) Comparisons
    if "MSA+BI" in pipes:
        grid_bi = get_grid("MSA+BI")
        if "MSA+ML" in pipes:
            delta = grid_ml - grid_bi
            _plot_heatmap(delta, "Likelihood vs Bayesian: MSA+ML vs MSA+BI (MrBayes)\n(ΔnRF = MSA+ML - MSA+BI)",
                          r"$\Delta nRF$ (>0: BI wins, <0: ML wins)", os.path.join(comp_dir, "comparison_ml_vs_bi.png"))
        if "PWA+NJ" in pipes:
            delta = grid_pwa - grid_bi
            _plot_heatmap(delta, "Regime Map: PWA+NJ vs MSA+BI (MrBayes)\n(ΔnRF = PWA+NJ - MSA+BI)",
                          r"$\Delta nRF$ (>0: BI wins, <0: PWA wins)", os.path.join(comp_dir, "comparison_pwa_vs_bi.png"))
        if "GS" in pipes:
            delta = grid_gs - grid_bi
            _plot_heatmap(delta, "Graph Splitting vs Bayesian: GS vs MSA+BI (MrBayes)\n(ΔnRF = GS - MSA+BI)",
                          r"$\Delta nRF$ (>0: BI wins, <0: GS wins)", os.path.join(comp_dir, "comparison_gs_vs_bi.png"))

    # RAxML (-f d) Comparisons (Paper syz049 benchmark)
    if "MSA+RAXML" in pipes:
        grid_raxml = get_grid("MSA+RAXML")
        if "MSA+NJ" in pipes:
            delta = grid_nj - grid_raxml
            _plot_heatmap(delta, "Reconstruction Engine: MSA+NJ vs MSA+RAXML\n(ΔnRF = MSA+NJ - MSA+RAXML)",
                          r"$\Delta nRF$ (>0: RAxML wins, <0: NJ wins)", os.path.join(comp_dir, "comparison_raxml_vs_nj.png"))
        if "PWA+NJ" in pipes:
            delta = grid_pwa - grid_raxml
            _plot_heatmap(delta, "Pipeline Comparison: PWA+NJ vs MSA+RAXML\n(ΔnRF = PWA+NJ - MSA+RAXML)",
                          r"$\Delta nRF$ (>0: RAxML wins, <0: PWA wins)", os.path.join(comp_dir, "comparison_pwa_vs_raxml.png"))
        if "MSA+ML" in pipes:
            delta = grid_ml - grid_raxml
            _plot_heatmap(delta, "ML Search Engine: MSA+ML (IQ-TREE) vs MSA+RAXML (RAxML -f d)\n(ΔnRF = IQ-TREE - RAxML)",
                          r"$\Delta nRF$ (>0: RAxML wins, <0: IQ-TREE wins)", os.path.join(comp_dir, "comparison_iqtree_vs_raxml.png"))
        if "GS" in pipes:
            delta = grid_raxml - grid_gs
            _plot_heatmap(delta, "Regime Map: MSA+RAXML vs GS (Graph Splitting)\n(ΔnRF = MSA+RAXML - GS)",
                          r"$\Delta nRF$ (>0: GS wins, <0: RAxML wins)", os.path.join(comp_dir, "comparison_raxml_vs_gs.png"))
        if "TRUE_MSA+RAXML" in pipes:
            grid_true_raxml = get_grid("TRUE_MSA+RAXML")
            delta = grid_raxml - grid_true_raxml
            _plot_heatmap(delta, "Alignment Error Penalty (RAxML): MSA+RAXML vs TRUE_MSA+RAXML\n(ΔnRF = MSA+RAXML - TRUE_MSA+RAXML)",
                          r"$\Delta nRF$ (>0: Alignment Error Degraded Tree)", os.path.join(comp_dir, "comparison_alignment_loss_raxml.png"))


    # -------------------------------------------------------------
    # 3. True Patristic Distance Matrix Benchmark Comparisons
    # -------------------------------------------------------------
    if "TRUE_DIST+NJ" in pipes:
        grid_true_dist = get_grid("TRUE_DIST+NJ")

        # PWA+NJ vs TRUE_DIST+NJ
        if "PWA+NJ" in pipes:
            delta = get_grid("PWA+NJ") - grid_true_dist
            _plot_heatmap(delta, "PWA vs True Distance: PWA+NJ vs TRUE_DIST+NJ\n(ΔnRF = PWA+NJ - TRUE_DIST+NJ)",
                          r"$\Delta nRF$ (>0: True Dist wins, <0: PWA wins)", os.path.join(comp_dir, "comparison_pwa_vs_true_dist.png"))

        # MSA+NJ vs TRUE_DIST+NJ
        if "MSA+NJ" in pipes:
            delta = get_grid("MSA+NJ") - grid_true_dist
            _plot_heatmap(delta, "MSA vs True Distance: MSA+NJ vs TRUE_DIST+NJ\n(ΔnRF = MSA+NJ - TRUE_DIST+NJ)",
                          r"$\Delta nRF$ (>0: True Dist wins, <0: MSA wins)", os.path.join(comp_dir, "comparison_msa_nj_vs_true_dist.png"))

        # TRUE_MSA+NJ vs TRUE_DIST+NJ (Sequence sampling & distance estimation error)
        if "TRUE_MSA+NJ" in pipes:
            delta = get_grid("TRUE_MSA+NJ") - grid_true_dist
            _plot_heatmap(delta, "Distance Estimation Loss: TRUE_MSA+NJ vs TRUE_DIST+NJ\n(ΔnRF = TRUE_MSA+NJ - TRUE_DIST+NJ)",
                          r"$\Delta nRF$ (>0: Sampling/Distance Error)", os.path.join(comp_dir, "comparison_dist_estimation_loss_nj.png"))

    # -------------------------------------------------------------
    # 4. Absolute Performance Heatmap per Pipeline
    # -------------------------------------------------------------
    for pname in pipes:
        grid_abs = get_grid(pname)
        if grid_abs is not None:
            clean_name = pname.replace("+", "_").replace(" ", "_")
            _plot_absolute_heatmap(grid_abs, f"Absolute Accuracy: {pname} (Mean nRF)",
                                  os.path.join(comp_dir, f"absolute_nrf_{clean_name}.png"))

def _plot_absolute_heatmap(grid, title, out_path):
    """Plots absolute nRF heatmap (0.0 to 1.0) using YlOrRd color palette."""
    plt.figure(figsize=(9, 7))
    sns.set_theme(style="whitegrid", font_scale=1.1)

    sns.heatmap(
        grid,
        annot=True,
        fmt=".3f",
        cmap="YlOrRd",
        vmin=0.0,
        vmax=1.0,
        cbar_kws={'label': 'Normalized RF Distance (0 = Perfect, 1 = Completely Disjoint)'}
    )
    plt.title(title, fontsize=13, fontweight="bold", pad=15)
    plt.xlabel("Evolutionary Distance D (substitutions/site)", fontsize=12)
    plt.ylabel("Sequence Length L (amino acids)", fontsize=12)
    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    plt.close()
    print(f"Generated: {out_path}")

def _plot_heatmap(delta_grid, title, cbar_label, out_path):
    plt.figure(figsize=(9, 7))
    sns.set_theme(style="whitegrid", font_scale=1.1)
    cmap = sns.diverging_palette(240, 10, as_cmap=True)

    sns.heatmap(
        delta_grid,
        annot=True,
        fmt=".3f",
        cmap=cmap,
        center=0,
        cbar_kws={'label': cbar_label}
    )
    plt.title(title, fontsize=13, fontweight="bold", pad=15)
    plt.xlabel("Evolutionary Distance D (substitutions/site)", fontsize=12)
    plt.ylabel("Sequence Length L (amino acids)", fontsize=12)
    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    plt.close()
    print(f"Generated: {out_path}")

def generate_boxplots(df, outdir):
    """Generates grid of nRF distribution boxplots with mean annotations across all conditions."""
    distances = sorted(df["distance"].unique())
    lengths = sorted(df["length"].unique())
    
    # Preferred ordering of known pipelines
    preferred_order = [
        "PWA+NJ", "MSA+NJ", "MSA+ML", "MSA+RAXML", "GS", "MSA+BI",
        "TRUE_PWA+NJ", "TRUE_MSA+NJ", "TRUE_MSA+ML", "TRUE_MSA+RAXML", "TRUE_MSA+BI", "TRUE_DIST+NJ",
        "PWA+FastME", "MSA+FastME", "PWA+FastME_SPR", "MSA+FastME_LG_G"
    ]
    present_pipes = list(df["pipeline"].unique())
    pipes = [p for p in preferred_order if p in present_pipes] + [p for p in present_pipes if p not in preferred_order]

    fig, axes = plt.subplots(
        nrows=len(distances),
        ncols=len(lengths),
        figsize=(max(15, len(lengths) * 4), max(12, len(distances) * 3)),
        sharey=True
    )

    if len(distances) == 1 and len(lengths) == 1:
        axes = np.array([[axes]])
    elif len(distances) == 1:
        axes = np.array([axes])
    elif len(lengths) == 1:
        axes = np.array([[ax] for ax in axes])

    # Dynamic palette fallback to avoid any missing keys
    default_colors = sns.color_palette("tab10", len(pipes))
    palette = {p: PIPELINE_COLORS.get(p, default_colors[i % len(default_colors)]) for i, p in enumerate(pipes)}

    for r, d in enumerate(distances):
        for c, l in enumerate(lengths):
            ax = axes[r, c]
            sub = df[(df["distance"] == d) & (df["length"] == l)]

            if not sub.empty:
                sns.boxplot(
                    data=sub,
                    x="pipeline",
                    y="nrf_distance",
                    hue="pipeline",
                    legend=False,
                    order=pipes,
                    palette=palette,
                    ax=ax,
                    width=0.6,
                    fliersize=2
                )
                means = sub.groupby("pipeline")["nrf_distance"].mean()
                for idx, p in enumerate(pipes):
                    if p in means:
                        m_val = means[p]
                        ax.text(idx, min(0.95, m_val + 0.05), f"μ={m_val:.2f}",
                                ha="center", va="bottom", fontsize=8, fontweight="bold", color="#111111")

            ax.set_title(f"D={d}, L={l}", fontsize=10, fontweight="bold")
            ax.set_xlabel("")
            if c == 0:
                ax.set_ylabel("nRF Distance", fontsize=10)
            else:
                ax.set_ylabel("")
            ax.set_ylim(-0.05, 1.05)
            ax.tick_params(axis='x', rotation=30, labelsize=8)

    plt.suptitle("Normalized RF Distance Distributions Across All Conditions", fontsize=15, fontweight="bold", y=0.995)
    plt.tight_layout()
    out_path = os.path.join(outdir, "nrf_boxplots.png")
    plt.savefig(out_path, dpi=300)
    plt.close()
    print(f"Generated: {out_path}")

def generate_summary_table(df, outdir):
    """Generates aggregated summary statistics table."""
    summary = df.groupby(["distance", "length", "pipeline"])["nrf_distance"].agg(
        count="count",
        mean="mean",
        std="std",
        median="median",
        min="min",
        max="max"
    ).reset_index()

    out_csv = os.path.join(outdir, "summary_statistics.csv")
    summary.to_csv(out_csv, index=False)
    print(f"Generated: {out_csv}")


# ==========================================
# 1.5 Sequence Similarity Score (SSS) Benchmark Analysis
# ==========================================

def generate_sss_analysis(df, outdir):
    """
    Performs Sequence Similarity Score (SSS w_bar) benchmark analysis matching Matsui & Iwasaki (2020) Fig. 3a:
    1. nrf_vs_sss_curves.png: LOWESS smoothed curves of topological error (nRF) & accuracy (1-nRF) vs SSS.
    2. sss_distribution_by_distance.png: Distribution of SSS across distances, highlighting extreme breakdown rates.
    3. sss_breakdown_report.csv: Stratified performance metrics across SSS regimes (<=0.03, 0.03-0.06, 0.06-0.15, >=0.15).
    """
    if "sss_mean" not in df.columns or df["sss_mean"].dropna().empty:
        print("Note: sss_mean column not present or empty in dataframe. Skipping SSS analysis.")
        return

    # 1. Fig 3a Replication: nRF & Accuracy vs SSS
    preferred_order = [
        "PWA+NJ", "MSA+NJ", "MSA+ML", "MSA+RAXML", "GS", "MSA+BI",
        "TRUE_PWA+NJ", "TRUE_MSA+NJ", "TRUE_MSA+ML", "TRUE_MSA+RAXML", "TRUE_MSA+BI", "TRUE_DIST+NJ"
    ]
    pipes = [p for p in preferred_order if p in df["pipeline"].unique()] + [p for p in df["pipeline"].unique() if p not in preferred_order]

    fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(16, 7), sharex=True)

    for pipe in pipes:
        sub = df[(df["pipeline"] == pipe) & df["sss_mean"].notna() & df["nrf_distance"].notna()].copy()
        if len(sub) < 3:
            continue
        sub = sub.sort_values("sss_mean")
        x = sub["sss_mean"].values
        y_err = sub["nrf_distance"].values
        y_acc = 1.0 - y_err

        color = PIPELINE_COLORS.get(pipe, "#333333")

        # Plot LOWESS smooth curve if statsmodels is available
        plotted_smooth = False
        try:
            from statsmodels.nonparametric.smoothers_lowess import lowess
            frac_param = 0.4 if len(x) > 20 else 0.6
            sm_err = lowess(y_err, x, frac=frac_param, it=3)
            sm_acc = lowess(y_acc, x, frac=frac_param, it=3)
            ax1.plot(sm_err[:, 0], sm_err[:, 1], color=color, linewidth=2.5, label=pipe)
            ax2.plot(sm_acc[:, 0], sm_acc[:, 1], color=color, linewidth=2.5, label=pipe)
            plotted_smooth = True
        except Exception:
            pass

        if not plotted_smooth:
            # Fallback to rolling mean
            sub["smooth_err"] = sub["nrf_distance"].rolling(window=max(3, len(sub)//10), min_periods=1).mean()
            sub["smooth_acc"] = (1.0 - sub["nrf_distance"]).rolling(window=max(3, len(sub)//10), min_periods=1).mean()
            ax1.plot(sub["sss_mean"], sub["smooth_err"], color=color, linewidth=2.5, label=pipe)
            ax2.plot(sub["sss_mean"], sub["smooth_acc"], color=color, linewidth=2.5, label=pipe)

        # Plot semitransparent scatter points
        ax1.scatter(x, y_err, color=color, alpha=0.10, s=14, edgecolors="none")
        ax2.scatter(x, y_acc, color=color, alpha=0.10, s=14, edgecolors="none")

    # Threshold annotations
    for ax, title, ylabel in [
        (ax1, "Topological Error (nRF) vs Sequence Similarity Score (SSS)", "Normalized RF Distance (Lower is Better)"),
        (ax2, "Topological Accuracy (1 - nRF) vs Sequence Similarity Score (SSS)\n[Matsui & Iwasaki (2020) Fig. 3a Replication]", "Correct Topology Ratio (1 - nRF, Higher is Better)")
    ]:
        ax.axvline(0.06, color="#d62728", linestyle="--", linewidth=1.5, alpha=0.85, label=r"Crossover Threshold ($SSS = 0.06$)")
        ax.axvline(0.03, color="#8b0000", linestyle=":", linewidth=1.8, alpha=0.85, label=r"Extreme Breakdown ($SSS \leq 0.03$)")
        ax.set_title(title, fontsize=12, fontweight="bold", pad=12)
        ax.set_xlabel(r"Average Sequence Similarity Score ($\bar{w}$)", fontsize=11)
        ax.set_ylabel(ylabel, fontsize=11)
        ax.set_xlim(-0.01, min(0.65, max(0.2, df["sss_mean"].max() + 0.02)))
        ax.grid(True, linestyle="--", alpha=0.5)
        ax.legend(fontsize=9, loc="best", framealpha=0.9)

    plt.tight_layout()
    plot_path1 = os.path.join(outdir, "nrf_vs_sss_curves.png")
    plt.savefig(plot_path1, dpi=300)
    plt.close()
    print(f"Generated: {plot_path1}")

    # 2. SSS Distribution across Distances
    plt.figure(figsize=(11, 6))
    rep_level = df.drop_duplicates(subset=["distance", "length", "replicate"]).copy()
    distances = sorted(rep_level["distance"].unique())

    ax = sns.boxplot(data=rep_level, x="distance", y="sss_mean", color="#9ecae1", width=0.5, fliersize=2)
    sns.stripplot(data=rep_level, x="distance", y="sss_mean", color="#3182bd", alpha=0.35, size=4, jitter=0.2, ax=ax)

    ax.axhline(0.06, color="#d62728", linestyle="--", linewidth=1.5, label=r"syz049 Threshold ($SSS = 0.06$)")
    ax.axhline(0.03, color="#8b0000", linestyle=":", linewidth=1.8, label=r"Extreme Breakdown ($SSS \leq 0.03$)")

    # Annotate fraction below 0.06 and 0.03 per distance
    for idx, d in enumerate(distances):
        sub_d = rep_level[rep_level["distance"] == d]["sss_mean"]
        n_tot = len(sub_d)
        if n_tot > 0:
            n_sub06 = np.count_nonzero(sub_d < 0.06)
            n_sub03 = np.count_nonzero(sub_d <= 0.03)
            pct06 = (n_sub06 / n_tot) * 100.0
            pct03 = (n_sub03 / n_tot) * 100.0
            ax.text(idx, -0.015, f"<0.06: {pct06:.1f}%\n≤0.03: {pct03:.1f}%",
                    ha="center", va="top", fontsize=8, color="#8b0000", fontweight="bold")

    ax.set_title("Distribution of Sequence Similarity Score (SSS) across Evolutionary Distances", fontsize=13, fontweight="bold", pad=15)
    ax.set_xlabel("Evolutionary Distance D (substitutions/site)", fontsize=11)
    ax.set_ylabel(r"Average Sequence Similarity Score ($\bar{w}$)", fontsize=11)
    ax.set_ylim(-0.06, max(0.5, rep_level["sss_mean"].max() + 0.05))
    ax.legend(loc="upper right", fontsize=9)
    plt.tight_layout()
    plot_path2 = os.path.join(outdir, "sss_distribution_by_distance.png")
    plt.savefig(plot_path2, dpi=300)
    plt.close()
    print(f"Generated: {plot_path2}")

    # 3. Stratified SSS Regime Report Table
    def categorize_sss(val):
        if pd.isna(val):
            return "Unknown"
        if val <= 0.03:
            return "1. Extreme (SSS <= 0.03)"
        elif val < 0.06:
            return "2. Severe (0.03 < SSS < 0.06)"
        elif val < 0.15:
            return "3. Moderate (0.06 <= SSS < 0.15)"
        else:
            return "4. Mild (SSS >= 0.15)"

    df_strat = df.copy()
    df_strat["sss_regime"] = df_strat["sss_mean"].apply(categorize_sss)
    report_rows = []
    for (regime, pipe), grp in df_strat.groupby(["sss_regime", "pipeline"]):
        report_rows.append({
            "sss_regime": regime,
            "pipeline": pipe,
            "replicate_count": len(grp),
            "mean_nrf": round(float(grp["nrf_distance"].mean()), 4),
            "std_nrf": round(float(grp["nrf_distance"].std()), 4),
            "median_nrf": round(float(grp["nrf_distance"].median()), 4),
            "accuracy_mean": round(1.0 - float(grp["nrf_distance"].mean()), 4),
            "mean_sss": round(float(grp["sss_mean"].mean()), 4),
        })
    if report_rows:
        rep_df = pd.DataFrame(report_rows).sort_values(["sss_regime", "mean_nrf"])
        csv_path = os.path.join(outdir, "sss_breakdown_report.csv")
        rep_df.to_csv(csv_path, index=False)
        print(f"Generated: {csv_path}")


# ==========================================
# 2. Sequence Length & Retention Analysis
# ==========================================

def analyze_sequence_lengths(scan_dir, outdir):
    """
    Scans seqs.fasta across replication directories or work/ and generates length distribution reports.
    """
    from Bio import SeqIO
    records = []

    # 1. Try scanning replication directories D*_L*_rep*
    rep_dirs = glob.glob(os.path.join(scan_dir, "**", "D*_L*_rep*"), recursive=True)
    if rep_dirs:
        print(f"Analyzing sequence lengths from {len(rep_dirs)} replication directories...")
        for cdir in rep_dirs:
            dname = os.path.basename(cdir)
            parts = dname.split("_")
            if len(parts) != 3:
                continue
            try:
                dist = float(parts[0].replace("D", ""))
                length = int(parts[1].replace("L", ""))
                rep = int(parts[2].replace("rep", ""))
            except ValueError:
                continue

            fasta_path = os.path.join(cdir, "seqs.fasta")
            if not os.path.exists(fasta_path):
                continue

            for rec in SeqIO.parse(fasta_path, "fasta"):
                clean_seq = str(rec.seq).replace("-", "").upper()
                seq_len = len(clean_seq)
                records.append({
                    "distance": dist,
                    "initial_length": length,
                    "replicate": rep,
                    "taxon": rec.id,
                    "actual_length": seq_len,
                    "retention_rate": (seq_len / length) * 100.0 if length > 0 else 0.0
                })

    # 2. Fallback: Scan work/ directories using .command.sh
    if not records:
        fasta_files = glob.glob(os.path.join(scan_dir, "**", "seqs.fasta"), recursive=True)
        if fasta_files:
            print(f"Analyzing sequence lengths from {len(fasta_files)} work task fasta files...")
            for fpath in fasta_files:
                task_dir = os.path.dirname(fpath)
                cmd_sh = os.path.join(task_dir, ".command.sh")
                if not os.path.exists(cmd_sh):
                    continue
                try:
                    with open(cmd_sh, "r") as f:
                        content = f.read()
                    dist_m = re.search(r"--distance\s+([\d\.]+)", content)
                    len_m  = re.search(r"--length\s+(\d+)", content)
                    rep_m  = re.search(r"--seed\s+(\d+)", content)
                    if not (dist_m and len_m and rep_m):
                        continue
                    dist = float(dist_m.group(1))
                    length = int(len_m.group(1))
                    rep = int(rep_m.group(1))

                    for rec in SeqIO.parse(fpath, "fasta"):
                        clean_seq = str(rec.seq).replace("-", "").upper()
                        seq_len = len(clean_seq)
                        records.append({
                            "distance": dist,
                            "initial_length": length,
                            "replicate": rep,
                            "taxon": rec.id,
                            "actual_length": seq_len,
                            "retention_rate": (seq_len / length) * 100.0 if length > 0 else 0.0
                        })
                except Exception:
                    continue

    if not records:
        print("Note: No sequence fasta files found for length analysis.")
        return

    len_df = pd.DataFrame(records)
    len_outdir = os.path.join(outdir, "length_analysis")
    os.makedirs(len_outdir, exist_ok=True)

    # Statistics table
    stat_df = len_df.groupby(["initial_length", "distance"])["actual_length"].agg(
        taxa_count="count",
        mean_length="mean",
        std_length="std",
        min_length="min",
        median_length="median",
        max_length="max"
    ).reset_index()
    stat_csv = os.path.join(len_outdir, "sequence_length_statistics.csv")
    stat_df.to_csv(stat_csv, index=False)
    print(f"Generated: {stat_csv}")

    # Condition summary plot
    lengths = sorted(len_df["initial_length"].unique())
    fig, axes = plt.subplots(1, len(lengths), figsize=(max(18, len(lengths) * 4), 5.5), sharey=False)
    if len(lengths) == 1:
        axes = [axes]

    for i, l in enumerate(lengths):
        ax = axes[i]
        sub = len_df[len_df["initial_length"] == l]
        sns.boxplot(data=sub, x="distance", y="actual_length", hue="distance", legend=False, palette="Blues", ax=ax)
        ax.axhline(l, color="red", linestyle="--", linewidth=1.2, label=f"Target L={l}")
        ax.set_title(f"Initial L = {l} aa", fontsize=11, fontweight="bold")
        ax.set_xlabel("Distance D", fontsize=10)
        ax.set_ylabel("Actual Length (aa)" if i == 0 else "")

    plt.suptitle("Remaining Amino Acid Sequence Length vs Evolutionary Distance", fontsize=13, fontweight="bold", y=1.02)
    plt.tight_layout()
    len_plot = os.path.join(len_outdir, "sequence_lengths_all_conditions.png")
    plt.savefig(len_plot, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Generated: {len_plot}")

# ==========================================
# 3. Replication Artifact Harvesting & Visual Reports (1 Replicate per Condition)
# ==========================================

TARGET_REPLICATE_FILES = {
    # Sequences & Alignments
    "seqs.fasta": ["seqs_{rep}.fasta", "seqs.fasta", "unaligned.fasta", "sim_{rep}.unaligned.fa"],
    "true_msa.fasta": ["true_msa_{rep}.fasta", "true_msa.fasta", "sim_{rep}.fa"],
    "msa.fasta": ["msa_{rep}.fasta", "msa.fasta", "mafft_msa.fasta"],

    # Distance Matrices
    "pwa_matrix.phylip": ["pwa_matrix_{rep}.phylip", "pwa_matrix.phylip"],
    "msa_matrix.phylip": ["msa_matrix_{rep}.phylip", "msa_matrix.phylip"],
    "true_matrix.phylip": ["true_matrix_{rep}.phylip", "true_matrix.phylip"],
    "true_pwa_matrix.phylip": ["true_pwa_matrix_{rep}.phylip", "true_pwa_matrix.phylip"],
    "true_msa_matrix.phylip": ["true_msa_matrix_{rep}.phylip", "true_msa_matrix.phylip"],

    # Phylogenetic Trees
    "true_tree.nwk": ["true_tree_{rep}.nwk", "true_tree.nwk", "sim_{rep}.treefile", "sim_{rep}.tree"],
    "pwa_nj.nwk": ["pwa_nj_{rep}.nwk", "pwa_nj.nwk"],
    "msa_nj.nwk": ["msa_nj_{rep}.nwk", "msa_nj.nwk"],
    "msa_ml.nwk": ["msa_ml_{rep}.nwk", "msa_ml.nwk"],
    "msa_bi.nwk": ["msa_bi_{rep}.nwk", "msa_bi.nwk"],
    "gs.nwk": ["gs_{rep}.nwk", "gs.nwk"],
    "true_pwa_nj.nwk": ["true_pwa_nj_{rep}.nwk", "true_pwa_nj.nwk"],
    "true_msa_nj.nwk": ["true_msa_nj_{rep}.nwk", "true_msa_nj.nwk"],
    "true_msa_ml.nwk": ["true_msa_ml_{rep}.nwk", "true_msa_ml.nwk"],
    "true_msa_bi.nwk": ["true_msa_bi_{rep}.nwk", "true_msa_bi.nwk"],
    "true_dist_nj.nwk": ["true_dist_nj_{rep}.nwk", "true_dist_nj.nwk"],

    # Metadata
    "msa_ml_meta.json": ["msa_ml_meta_{rep}.json", "msa_ml_meta.json"],
    "true_msa_ml_meta.json": ["true_msa_ml_meta_{rep}.json", "true_msa_ml_meta.json"]
}


def extract_condition_info(task_dir: str) -> Optional[Dict[str, Any]]:
    """
    Extracts experimental condition parameters (distance D, length L, taxa N, chunk, replicate)
    from a Nextflow task directory.

    Parameters
    ----------
    task_dir : str
        Path to the Nextflow task execution directory inside work/.

    Returns
    -------
    Optional[Dict[str, Any]]
        Dictionary containing extracted condition metadata ('cond_key', 'dist', 'length', 'taxa', 'chunk', 'rep'),
        or None if no condition metadata could be determined.
    """
    # 1. Inspect CSV filenames in task directory
    csv_files = glob.glob(os.path.join(task_dir, "*.csv"))
    for cf in csv_files:
        bname = os.path.basename(cf)
        if bname.startswith(".") or "benchmark" in bname:
            continue
        m_taxa_chk = re.search(r"_N(\d+)_D([\d\.]+)_L(\d+)_chk(\d+)\.csv$", bname)
        if m_taxa_chk:
            taxa, dist, length, chk = m_taxa_chk.groups()
            return {"cond_key": f"N{taxa}_D{dist}_L{length}", "dist": float(dist), "length": int(length), "taxa": int(taxa), "chunk": int(chk), "rep": None}
        m_std_chk = re.search(r"_D([\d\.]+)_L(\d+)_chk(\d+)\.csv$", bname)
        if m_std_chk:
            dist, length, chk = m_std_chk.groups()
            return {"cond_key": f"D{dist}_L{length}", "dist": float(dist), "length": int(length), "taxa": None, "chunk": int(chk), "rep": None}
        m_taxa_rep = re.search(r"_N(\d+)_D([\d\.]+)_L(\d+)_rep(\d+)\.csv$", bname)
        if m_taxa_rep:
            taxa, dist, length, rep = m_taxa_rep.groups()
            return {"cond_key": f"N{taxa}_D{dist}_L{length}", "dist": float(dist), "length": int(length), "taxa": int(taxa), "chunk": None, "rep": int(rep)}
        m_std_rep = re.search(r"_D([\d\.]+)_L(\d+)_rep(\d+)\.csv$", bname)
        if m_std_rep:
            dist, length, rep = m_std_rep.groups()
            return {"cond_key": f"D{dist}_L{length}", "dist": float(dist), "length": int(length), "taxa": None, "chunk": None, "rep": int(rep)}

    # 2. Inspect .command.run
    cmd_run = os.path.join(task_dir, ".command.run")
    if os.path.exists(cmd_run):
        try:
            with open(cmd_run, "r", errors="ignore") as rf:
                txt = rf.read()
            m_tag = re.search(r"(?:NXF_TASK_TAG|# NEXTFLOW TASK:[^(\n]*\()\s*['\"]?(?:N=(\d+)_)?D=([\d\.]+)_L=(\d+)(?:_chk=(\d+))?(?:\[(\d+)\.\.(\d+)\])?(?:_rep=(\d+))?['\"]?", txt)
            if m_tag:
                taxa, dist, length, chk, r_start, r_end, rep = m_tag.groups()
                rep_val = int(rep) if rep else (int(r_start) if r_start else None)
                cond_key = f"N{taxa}_D{dist}_L{length}" if taxa else f"D{dist}_L{length}"
                return {
                    "cond_key": cond_key,
                    "dist": float(dist),
                    "length": int(length),
                    "taxa": int(taxa) if taxa else None,
                    "chunk": int(chk) if chk else None,
                    "rep": rep_val
                }
        except Exception:
            pass

    # 3. Inspect .command.sh
    cmd_sh = os.path.join(task_dir, ".command.sh")
    if os.path.exists(cmd_sh):
        try:
            with open(cmd_sh, "r", errors="ignore") as sf:
                txt = sf.read()
            dist_m = re.search(r"--(?:distance|branch-scale)\s+([\d\.]+)", txt)
            len_m = re.search(r"--length\s+(\d+)", txt)
            taxa_m = re.search(r"(?:--num_taxa\s+(\d+)|/(\d+)\}\")", txt)
            rep_m = re.search(r"--(?:seed|replicate)\s+(\d+)|-seed\s+(\d+)", txt)
            if dist_m and len_m:
                dist = float(dist_m.group(1))
                length = int(len_m.group(1))
                taxa = int(taxa_m.group(1) or taxa_m.group(2)) if taxa_m else None
                rep_val = int(rep_m.group(1) or rep_m.group(2)) if rep_m else None
                cond_key = f"N{taxa}_D{dist}_L{length}" if taxa else f"D{dist}_L{length}"
                return {
                    "cond_key": cond_key,
                    "dist": dist,
                    "length": length,
                    "taxa": taxa,
                    "chunk": None,
                    "rep": rep_val
                }
        except Exception:
            pass

    return None


def scan_work_condition_replicates(
    work_dir: str,
    target_rep: int = 1
) -> Dict[str, Dict[str, Any]]:
    """
    Scans Nextflow work/ directory and indexes all artifact files for target replicate per condition.

    Parameters
    ----------
    work_dir : str
        Path to Nextflow work/ directory.
    target_rep : int, optional
        Target replicate number to harvest for each condition (default: 1).

    Returns
    -------
    Dict[str, Dict[str, Any]]
        Mapping of condition key to dictionary with metadata and matched filepaths.
    """
    print(f"Scanning Nextflow work directory '{work_dir}' for replicate {target_rep} artifacts across conditions...")
    task_dirs = set()
    for marker in [".command.sh", ".command.run"]:
        for p in glob.glob(os.path.join(work_dir, "**", marker), recursive=True):
            task_dirs.add(os.path.dirname(p))

    conditions: Dict[str, Dict[str, Any]] = defaultdict(lambda: {"metadata": {}, "files": {}})

    for t_dir in sorted(task_dirs):
        exit_file = os.path.join(t_dir, ".exitcode")
        if os.path.exists(exit_file):
            try:
                with open(exit_file, "r") as ef:
                    if ef.read().strip() != "0":
                        continue
            except Exception:
                pass

        info = extract_condition_info(t_dir)
        if not info:
            continue

        cond_key = info["cond_key"]
        if not conditions[cond_key]["metadata"]:
            conditions[cond_key]["metadata"] = info

        # Search for target replicate artifacts in this task directory
        for canonical_name, patterns in TARGET_REPLICATE_FILES.items():
            if canonical_name in conditions[cond_key]["files"]:
                continue
            for pat in patterns:
                fname = pat.format(rep=target_rep)
                cand_path = os.path.join(t_dir, fname)
                if os.path.exists(cand_path) and os.path.getsize(cand_path) > 0:
                    conditions[cond_key]["files"][canonical_name] = cand_path
                    break

    print(f"Identified {len(conditions)} unique experimental conditions in '{work_dir}'.")
    return conditions


def compute_tree_nrf(true_tree_file: str, est_tree_file: str) -> Tuple[Optional[float], Optional[float]]:
    """
    Computes Robinson-Foulds (RF) and normalized Robinson-Foulds (nRF) distance using DendroPy.

    Parameters
    ----------
    true_tree_file : str
        Path to ground-truth Newick tree.
    est_tree_file : str
        Path to estimated Newick tree.

    Returns
    -------
    Tuple[Optional[float], Optional[float]]
        (rf_distance, nrf_distance). Returns (None, None) if calculation fails.
    """
    try:
        import dendropy
        from dendropy.calculate import treecompare
        tns = dendropy.TaxonNamespace()
        t_true = dendropy.Tree.get(path=true_tree_file, schema="newick", taxon_namespace=tns, preserve_underscores=True)
        t_est = dendropy.Tree.get(path=est_tree_file, schema="newick", taxon_namespace=tns, preserve_underscores=True)
        t_true.is_rooted = False
        t_true.deroot()
        t_est.is_rooted = False
        t_est.deroot()
        t_true.encode_bipartitions()
        t_est.encode_bipartitions()
        rf = treecompare.symmetric_difference(t_true, t_est)
        num_taxa = len(tns)
        max_rf = 2 * (num_taxa - 3)
        nrf = rf / max_rf if max_rf > 0 else 0.0
        return float(rf), float(nrf)
    except Exception:
        return None, None


def parse_phylip_matrix(matrix_file: str) -> Tuple[List[str], np.ndarray]:
    """
    Parses a PHYLIP format distance matrix.

    Parameters
    ----------
    matrix_file : str
        Path to PHYLIP matrix file.

    Returns
    -------
    Tuple[List[str], np.ndarray]
        Taxa names list and NxN numpy distance matrix.
    """
    with open(matrix_file, "r") as f:
        lines = [l.strip() for l in f if l.strip()]
    if not lines:
        raise ValueError(f"Empty PHYLIP matrix file: {matrix_file}")

    n_taxa = int(lines[0].split()[0])
    names: List[str] = []
    mat = np.zeros((n_taxa, n_taxa), dtype=float)
    for i, line in enumerate(lines[1:n_taxa + 1]):
        tokens = line.split()
        names.append(tokens[0])
        mat[i, :] = [float(x) for x in tokens[1:n_taxa + 1]]
    return names, mat


def compute_patristic_matrix_from_tree(tree_file: str, taxon_order: List[str]) -> Optional[np.ndarray]:
    """
    Computes pairwise patristic distance matrix from a Newick tree for specified taxon order.

    Parameters
    ----------
    tree_file : str
        Path to Newick tree.
    taxon_order : List[str]
        Ordered list of taxon names.

    Returns
    -------
    Optional[np.ndarray]
        NxN patristic distance matrix, or None if calculation fails.
    """
    try:
        import dendropy
        tns = dendropy.TaxonNamespace(taxon_order)
        tree = dendropy.Tree.get(path=tree_file, schema="newick", taxon_namespace=tns, preserve_underscores=True)
        pdm = tree.phylogenetic_distance_matrix()
        N = len(taxon_order)
        mat = np.zeros((N, N), dtype=float)
        for i, t1 in enumerate(taxon_order):
            tx1 = tns.get_taxon(t1)
            for j, t2 in enumerate(taxon_order):
                if i != j:
                    tx2 = tns.get_taxon(t2)
                    mat[i, j] = pdm.distance(tx1, tx2)
        return mat
    except Exception:
        return None


def plot_replicate_trees(
    tree_files: Dict[str, str],
    out_path: str,
    true_tree_file: Optional[str] = None
) -> None:
    """
    Draws a multi-panel visual comparison of phylogenetic trees for a single replicate.

    Parameters
    ----------
    tree_files : Dict[str, str]
        Dictionary of {method_label: nwk_filepath}.
    out_path : str
        Output PNG image path.
    true_tree_file : Optional[str], optional
        Path to the true tree Newick file for nRF comparison.
    """
    try:
        from Bio import Phylo
    except ImportError:
        print("Warning: Bio.Phylo is not available. Skipping tree plot.")
        return

    method_order = [
        ("true_tree.nwk", "True Tree (Reference)"),
        ("pwa_nj.nwk", "PWA + NJ"),
        ("msa_nj.nwk", "MSA + NJ"),
        ("msa_ml.nwk", "MSA + ML"),
        ("msa_bi.nwk", "MSA + BI (MrBayes)"),
        ("gs.nwk", "Graph Splitting (GS)"),
        ("true_pwa_nj.nwk", "True PWA + NJ"),
        ("true_msa_nj.nwk", "True MSA + NJ"),
        ("true_msa_ml.nwk", "True MSA + ML"),
        ("true_msa_bi.nwk", "True MSA + BI"),
        ("true_dist_nj.nwk", "True Dist + NJ")
    ]

    selected = []
    for k, lbl in method_order:
        if k in tree_files and os.path.exists(tree_files[k]):
            selected.append((k, lbl, tree_files[k]))

    for k, p in tree_files.items():
        if k.endswith(".nwk") and k not in [m[0] for m in method_order] and os.path.exists(p):
            selected.append((k, k.replace(".nwk", "").replace("_", " ").upper(), p))

    if not selected:
        return

    n_trees = len(selected)
    ncols = min(4, n_trees)
    nrows = int(np.ceil(n_trees / ncols))

    fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 5.0, nrows * 4.5), squeeze=False)
    ax_flat = axes.flatten()

    ref_tree_path = true_tree_file or tree_files.get("true_tree.nwk")

    for idx, (key, label, path) in enumerate(selected):
        ax = ax_flat[idx]
        try:
            tree = Phylo.read(path, "newick")
            title = label
            if key != "true_tree.nwk" and ref_tree_path and os.path.exists(ref_tree_path):
                _, nrf = compute_tree_nrf(ref_tree_path, path)
                if nrf is not None:
                    title += f"\n(nRF = {nrf:.3f})"
            Phylo.draw(tree, axes=ax, do_show=False)
            ax.set_title(title, fontsize=11, fontweight="bold", pad=10)
            ax.set_xlabel("")
            ax.set_ylabel("")
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.spines["left"].set_visible(False)
            ax.spines["bottom"].set_visible(False)
            ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
        except Exception as e:
            ax.text(0.5, 0.5, f"Could not render:\n{label}\n({e})", ha="center", va="center", fontsize=9)
            ax.set_title(label, fontsize=10, fontweight="bold")
            ax.axis("off")

    for idx in range(n_trees, len(ax_flat)):
        ax_flat[idx].axis("off")

    plt.suptitle("Phylogenetic Trees Comparison (1 Replicate Sample)", fontsize=14, fontweight="bold", y=0.99)
    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()


def plot_replicate_distance_matrices(
    matrix_files: Dict[str, str],
    out_path: str,
    true_tree_file: Optional[str] = None
) -> None:
    """
    Plots side-by-side heatmaps of pairwise distance matrices for a single replicate.

    Parameters
    ----------
    matrix_files : Dict[str, str]
        Dictionary of {matrix_name: phylip_filepath}.
    out_path : str
        Output PNG image path.
    true_tree_file : Optional[str], optional
        Optional path to true tree to compute true patristic distance matrix if not in files.
    """
    parsed: Dict[str, Tuple[List[str], np.ndarray]] = {}

    for k, p in matrix_files.items():
        if k.endswith(".phylip") and os.path.exists(p):
            try:
                names, mat = parse_phylip_matrix(p)
                label = k.replace("_matrix.phylip", "").replace(".phylip", "").upper()
                parsed[label] = (names, mat)
            except Exception:
                pass

    # Compute true patristic distance matrix if true_tree is provided
    if "TRUE" not in parsed and true_tree_file and os.path.exists(true_tree_file):
        base_names = next(iter(parsed.values()))[0] if parsed else None
        if base_names:
            true_mat = compute_patristic_matrix_from_tree(true_tree_file, base_names)
            if true_mat is not None:
                parsed["TRUE PATRISTIC"] = (base_names, true_mat)

    if not parsed:
        return

    # Add difference panel if PWA and MSA are both present
    diff_mat = None
    diff_taxa = None
    if "PWA" in parsed and "MSA" in parsed:
        pwa_names, pwa_mat = parsed["PWA"]
        msa_names, msa_mat = parsed["MSA"]
        if pwa_names == msa_names:
            diff_mat = pwa_mat - msa_mat
            diff_taxa = pwa_names

    panels = list(parsed.items())
    if diff_mat is not None and diff_taxa is not None:
        panels.append(("DIFFERENCE (PWA - MSA)", (diff_taxa, diff_mat)))

    n_panels = len(panels)
    fig, axes = plt.subplots(1, n_panels, figsize=(5.2 * n_panels, 4.8), squeeze=False)
    ax_flat = axes.flatten()

    # Calculate common max for distance matrices
    max_d = max(np.max(mat) for lbl, (_, mat) in panels if "DIFFERENCE" not in lbl)
    max_d = max(0.1, max_d)

    for idx, (label, (taxa, mat)) in enumerate(panels):
        ax = ax_flat[idx]
        if "DIFFERENCE" in label:
            cmap = sns.diverging_palette(240, 10, as_cmap=True)
            vlim = max(0.01, np.max(np.abs(mat)))
            sns.heatmap(mat, cmap=cmap, vmin=-vlim, vmax=vlim, center=0,
                        xticklabels=taxa, yticklabels=taxa, ax=ax,
                        cbar_kws={'label': 'Δ Distance (PWA - MSA)'})
        else:
            sns.heatmap(mat, cmap="viridis", vmin=0.0, vmax=max_d,
                        xticklabels=taxa, yticklabels=taxa, ax=ax,
                        cbar_kws={'label': 'Substitutions / site'})
        ax.set_title(f"{label} Matrix", fontsize=11, fontweight="bold", pad=10)
        ax.tick_params(axis='both', labelsize=8)

    plt.suptitle("Pairwise Distance Matrix Heatmaps (1 Replicate Sample)", fontsize=13, fontweight="bold", y=1.02)
    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()


def plot_replicate_msa_overview(
    fasta_files: Dict[str, str],
    out_path: str
) -> None:
    """
    Plots sequence length distribution and alignment gap density across positions.

    Parameters
    ----------
    fasta_files : Dict[str, str]
        Dictionary containing paths to 'seqs.fasta', 'msa.fasta', and optionally 'true_msa.fasta'.
    out_path : str
        Output PNG image path.
    """
    try:
        from Bio import SeqIO
    except ImportError:
        print("Warning: Bio.SeqIO is not available. Skipping MSA overview plot.")
        return

    seqs_path = fasta_files.get("seqs.fasta")
    msa_path = fasta_files.get("msa.fasta")
    true_msa_path = fasta_files.get("true_msa.fasta")

    if not seqs_path and not msa_path:
        return

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))

    # 1. Sequence Length per Taxon
    if seqs_path and os.path.exists(seqs_path):
        seq_records = list(SeqIO.parse(seqs_path, "fasta"))
        taxa = [r.id for r in seq_records]
        lengths = [len(str(r.seq).replace("-", "")) for r in seq_records]

        y_pos = np.arange(len(taxa))
        ax1.barh(y_pos, lengths, color="#1f77b4", alpha=0.8, edgecolor="#0d47a1")
        ax1.set_yticks(y_pos)
        ax1.set_yticklabels(taxa, fontsize=8)
        ax1.set_xlabel("Remaining Sequence Length (aa)", fontsize=10, fontweight="bold")
        ax1.set_title("Unaligned Sequence Length per Taxon", fontsize=11, fontweight="bold")
        ax1.grid(axis="x", linestyle="--", alpha=0.7)
    else:
        ax1.axis("off")

    # 2. Alignment Gap Fraction Profile
    has_msa = msa_path and os.path.exists(msa_path)
    has_true_msa = true_msa_path and os.path.exists(true_msa_path)

    if has_msa or has_true_msa:
        if has_msa:
            msa_recs = list(SeqIO.parse(msa_path, "fasta"))
            if msa_recs:
                n_tax = len(msa_recs)
                aln_len = len(msa_recs[0].seq)
                gap_counts = np.zeros(aln_len)
                for r in msa_recs:
                    s = str(r.seq)
                    for col_idx in range(min(aln_len, len(s))):
                        if s[col_idx] == "-":
                            gap_counts[col_idx] += 1
                gap_pct = (gap_counts / n_tax) * 100.0
                ax2.plot(range(1, aln_len + 1), gap_pct, label=f"MAFFT Inferred MSA (len={aln_len})", color="#ff7f0e", lw=1.2)

        if has_true_msa:
            true_recs = list(SeqIO.parse(true_msa_path, "fasta"))
            if true_recs:
                n_tax = len(true_recs)
                aln_len = len(true_recs[0].seq)
                gap_counts = np.zeros(aln_len)
                for r in true_recs:
                    s = str(r.seq)
                    for col_idx in range(min(aln_len, len(s))):
                        if s[col_idx] == "-":
                            gap_counts[col_idx] += 1
                gap_pct = (gap_counts / n_tax) * 100.0
                ax2.plot(range(1, aln_len + 1), gap_pct, label=f"True MSA (len={aln_len})", color="#2ca02c", lw=1.2, linestyle="--")

        ax2.set_xlabel("Alignment Position (site index)", fontsize=10, fontweight="bold")
        ax2.set_ylabel("Gap Proportion (%)", fontsize=10, fontweight="bold")
        ax2.set_ylim(-2, 102)
        ax2.set_title("Alignment Gap Profile Across Sites", fontsize=11, fontweight="bold")
        ax2.legend(loc="upper right", fontsize=9)
        ax2.grid(True, linestyle="--", alpha=0.7)
    else:
        ax2.axis("off")

    plt.suptitle("Multiple Sequence Alignment & Indel Overview", fontsize=13, fontweight="bold", y=1.02)
    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()


def create_condition_replications_from_work(
    work_dir: str,
    rep_outdir: str,
    target_rep: int = 1,
    generate_plots: bool = True,
    mode: str = "copy",
    threads: int = 16
) -> None:
    """
    Extracts 1 replicate (MSA, distance matrices, trees) per condition from Nextflow work/
    directory, saves them into structured folders under results/replications/, and generates
    visual comparison reports.

    Parameters
    ----------
    work_dir : str
        Path to Nextflow work/ directory.
    rep_outdir : str
        Output destination directory (e.g., results/replications/).
    target_rep : int, optional
        Replicate ID to extract (default: 1).
    generate_plots : bool, optional
        Whether to generate visual comparison figures (default: True).
    mode : str, optional
        File transfer mode: 'copy' or 'hardlink' (default: 'copy').
    threads : int, optional
        Number of worker threads (default: 16).
    """
    cond_data = scan_work_condition_replicates(work_dir, target_rep=target_rep)
    if not cond_data:
        print(f"Warning: No valid condition artifacts found in '{work_dir}'.", file=sys.stderr)
        return

    os.makedirs(rep_outdir, exist_ok=True)
    summary_records = []

    print(f"\nProcessing {len(cond_data)} conditions for replicate {target_rep} into '{rep_outdir}'...")

    for cond_key, item in sorted(cond_data.items()):
        meta = item["metadata"]
        files = item["files"]
        if not files:
            continue

        cond_folder = f"{cond_key}_rep{target_rep}"
        dest_dir = os.path.join(rep_outdir, cond_folder)
        os.makedirs(dest_dir, exist_ok=True)

        copied_paths: Dict[str, str] = {}
        for canonical_name, src_file in files.items():
            dst_file = os.path.join(dest_dir, canonical_name)
            if os.path.exists(dst_file):
                try:
                    os.remove(dst_file)
                except OSError:
                    pass
            if mode == "hardlink":
                try:
                    os.link(src_file, dst_file)
                    copied_paths[canonical_name] = dst_file
                    continue
                except OSError:
                    pass
            shutil.copy2(src_file, dst_file)
            copied_paths[canonical_name] = dst_file

        # Evaluate nRF distances for summary
        true_tree_file = copied_paths.get("true_tree.nwk")
        nrf_results: Dict[str, float] = {}
        if true_tree_file and os.path.exists(true_tree_file):
            for fname, fpath in copied_paths.items():
                if fname.endswith(".nwk") and fname != "true_tree.nwk":
                    pipe_name = fname.replace(".nwk", "").upper()
                    _, nrf = compute_tree_nrf(true_tree_file, fpath)
                    if nrf is not None:
                        nrf_results[pipe_name] = nrf

        # Visualizations
        if generate_plots:
            # 1. Trees comparison
            tree_files = {k: v for k, v in copied_paths.items() if k.endswith(".nwk")}
            if tree_files:
                tree_plot_path = os.path.join(dest_dir, "trees_comparison.png")
                plot_replicate_trees(tree_files, tree_plot_path, true_tree_file=true_tree_file)

            # 2. Distance matrices
            mat_files = {k: v for k, v in copied_paths.items() if k.endswith(".phylip")}
            if mat_files:
                mat_plot_path = os.path.join(dest_dir, "distance_matrices.png")
                plot_replicate_distance_matrices(mat_files, mat_plot_path, true_tree_file=true_tree_file)

            # 3. MSA overview
            fasta_files = {k: v for k, v in copied_paths.items() if k.endswith(".fasta")}
            if fasta_files:
                msa_plot_path = os.path.join(dest_dir, "msa_overview.png")
                plot_replicate_msa_overview(fasta_files, msa_plot_path)

        # Write replicate summary JSON
        summary_dict = {
            "condition": cond_key,
            "replicate": target_rep,
            "distance": meta.get("dist"),
            "initial_length": meta.get("length"),
            "taxa": meta.get("taxa"),
            "nrf_distances": nrf_results,
            "harvested_files": list(copied_paths.keys())
        }
        with open(os.path.join(dest_dir, "summary.json"), "w") as jf:
            json.dump(summary_dict, jf, indent=2)

        record = {
            "condition": cond_key,
            "distance": meta.get("dist"),
            "length": meta.get("length"),
            "taxa": meta.get("taxa"),
            "replicate": target_rep,
            **{f"nrf_{k.lower()}": v for k, v in nrf_results.items()}
        }
        summary_records.append(record)

    if summary_records:
        summary_df = pd.DataFrame(summary_records)
        summary_csv = os.path.join(rep_outdir, "replications_summary.csv")
        summary_df.to_csv(summary_csv, index=False)
        print(f"Generated replications summary table: {summary_csv}")

    print(f"Successfully harvested and visualized 1 replicate per condition in '{rep_outdir}'!")


def organize_replications_from_work(work_dir: str, rep_outdir: str, mode: str = "copy", threads: int = 16) -> None:
    """
    Harvests all tree files, FASTA, MSAs, and matrices across all replicates from work/ into results/replications/
    """
    try:
        from bin.organize_replications import organize_replications
        organize_replications(work_dir, rep_outdir, mode=mode, threads=threads)
    except ImportError:
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        try:
            from organize_replications import organize_replications
            organize_replications(work_dir, rep_outdir, mode=mode, threads=threads)
        except ImportError:
            print("Notice: Fallback to condition-wise single replicate organizer.")
            create_condition_replications_from_work(work_dir, rep_outdir, target_rep=1, mode=mode, threads=threads)


# ==========================================
# Main CLI Entry Point
# ==========================================

def main():
    parser = argparse.ArgumentParser(description="Integrated Benchmark Reporting, Sequence Length Analysis, and Replication Organizer")
    parser.add_argument("--csv", help="Path to aggregated benchmark_summary.csv")
    parser.add_argument("--workdir", help="Path to Nextflow work/ directory (optional: enables auto-aggregation & replication harvesting)")
    parser.add_argument("--repdir", help="Path to results/replications/ directory (optional: for sequence length analysis)")
    parser.add_argument("--outdir", default="results", help="Output directory for reports and figures (default: results)")
    parser.add_argument("--organize_replications", action="store_true", help="Harvest all replication files from work/ into results/replications/")
    parser.add_argument("--create_replications", action="store_true", default=None,
                        help="Extract 1 replicate per condition (MSA, distance matrices, trees) and generate visual reports into outdir/replications/")
    parser.add_argument("--target_rep", type=int, default=1, help="Replicate ID to extract per condition (default: 1)")
    parser.add_argument("--no_rep_plots", action="store_true", help="Skip generating visual comparison plots in replication folders")
    parser.add_argument("--mode", choices=["copy", "hardlink"], default="copy", help="Transfer mode for replication files (default: copy)")
    parser.add_argument("--threads", type=int, default=16, help="Number of worker threads (default: 16)")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    rep_dst = args.repdir or os.path.join(args.outdir, "replications")

    # 1. Phylogenetic Benchmark Analysis
    csv_file = args.csv
    if not csv_file and args.workdir:
        cand = os.path.join(args.outdir, "benchmark_summary.csv")
        if os.path.exists(cand):
            csv_file = cand

    if csv_file and os.path.exists(csv_file) and os.path.getsize(csv_file) > 0:
        df = pd.read_csv(csv_file)
        # Attempt to auto-merge SSS summary if sss_mean is not in df
        if "sss_mean" not in df.columns:
            sss_cand = os.path.join(args.outdir, "sss_summary.csv")
            if not os.path.exists(sss_cand):
                sss_cand = os.path.join(os.path.dirname(csv_file), "sss_summary.csv")
            if os.path.exists(sss_cand):
                df_sss = pd.read_csv(sss_cand).drop_duplicates(subset=["distance", "length", "replicate"])
                df = pd.merge(df, df_sss, on=["distance", "length", "replicate"], how="left")
                print(f"Merged SSS data from '{sss_cand}' into benchmark dataframe.")

        if not df.empty and "nrf_distance" in df.columns:
            print(f"Generating phylogenetic benchmark reports from '{csv_file}' ({len(df)} records)...")
            generate_regime_map(df, args.outdir)
            generate_method_comparisons(df, args.outdir)
            generate_boxplots(df, args.outdir)
            generate_summary_table(df, args.outdir)
            if "sss_mean" in df.columns and df["sss_mean"].notna().any():
                generate_sss_analysis(df, args.outdir)

    # 2. Sequence Length Analysis
    scan_source = args.repdir or (rep_dst if os.path.exists(rep_dst) else args.workdir)
    if scan_source and os.path.exists(scan_source):
        analyze_sequence_lengths(scan_source, args.outdir)

    # 3. Condition-wise 1 Replicate Harvesting & Visual Report Generation
    # Triggered if:
    #   - args.create_replications is explicitly set to True, OR
    #   - args.workdir is provided and exists (and args.create_replications is not False and not args.organize_replications)
    should_create_reps = (
        args.create_replications is True or
        (args.workdir and os.path.exists(args.workdir) and args.create_replications is not False and not args.organize_replications)
    )

    if should_create_reps and args.workdir and os.path.exists(args.workdir):
        create_condition_replications_from_work(
            work_dir=args.workdir,
            rep_outdir=rep_dst,
            target_rep=args.target_rep,
            generate_plots=not args.no_rep_plots,
            mode=args.mode,
            threads=args.threads
        )
    elif args.organize_replications and args.workdir and os.path.exists(args.workdir):
        organize_replications_from_work(args.workdir, rep_dst, mode=args.mode, threads=args.threads)

    print(f"\nAll post-processing reports and analyses completed in '{args.outdir}'!")


if __name__ == "__main__":
    main()

