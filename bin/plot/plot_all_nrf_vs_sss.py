#!/usr/bin/env python3
"""
Publication-Quality SSS Plotting Script: plot_all_nrf_vs_sss.py
Generates comprehensive phylogenetic method benchmark plots with SSS on the horizontal axis
across all 14 result directories, directly replicating Matsui & Iwasaki (2020) Fig. 3a.

For each directory under results/:
1. nrf_vs_sss_curves.png: 2-panel comparison (Log-scale SSS: nRF vs SSS, Accuracy vs SSS) with scatter & LOWESS fits.
2. nrf_vs_sss_linear.png: 2-panel comparison (Linear-scale SSS).
3. nrf_vs_sss_by_condition.png: Stratified multi-panel plots (by Alpha, Taxa, ICS, or Length).
4. sss_breakdown_report.csv: Aggregated performance metrics across SSS regimes.
"""

import os
import sys
import argparse
from pathlib import Path
from typing import Dict, List, Tuple, Any, Optional
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from statsmodels.nonparametric.smoothers_lowess import lowess

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
RESULTS_DIR = PROJECT_ROOT / "results"

# Consistent styling matching Matsui & Iwasaki (2020) and project visual guidelines
plt.rcParams.update({
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "font.family": "sans-serif",
    "axes.edgecolor": "#333333",
    "axes.linewidth": 1.0,
    "grid.color": "#e0e0e0",
    "grid.linestyle": "--",
    "grid.alpha": 0.6,
    "savefig.dpi": 300,
    "figure.dpi": 300,
})

PIPELINE_CONFIGS: Dict[str, Dict[str, Any]] = {
    "MSA+ML": {"color": "#2ca02c", "label": "MSA+ML", "linestyle": "-", "linewidth": 2.5, "zorder": 10},
    "MSA+RAXML": {"color": "#006400", "label": "MSA+RAXML", "linestyle": "-", "linewidth": 2.2, "zorder": 9},
    "MSA+NJ": {"color": "#ff7f0e", "label": "MSA+NJ", "linestyle": "-", "linewidth": 2.4, "zorder": 8},
    "PWA+NJ": {"color": "#1f77b4", "label": "PSA+NJ", "linestyle": "-", "linewidth": 2.4, "zorder": 7},
    "GS": {"color": "#9467bd", "label": "GS", "linestyle": "-", "linewidth": 2.6, "zorder": 11},
    "PWA+FastME": {"color": "#e377c2", "label": "PSA+FastME", "linestyle": "-", "linewidth": 2.2, "zorder": 8},
    "MSA+FastME": {"color": "#bcbd22", "label": "MSA+FastME", "linestyle": "-", "linewidth": 2.2, "zorder": 8},
    "PWA+FastME_SPR": {"color": "#e377c2", "label": "PSA+FastME_SPR", "linestyle": "-", "linewidth": 2.4, "zorder": 9},
    "MSA+FastME_LG_G": {"color": "#bcbd22", "label": "MSA+FastME_LG_G", "linestyle": "-", "linewidth": 2.4, "zorder": 9},
    "TRUE_MSA+ML": {"color": "#8c564b", "label": "TRUE_MSA+ML", "linestyle": "--", "linewidth": 1.8, "zorder": 4},
    "TRUE_MSA+RAXML": {"color": "#c49c94", "label": "TRUE_MSA+RAXML", "linestyle": "-.", "linewidth": 1.8, "zorder": 4},
    "TRUE_MSA+NJ": {"color": "#7f7f7f", "label": "TRUE_MSA+NJ", "linestyle": ":", "linewidth": 1.8, "zorder": 3},
    "TRUE_PWA+NJ": {"color": "#17becf", "label": "TRUE_PSA+NJ", "linestyle": "-.", "linewidth": 2.0, "zorder": 5},
    "TRUE_DIST+NJ": {"color": "#333333", "label": "TRUE_DIST+NJ", "linestyle": ":", "linewidth": 1.5, "zorder": 2},
}

EXPERIMENT_FILES = [
    ("results_alpha", "benchmark_alpha_summary.csv", "alpha"),
    ("results_fastme_NoOption", "benchmark_fastme_summary.csv", "length"),
    ("results_fastme_options", "benchmark_fastme_options_summary.csv", "length"),
    ("results_gamma", "benchmark_gamma_summary.csv", "alpha"),
    ("results_high_dist", "benchmark_high_dist_summary.csv", "length"),
    ("results_high_gap", "benchmark_summary.csv", "length"),
    ("results_ics", "benchmark_ics_summary.csv", "ics_prop"),
    ("results_ics_full", "benchmark_ics_full_summary.csv", "length"),
    ("results_new", "benchmark_summary.csv", "length"),
    ("results_paper_tree", "benchmark_summary.csv", "length"),
    ("results_paper_tree_pre_fix", "benchmark_summary.csv", "length"),
    ("results_zipfian_indel", "benchmark_summary.csv", "length"),
    ("results_simple", "benchmark_high_dist_summary.csv", "length"),
    ("results_taxon", "benchmark_taxon_summary.csv", "num_taxa"),
    ("results_true_pwa", "benchmark_true_pwa_summary.csv", "length"),
]


def compute_lowess_curve(x: np.ndarray, y: np.ndarray, frac: float = 0.3, log_x: bool = True) -> Tuple[np.ndarray, np.ndarray]:
    """Fits LOWESS smoother over x and y, properly handling log scale and MAD=0 robust collapse edge cases."""
    valid = np.isfinite(x) & np.isfinite(y) & (x > 0 if log_x else True)
    x_v, y_v = x[valid], y[valid]
    if len(x_v) < 5:
        return np.array([]), np.array([])
    sorted_idx = np.argsort(x_v)
    x_s, y_s = x_v[sorted_idx], y_v[sorted_idx]
    if len(x_s) < 10:
        return x_s, y_s

    # Perform smoothing in log10 space for x if log scale is requested
    x_in = np.log10(x_s) if log_x else x_s

    # Attempt standard robust lowess (it=3). If residuals collapse (MAD=0 e.g. many y=0)
    # causing lowess to output raw jagged points, fallback to it=0 (local linear regression).
    lw = lowess(y_s, x_in, frac=frac, it=3, return_sorted=True)
    if np.allclose(lw[:, 1], y_s) or np.max(np.abs(np.diff(lw[:, 1]))) > 0.08:
        lw = lowess(y_s, x_in, frac=frac, it=0, return_sorted=True)

    out_x = np.power(10.0, lw[:, 0]) if log_x else lw[:, 0]
    out_y = lw[:, 1]
    return out_x, out_y


def plot_main_nrf_vs_sss(df: pd.DataFrame, out_path: Path, exp_title: str, log_x: bool = True):
    """Generates the main 2-panel figure: nRF vs SSS and Accuracy (1 - nRF) vs SSS."""
    fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(15.5, 6.2))

    pipelines = df["pipeline"].unique()
    for pipe in pipelines:
        sub = df[df["pipeline"] == pipe].dropna(subset=["sss_mean", "nrf_distance"])
        if len(sub) == 0:
            continue

        cfg = PIPELINE_CONFIGS.get(pipe, {
            "color": "#888888", "label": pipe.replace("PWA", "PSA"), "linestyle": "-", "linewidth": 2.0, "zorder": 5
        })

        x = sub["sss_mean"].values
        y_nrf = sub["nrf_distance"].values
        y_acc = 1.0 - y_nrf

        # Scatter
        is_ctrl = "TRUE" in pipe or "Control" in cfg["label"]
        alpha_val = 0.08 if is_ctrl else 0.18
        ax1.scatter(x, y_nrf, color=cfg["color"], alpha=alpha_val, s=12, edgecolors="none", zorder=cfg["zorder"] - 1)
        ax2.scatter(x, y_acc, color=cfg["color"], alpha=alpha_val, s=12, edgecolors="none", zorder=cfg["zorder"] - 1)

        # LOWESS fit
        lx_nrf, ly_nrf = compute_lowess_curve(x, y_nrf, frac=0.3, log_x=log_x)
        lx_acc, ly_acc = compute_lowess_curve(x, y_acc, frac=0.3, log_x=log_x)

        if len(lx_nrf) > 0:
            ax1.plot(lx_nrf, ly_nrf, color=cfg["color"], linestyle=cfg["linestyle"], linewidth=cfg["linewidth"],
                     label=cfg["label"], zorder=cfg["zorder"])
        if len(lx_acc) > 0:
            ax2.plot(lx_acc, ly_acc, color=cfg["color"], linestyle=cfg["linestyle"], linewidth=cfg["linewidth"],
                     label=cfg["label"], zorder=cfg["zorder"])

    min_x = max(0.003, df["sss_mean"].min() * 0.8) if df["sss_mean"].min() > 0 else 0.005
    max_x = min(1.02, max(0.2, df["sss_mean"].max() * 1.05))

    for ax in [ax1, ax2]:
        if log_x:
            ax.set_xscale("log")
            ax.set_xlim(min_x, max_x)
            ax.set_xlabel(r"Sequence Similarity Score: SSS (log scale)", fontweight="bold", fontsize=11)
        else:
            ax.set_xlim(-0.01, max_x)
            ax.set_xlabel(r"Sequence Similarity Score: SSS (linear scale)", fontweight="bold", fontsize=11)

    scale_str = "Log Scale" if log_x else "Linear Scale"
    ax1.set_title(f"A. Topological Error: nRF vs SSS\n{exp_title} - {scale_str}", fontweight="bold", pad=10)
    ax1.set_ylabel(r"Normalized Robinson-Foulds Distance: nRF", fontweight="bold", fontsize=11)
    ax1.set_ylim(-0.02, 1.0)
    ax1.legend(loc="upper right", frameon=True, framealpha=0.9, fontsize=8.5)

    ax2.set_title(f"B. Topological Accuracy: 1 - nRF vs SSS\nMatsui & Iwasaki (2020) Fig. 3a Replication - {scale_str}", fontweight="bold", pad=10)
    ax2.set_ylabel(r"Topological Accuracy: 1 - nRF", fontweight="bold", fontsize=11)
    ax2.set_ylim(-0.02, 1.02)
    ax2.legend(loc="lower right", frameon=True, framealpha=0.9, fontsize=8.5)

    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    plt.close()


def plot_condition_stratified(df: pd.DataFrame, strat_col: str, out_path: Path, exp_title: str):
    """Generates multi-panel figure stratified by experimental parameter (e.g. Length, Alpha, Taxa, ICS)."""
    if strat_col not in df.columns or df[strat_col].nunique() <= 1:
        return

    vals = sorted(df[strat_col].unique())
    ncols = len(vals)
    if ncols > 5:
        nrows = 2
        ncols = (len(vals) + 1) // 2
    else:
        nrows = 1

    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(4.8 * ncols, 4.5 * nrows), sharey=True)
    axes_flat = axes.flatten() if hasattr(axes, "flatten") else [axes]

    strat_label_map = {
        "length": "Sequence Length L = {} aa",
        "alpha": r"Gamma Shape $\alpha = {}$",
        "num_taxa": "Taxon Count N = {}",
        "ics_prop": "ICS Proportion = {}",
    }
    label_fmt = strat_label_map.get(strat_col, f"{strat_col} = {{}}")

    for idx, v in enumerate(vals):
        ax = axes_flat[idx]
        sub_v = df[df[strat_col] == v]

        for pipe in sub_v["pipeline"].unique():
            sub_p = sub_v[sub_v["pipeline"] == pipe].dropna(subset=["sss_mean", "nrf_distance"])
            if len(sub_p) == 0:
                continue
            cfg = PIPELINE_CONFIGS.get(pipe, {
                "color": "#888888", "label": pipe.replace("PWA", "PSA"), "linestyle": "-", "linewidth": 1.8, "zorder": 5
            })

            x = sub_p["sss_mean"].values
            y = sub_p["nrf_distance"].values
            ax.scatter(x, y, color=cfg["color"], alpha=0.12, s=10, edgecolors="none")

            lx, ly = compute_lowess_curve(x, y, frac=0.35)
            if len(lx) > 0:
                ax.plot(lx, ly, color=cfg["color"], linestyle=cfg["linestyle"], linewidth=cfg["linewidth"],
                        label=cfg["label"], zorder=cfg["zorder"])

        ax.set_title(label_fmt.format(v), fontweight="bold", fontsize=11)
        ax.set_xlabel(r"$SSS$", fontweight="bold")
        ax.set_ylim(-0.02, 0.95)
        ax.set_xscale("log")
        if idx % ncols == 0:
            ax.set_ylabel(r"Topological Error: nRF", fontweight="bold")
        if idx == ncols - 1:
            ax.legend(loc="upper right", fontsize=8, frameon=True, framealpha=0.88)

    # Hide extra subplots if any
    for j in range(idx + 1, len(axes_flat)):
        axes_flat[j].set_visible(False)

    plt.suptitle(f"{exp_title} - Stratified by {strat_col.capitalize()} - nRF vs SSS", fontweight="bold", fontsize=13, y=1.02)
    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()


def generate_sss_breakdown_table(df: pd.DataFrame, out_path: Path):
    """Categorizes SSS into 4 standard regimes and generates summary table."""
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

    df_strat = df.dropna(subset=["sss_mean", "nrf_distance"]).copy()
    df_strat["sss_regime"] = df_strat["sss_mean"].apply(categorize_sss)

    rows = []
    for (regime, pipe), grp in df_strat.groupby(["sss_regime", "pipeline"]):
        nrf_m = float(grp["nrf_distance"].mean())
        rows.append({
            "sss_regime": regime,
            "pipeline": pipe,
            "replicate_count": len(grp),
            "mean_nrf": round(nrf_m, 4),
            "std_nrf": round(float(grp["nrf_distance"].std()), 4) if len(grp) > 1 else 0.0,
            "median_nrf": round(float(grp["nrf_distance"].median()), 4),
            "accuracy_mean": round(1.0 - nrf_m, 4),
            "mean_sss": round(float(grp["sss_mean"].mean()), 4),
        })

    if rows:
        rep_df = pd.DataFrame(rows).sort_values(["sss_regime", "mean_nrf"])
        rep_df.to_csv(out_path, index=False)


def process_all_plots():
    print("=" * 70)
    print("Generating SSS-based Publication Figures across all 14 Experiments...")
    print("=" * 70)

    for exp_dir_name, csv_name, strat_col in EXPERIMENT_FILES:
        exp_dir = RESULTS_DIR / exp_dir_name
        csv_path = exp_dir / csv_name
        if not csv_path.exists():
            print(f"[Warning] {csv_path} not found. Skipping.")
            continue

        df = pd.read_csv(csv_path)
        if "sss_mean" not in df.columns:
            print(f"[{exp_dir_name}] 'sss_mean' column missing in {csv_name}. Skipping plot generation.")
            continue

        # Exclude sequence length 100 as requested
        if "length" in df.columns:
            n_before = len(df)
            df = df[pd.to_numeric(df["length"], errors="coerce") != 100].copy()
            if len(df) < n_before:
                print(f"  [Filtered] Excluded length == 100 ({n_before:,} -> {len(df):,} rows)")

        exp_title = exp_dir_name.replace("results_", "").replace("_", " ").title()
        print(f"Processing plots for: {exp_dir_name} ({len(df):,} rows)...")

        out_log = exp_dir / "nrf_vs_sss_curves.png"
        out_lin = exp_dir / "nrf_vs_sss_linear.png"
        out_strat = exp_dir / "nrf_vs_sss_by_condition.png"
        out_rep = exp_dir / "sss_breakdown_report.csv"

        # 1. Main comparison (Log Scale)
        plot_main_nrf_vs_sss(df, out_log, exp_title, log_x=True)
        # 2. Main comparison (Linear Scale)
        plot_main_nrf_vs_sss(df, out_lin, exp_title, log_x=False)
        # 3. Stratified by experimental condition
        plot_condition_stratified(df, strat_col, out_strat, exp_title)
        # 4. SSS breakdown report
        generate_sss_breakdown_table(df, out_rep)

        print(f"  -> Generated: {out_log.name}, {out_lin.name}, {out_strat.name}, {out_rep.name}")

    print("=" * 70)
    print("All SSS figures and breakdown reports generated successfully!")
    print("=" * 70)


if __name__ == "__main__":
    process_all_plots()
