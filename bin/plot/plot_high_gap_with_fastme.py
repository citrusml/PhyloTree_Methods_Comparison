#!/usr/bin/env python3
"""
plot_high_gap_with_fastme.py

Integrates FastME results (from results_fastme_options and results_fastme_NoOption)
into results_high_gap and generates high-resolution comparison plots:
1. nrf_vs_sss_by_condition.png (Stratified by Length: 300, 500, 1000, 1500 aa)
2. nrf_vs_sss_curves.png (Overall Log Scale: nRF & Accuracy vs SSS)
3. nrf_vs_sss_linear.png (Overall Linear Scale)
"""

import os
import sys
from pathlib import Path
from typing import Dict, Any, Tuple
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from statsmodels.nonparametric.smoothers_lowess import lowess

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
RESULTS_DIR = PROJECT_ROOT / "results"
HIGH_GAP_DIR = RESULTS_DIR / "results_high_gap"

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
    "MSA+NJ": {"color": "#ff7f0e", "label": "MSA+NJ", "linestyle": "-", "linewidth": 2.4, "zorder": 8},
    "PWA+NJ": {"color": "#1f77b4", "label": "PSA+NJ", "linestyle": "-", "linewidth": 2.4, "zorder": 7},
    # FastME (SPR + LG+G)
    "MSA+FastME_LG_G": {"color": "#bcbd22", "label": "MSA+FastME (LG+G)", "linestyle": "-", "linewidth": 2.2, "zorder": 9},
    "PWA+FastME_SPR": {"color": "#e377c2", "label": "PSA+FastME (SPR)", "linestyle": "-", "linewidth": 2.2, "zorder": 9},
    # FastME (Default)
    "MSA+FastME": {"color": "#8c6d31", "label": "MSA+FastME (def)", "linestyle": ":", "linewidth": 1.8, "zorder": 7},
    "PWA+FastME": {"color": "#b5cf6b", "label": "PSA+FastME (def)", "linestyle": ":", "linewidth": 1.8, "zorder": 7},
    # Ground Truth Baselines
    "TRUE_MSA+ML": {"color": "#8c564b", "label": "TRUE_MSA+ML", "linestyle": "--", "linewidth": 1.8, "zorder": 4},
    "TRUE_MSA+NJ": {"color": "#7f7f7f", "label": "TRUE_MSA+NJ", "linestyle": ":", "linewidth": 1.8, "zorder": 3},
    "TRUE_PWA+NJ": {"color": "#17becf", "label": "TRUE_PSA+NJ", "linestyle": "-.", "linewidth": 2.0, "zorder": 5},
}


def compute_lowess_curve(x: np.ndarray, y: np.ndarray, frac: float = 0.3, log_x: bool = True) -> Tuple[np.ndarray, np.ndarray]:
    valid = np.isfinite(x) & np.isfinite(y) & (x > 0 if log_x else True)
    x_v, y_v = x[valid], y[valid]
    if len(x_v) < 5:
        return np.array([]), np.array([])
    sorted_idx = np.argsort(x_v)
    x_s, y_s = x_v[sorted_idx], y_v[sorted_idx]
    if len(x_s) < 10:
        return x_s, y_s

    x_in = np.log10(x_s) if log_x else x_s
    lw = lowess(y_s, x_in, frac=frac, it=3, return_sorted=True)
    if np.allclose(lw[:, 1], y_s) or np.max(np.abs(np.diff(lw[:, 1]))) > 0.08:
        lw = lowess(y_s, x_in, frac=frac, it=0, return_sorted=True)

    out_x = np.power(10.0, lw[:, 0]) if log_x else lw[:, 0]
    return out_x, lw[:, 1]


def plot_condition_stratified(df: pd.DataFrame, out_path: Path, include_default_fastme: bool = True):
    vals = sorted([v for v in df["length"].unique() if v != 100])
    ncols = len(vals)
    fig, axes = plt.subplots(nrows=1, ncols=ncols, figsize=(4.8 * ncols, 4.5), sharey=True)
    axes_flat = axes.flatten() if hasattr(axes, "flatten") else [axes]

    pipeline_order = [
        "MSA+ML", "MSA+NJ", "PWA+NJ",
        "MSA+FastME_LG_G", "PWA+FastME_SPR",
        "MSA+FastME", "PWA+FastME",
        "TRUE_MSA+ML", "TRUE_MSA+NJ", "TRUE_PWA+NJ"
    ]

    for idx, v in enumerate(vals):
        ax = axes_flat[idx]
        sub_v = df[df["length"] == v]

        for pipe in pipeline_order:
            if not include_default_fastme and pipe in ["MSA+FastME", "PWA+FastME"]:
                continue
            sub_p = sub_v[sub_v["pipeline"] == pipe].dropna(subset=["sss_mean", "nrf_distance"])
            if len(sub_p) == 0:
                continue
            cfg = PIPELINE_CONFIGS[pipe]
            x = sub_p["sss_mean"].values
            y = sub_p["nrf_distance"].values
            ax.scatter(x, y, color=cfg["color"], alpha=0.08, s=8, edgecolors="none")

            lx, ly = compute_lowess_curve(x, y, frac=0.35)
            if len(lx) > 0:
                ax.plot(lx, ly, color=cfg["color"], linestyle=cfg["linestyle"], linewidth=cfg["linewidth"],
                        label=cfg["label"], zorder=cfg["zorder"])

        ax.set_title(f"Sequence Length L = {v} aa", fontweight="bold", fontsize=11)
        ax.set_xlabel(r"$SSS$", fontweight="bold")
        ax.set_ylim(-0.02, 0.95)
        ax.set_xscale("log")
        if idx == 0:
            ax.set_ylabel(r"Topological Error: nRF", fontweight="bold")
        if idx == ncols - 1:
            ax.legend(loc="upper right", fontsize=7.5, frameon=True, framealpha=0.9)

    plt.suptitle("High Gap - Stratified by Length - nRF vs SSS (including FastME)", fontweight="bold", fontsize=13, y=1.03)
    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved stratified plot to: {out_path}")


def plot_main_nrf_vs_sss(df: pd.DataFrame, out_path: Path, log_x: bool = True, include_default_fastme: bool = True):
    fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(16.0, 6.2))

    pipeline_order = [
        "MSA+ML", "MSA+NJ", "PWA+NJ",
        "MSA+FastME_LG_G", "PWA+FastME_SPR",
        "MSA+FastME", "PWA+FastME",
        "TRUE_MSA+ML", "TRUE_MSA+NJ", "TRUE_PWA+NJ"
    ]

    for pipe in pipeline_order:
        if not include_default_fastme and pipe in ["MSA+FastME", "PWA+FastME"]:
            continue
        sub = df[df["pipeline"] == pipe].dropna(subset=["sss_mean", "nrf_distance"])
        if len(sub) == 0:
            continue
        cfg = PIPELINE_CONFIGS[pipe]
        x = sub["sss_mean"].values
        y_nrf = sub["nrf_distance"].values
        y_acc = 1.0 - y_nrf

        ax1.scatter(x, y_nrf, color=cfg["color"], alpha=0.06, s=8, edgecolors="none")
        ax2.scatter(x, y_acc, color=cfg["color"], alpha=0.06, s=8, edgecolors="none")

        lx_nrf, ly_nrf = compute_lowess_curve(x, y_nrf, frac=0.3)
        lx_acc, ly_acc = compute_lowess_curve(x, y_acc, frac=0.3)

        if len(lx_nrf) > 0:
            ax1.plot(lx_nrf, ly_nrf, color=cfg["color"], linestyle=cfg["linestyle"], linewidth=cfg["linewidth"],
                     label=cfg["label"], zorder=cfg["zorder"])
        if len(lx_acc) > 0:
            ax2.plot(lx_acc, ly_acc, color=cfg["color"], linestyle=cfg["linestyle"], linewidth=cfg["linewidth"],
                     label=cfg["label"], zorder=cfg["zorder"])

    min_x = max(0.003, df["sss_mean"].min() * 0.8)
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
    ax1.set_title(f"A. Topological Error: nRF vs SSS\nHigh Gap (with FastME) - {scale_str}", fontweight="bold", pad=10)
    ax1.set_ylabel(r"Normalized Robinson-Foulds Distance: nRF", fontweight="bold", fontsize=11)
    ax1.set_ylim(-0.02, 1.0)
    ax1.legend(loc="upper right", frameon=True, framealpha=0.9, fontsize=8.0)

    ax2.set_title(f"B. Topological Accuracy: 1 - nRF vs SSS\nMatsui & Iwasaki (2020) Fig. 3a Replication - {scale_str}", fontweight="bold", pad=10)
    ax2.set_ylabel(r"Topological Accuracy: 1 - nRF", fontweight="bold", fontsize=11)
    ax2.set_ylim(-0.02, 1.02)
    ax2.legend(loc="lower right", frameon=True, framealpha=0.9, fontsize=8.0)

    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    plt.close()
    print(f"Saved main plot to: {out_path}")


def main():
    df_hg = pd.read_csv(HIGH_GAP_DIR / "benchmark_summary.csv")
    df_opt = pd.read_csv(RESULTS_DIR / "results_fastme_options" / "benchmark_fastme_options_summary.csv")
    df_no = pd.read_csv(RESULTS_DIR / "results_fastme_NoOption" / "benchmark_fastme_summary.csv")
    df_no = df_no[df_no["length"] != 100].copy()

    # Combine all
    combined_df = pd.concat([df_hg, df_opt, df_no], ignore_index=True)
    # Filter length 100
    combined_df = combined_df[combined_df["length"] != 100].copy()

    # Save combined summary for future reference
    combined_csv = HIGH_GAP_DIR / "benchmark_summary_with_fastme.csv"
    combined_df.to_csv(combined_csv, index=False)
    print(f"Saved merged dataset: {combined_csv} ({len(combined_df)} rows)")

    # 1. Stratified plot (with all FastME variants) -> replaces nrf_vs_sss_by_condition.png
    plot_condition_stratified(combined_df, HIGH_GAP_DIR / "nrf_vs_sss_by_condition.png", include_default_fastme=True)

    # Also generate variant with only FastME (options) for clean comparison
    plot_condition_stratified(combined_df, HIGH_GAP_DIR / "nrf_vs_sss_by_condition_options_only.png", include_default_fastme=False)

    # 2. Main comparison curves
    plot_main_nrf_vs_sss(combined_df, HIGH_GAP_DIR / "nrf_vs_sss_curves.png", log_x=True, include_default_fastme=True)
    plot_main_nrf_vs_sss(combined_df, HIGH_GAP_DIR / "nrf_vs_sss_linear.png", log_x=False, include_default_fastme=True)


if __name__ == "__main__":
    main()
