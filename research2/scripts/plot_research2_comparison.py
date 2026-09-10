#!/usr/bin/env python3
"""
plot_research2_comparison.py

Generates publication-quality comparative plots and summary tables
evaluating AliSim vs INDELible across 4 conditions:
  - INDELible (rate=0.10/0.10)
  - INDELible (rate=0.05/0.05)
  - AliSim (rate=0.10/0.10)
  - AliSim (rate=0.05/0.05)
across sequence lengths [500, 1000, 1500] and distances [0.1, 0.2, 0.5, 1.0, 1.5, 2.0].
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# Style configuration
plt.rcParams.update({
    "font.size": 11,
    "font.family": "sans-serif",
    "axes.labelsize": 12,
    "axes.titlesize": 13,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 10,
    "figure.titlesize": 15,
    "lines.linewidth": 2.0,
    "lines.markersize": 6,
    "grid.alpha": 0.4,
    "grid.linestyle": "--"
})

# Curated high-contrast palette
CONDITION_STYLES = {
    "indelible_rate0.10": {
        "label": "INDELible (0.10)",
        "color": "#D32F2F",   # Strong Red
        "marker": "o",
        "linestyle": "-"
    },
    "indelible_rate0.05": {
        "label": "INDELible (0.05)",
        "color": "#F57C00",   # Orange
        "marker": "s",
        "linestyle": "--"
    },
    "alisim_rate0.10": {
        "label": "AliSim (0.10)",
        "color": "#1976D2",   # Strong Blue
        "marker": "^",
        "linestyle": "-"
    },
    "alisim_rate0.05": {
        "label": "AliSim (0.05)",
        "color": "#388E3C",   # Green
        "marker": "d",
        "linestyle": "--"
    }
}


def plot_fig1_geometric_properties(df, out_dir):
    """
    Figure 1: True MSA structural properties across evolutionary distance.
    4 panels: Alignment length, Gap ratio %, Mean sequence length, Indel block count.
    """
    lengths = sorted(df["root_length"].unique())
    fig, axes = plt.subplots(2, 2, figsize=(14, 11))
    
    # We plot with root_length = 1000 as representative, with annotations/variants
    target_len = 1000 if 1000 in lengths else lengths[0]
    df_sub = df[df["root_length"] == target_len]

    conds = list(CONDITION_STYLES.keys())
    dists = sorted(df_sub["distance"].unique())

    # (a) True MSA Length (Columns)
    ax = axes[0, 0]
    for cond in conds:
        d_c = df_sub[df_sub["cond_name"] == cond]
        if d_c.empty: continue
        means = d_c.groupby("distance")["true_msa_length"].mean()
        stds = d_c.groupby("distance")["true_msa_length"].std()
        cfg = CONDITION_STYLES[cond]
        ax.errorbar(means.index, means.values, yerr=stds.values, label=cfg["label"],
                     color=cfg["color"], marker=cfg["marker"], linestyle=cfg["linestyle"], capsize=3)
    ax.set_title(f"(a) True MSA Column Count (L={target_len})")
    ax.set_xlabel("Evolutionary Distance (D)")
    ax.set_ylabel("True MSA Length (columns)")
    ax.grid(True)
    ax.legend()

    # (b) Gap Percentage
    ax = axes[0, 1]
    for cond in conds:
        d_c = df_sub[df_sub["cond_name"] == cond]
        if d_c.empty: continue
        means = d_c.groupby("distance")["true_gap_pct"].mean()
        stds = d_c.groupby("distance")["true_gap_pct"].std()
        cfg = CONDITION_STYLES[cond]
        ax.errorbar(means.index, means.values, yerr=stds.values, label=cfg["label"],
                     color=cfg["color"], marker=cfg["marker"], linestyle=cfg["linestyle"], capsize=3)
    ax.set_title(f"(b) True MSA Gap Percentage (L={target_len})")
    ax.set_xlabel("Evolutionary Distance (D)")
    ax.set_ylabel("Gap Percentage (%)")
    ax.grid(True)
    ax.legend()

    # (c) Mean Sequence Length vs Initial Length
    ax = axes[1, 0]
    for cond in conds:
        d_c = df_sub[df_sub["cond_name"] == cond]
        if d_c.empty: continue
        means = d_c.groupby("distance")["mean_seq_len"].mean()
        stds = d_c.groupby("distance")["mean_seq_len"].std()
        cfg = CONDITION_STYLES[cond]
        ax.errorbar(means.index, means.values, yerr=stds.values, label=cfg["label"],
                     color=cfg["color"], marker=cfg["marker"], linestyle=cfg["linestyle"], capsize=3)
    ax.axhline(target_len, color="gray", linestyle=":", label="Initial Root Length")
    ax.set_title(f"(c) Mean Sequence Length (Residues, L={target_len})")
    ax.set_xlabel("Evolutionary Distance (D)")
    ax.set_ylabel("Mean Sequence Length (aa)")
    ax.grid(True)
    ax.legend()

    # (d) Number of Indel Blocks
    ax = axes[1, 1]
    for cond in conds:
        d_c = df_sub[df_sub["cond_name"] == cond]
        if d_c.empty: continue
        means = d_c.groupby("distance")["true_num_indels"].mean()
        stds = d_c.groupby("distance")["true_num_indels"].std()
        cfg = CONDITION_STYLES[cond]
        ax.errorbar(means.index, means.values, yerr=stds.values, label=cfg["label"],
                     color=cfg["color"], marker=cfg["marker"], linestyle=cfg["linestyle"], capsize=3)
    ax.set_title(f"(d) Indel Block Count (L={target_len})")
    ax.set_xlabel("Evolutionary Distance (D)")
    ax.set_ylabel("Number of Indel Blocks")
    ax.grid(True)
    ax.legend()

    fig.suptitle("Figure 1: Geometric Properties of True MSA (AliSim vs INDELible)", y=0.99)
    plt.tight_layout()
    out_path = os.path.join(out_dir, "fig1_msa_geometric_properties.png")
    plt.savefig(out_path, dpi=300)
    plt.close()
    print(f"Saved: {out_path}")


def plot_fig2_alignment_accuracy(df, out_dir):
    """
    Figure 2: Alignment Accuracy (Developer's score f_D) and SSS.
    (a) Distance vs f_D
    (b) SSS (log scale) vs f_D
    """
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    lengths = sorted(df["root_length"].unique())
    target_len = 1000 if 1000 in lengths else lengths[0]
    df_sub = df[df["root_length"] == target_len].dropna(subset=["mafft_fd_score"])

    conds = list(CONDITION_STYLES.keys())

    # (a) Distance vs f_D
    ax = axes[0]
    for cond in conds:
        d_c = df_sub[df_sub["cond_name"] == cond]
        if d_c.empty: continue
        means = d_c.groupby("distance")["mafft_fd_score"].mean()
        stds = d_c.groupby("distance")["mafft_fd_score"].std()
        cfg = CONDITION_STYLES[cond]
        ax.errorbar(means.index, means.values, yerr=stds.values, label=cfg["label"],
                     color=cfg["color"], marker=cfg["marker"], linestyle=cfg["linestyle"], capsize=3)
    ax.set_title("(a) Distance vs MAFFT Alignment Accuracy ($f_D$)")
    ax.set_xlabel("Evolutionary Distance (D)")
    ax.set_ylabel("Developer's Score ($f_D$)")
    ax.set_ylim(-0.02, 1.02)
    ax.grid(True)
    ax.legend()

    # (b) SSS (log scale) vs f_D
    ax = axes[1]
    has_sss = "sss_mean" in df_sub.columns and df_sub["sss_mean"].notna().any()
    if has_sss:
        for cond in conds:
            d_c = df_sub[df_sub["cond_name"] == cond].dropna(subset=["sss_mean"])
            if d_c.empty: continue
            cfg = CONDITION_STYLES[cond]
            # Scatter with small jitter or group mean
            means = d_c.groupby("distance")[["sss_mean", "mafft_fd_score"]].mean()
            ax.plot(means["sss_mean"], means["mafft_fd_score"], label=cfg["label"],
                    color=cfg["color"], marker=cfg["marker"], linestyle=cfg["linestyle"])
            ax.scatter(d_c["sss_mean"], d_c["mafft_fd_score"], color=cfg["color"], alpha=0.3, s=20)
        ax.set_xscale("log")
        ax.set_title("(b) Sequence Similarity (SSS, log) vs MAFFT $f_D$")
        ax.set_xlabel("Sequence Similarity Score (SSS, log scale)")
        ax.set_ylabel("Developer's Score ($f_D$)")
        ax.set_ylim(-0.02, 1.02)
        ax.grid(True, which="both")
        ax.legend()
    else:
        ax.text(0.5, 0.5, "SSS data pending", ha="center", va="center")

    fig.suptitle(f"Figure 2: Alignment Accuracy Evaluation (L={target_len})", y=1.02)
    plt.tight_layout()
    out_path = os.path.join(out_dir, "fig2_alignment_accuracy_vs_similarity.png")
    plt.savefig(out_path, dpi=300)
    plt.close()
    print(f"Saved: {out_path}")


def plot_fig3_tree_accuracy(df, out_dir):
    """
    Figure 3: Phylogenetic Tree Reconstruction Accuracy (RCRB = 1 - nRF).
    Panels by root sequence length (500, 1000, 1500).
    Compares MSA+ML, TRUE_MSA+ML, and PWA+NJ across the 4 conditions.
    """
    lengths = sorted(df["root_length"].unique())
    n_len = len(lengths)
    fig, axes = plt.subplots(1, n_len, figsize=(6 * n_len, 6), sharey=True)
    if n_len == 1:
        axes = [axes]

    conds = list(CONDITION_STYLES.keys())

    for idx, l in enumerate(lengths):
        ax = axes[idx]
        df_l = df[df["root_length"] == l]

        for cond in conds:
            d_c = df_l[df_l["cond_name"] == cond]
            if d_c.empty: continue
            cfg = CONDITION_STYLES[cond]
            
            # MSA+ML (Solid line)
            if "rcrb_msa_ml" in d_c.columns and d_c["rcrb_msa_ml"].notna().any():
                means = d_c.groupby("distance")["rcrb_msa_ml"].mean()
                ax.plot(means.index, means.values, label=f"{cfg['label']} MSA+ML",
                        color=cfg["color"], marker=cfg["marker"], linestyle="-")
                
            # TRUE_MSA+ML (Dotted line with same color)
            if "rcrb_true_msa_ml" in d_c.columns and d_c["rcrb_true_msa_ml"].notna().any():
                means_true = d_c.groupby("distance")["rcrb_true_msa_ml"].mean()
                ax.plot(means_true.index, means_true.values, label=f"{cfg['label']} TRUE+ML",
                        color=cfg["color"], marker=cfg["marker"], linestyle=":", alpha=0.6)

        ax.set_title(f"Sequence Length L = {l}")
        ax.set_xlabel("Evolutionary Distance (D)")
        if idx == 0:
            ax.set_ylabel("Tree Accuracy (RCRB = 1 - nRF)")
        ax.set_ylim(-0.02, 1.02)
        ax.grid(True)
        if idx == n_len - 1:
            ax.legend(bbox_to_anchor=(1.05, 1), loc="upper left", fontsize=9)

    fig.suptitle("Figure 3: Phylogenetic Tree Accuracy (MSA+ML vs True MSA+ML)", y=1.02)
    plt.tight_layout()
    out_path = os.path.join(out_dir, "fig3_tree_reconstruction_accuracy.png")
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved: {out_path}")


def plot_fig4_alignment_loss(df, out_dir):
    """
    Figure 4: Alignment Loss Delta = RCRB(TRUE_MSA+ML) - RCRB(MSA+ML).
    Direct quantification of error caused by alignment failure.
    """
    lengths = sorted(df["root_length"].unique())
    target_len = 1000 if 1000 in lengths else lengths[0]
    df_sub = df[df["root_length"] == target_len].dropna(subset=["alignment_loss_ml"])

    fig, ax = plt.subplots(figsize=(8, 6))
    conds = list(CONDITION_STYLES.keys())

    for cond in conds:
        d_c = df_sub[df_sub["cond_name"] == cond]
        if d_c.empty: continue
        means = d_c.groupby("distance")["alignment_loss_ml"].mean()
        stds = d_c.groupby("distance")["alignment_loss_ml"].std()
        cfg = CONDITION_STYLES[cond]
        ax.errorbar(means.index, means.values, yerr=stds.values, label=cfg["label"],
                     color=cfg["color"], marker=cfg["marker"], linestyle=cfg["linestyle"], capsize=3)

    ax.set_title(f"Figure 4: Alignment Loss Δ [True ML - MSA ML] (L={target_len})")
    ax.set_xlabel("Evolutionary Distance (D)")
    ax.set_ylabel("Alignment Loss (Accuracy Drop due to Alignment Error)")
    ax.grid(True)
    ax.legend()

    plt.tight_layout()
    out_path = os.path.join(out_dir, "fig4_alignment_loss_analysis.png")
    plt.savefig(out_path, dpi=300)
    plt.close()
    print(f"Saved: {out_path}")


def generate_summary_tables(df, out_dir):
    """Generates detailed summary table (CSV and Markdown)."""
    metrics = [
        "true_msa_length", "true_gap_pct", "mean_seq_len", "true_num_indels",
        "mafft_fd_score", "sss_mean", "rcrb_msa_ml", "rcrb_true_msa_ml", "alignment_loss_ml"
    ]
    avail_metrics = [m for m in metrics if m in df.columns]

    agg_dict = {m: ["mean", "std"] for m in avail_metrics}
    summary = df.groupby(["cond_name", "root_length", "distance"]).agg(agg_dict).round(4)
    
    csv_path = os.path.join(out_dir, "summary_table_research2.csv")
    summary.to_csv(csv_path)
    print(f"Summary table written to: {csv_path}")

    # Generate Markdown Report
    md_path = os.path.join(out_dir, "summary_report.md")
    with open(md_path, "w") as f:
        f.write("# AliSim vs INDELible Benchmark Summary Report (research2)\n\n")
        f.write(f"Total simulated runs analyzed: **{len(df)}**\n\n")
        f.write("### Key Highlights at D = 2.0 (Extreme Distance, L = 1000)\n\n")
        
        d2 = df[(df["distance"] == 2.0) & (df["root_length"] == 1000)]
        if not d2.empty:
            d2_summary = d2.groupby("cond_name").agg({
                "true_msa_length": "mean",
                "true_gap_pct": "mean",
                "mafft_fd_score": "mean",
                "rcrb_msa_ml": "mean",
                "alignment_loss_ml": "mean"
            }).round(3)
            f.write(d2_summary.to_markdown())
            f.write("\n\n")
        
        f.write("### Complete Aggregate Metrics\n\n")
        f.write("Please refer to `summary_table_research2.csv` for the full breakdown.\n")
    print(f"Summary report written to: {md_path}")


def main():
    parser = argparse.ArgumentParser(description="Plot comparative results for research2.")
    parser.add_argument("--csv", required=True, help="Input concatenated CSV file or pattern")
    parser.add_argument("--outdir", required=True, help="Output directory for plots")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    df = pd.read_csv(args.csv)
    print(f"Loaded {len(df)} rows from {args.csv}")

    plot_fig1_geometric_properties(df, args.outdir)
    plot_fig2_alignment_accuracy(df, args.outdir)
    plot_fig3_tree_accuracy(df, args.outdir)
    plot_fig4_alignment_loss(df, args.outdir)
    generate_summary_tables(df, args.outdir)


if __name__ == "__main__":
    main()
