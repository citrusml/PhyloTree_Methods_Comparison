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
import shutil
import argparse
from concurrent.futures import ThreadPoolExecutor
from collections import defaultdict
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
    "MSA+ML": "#2ca02c",       # Green
    "TRUE_MSA+NJ": "#9467bd",  # Purple
    "TRUE_MSA+ML": "#8c564b",  # Brown
    "TRUE_DIST+NJ": "#d62728"  # Red
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
    """Generates 5x5 grid of nRF distribution boxplots with mean annotations."""
    distances = sorted(df["distance"].unique())
    lengths = sorted(df["length"].unique())
    pipes = [p for p in ["PWA+NJ", "MSA+NJ", "MSA+ML", "TRUE_MSA+NJ", "TRUE_MSA+ML", "TRUE_DIST+NJ"] if p in df["pipeline"].unique()]

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

    palette = {p: PIPELINE_COLORS.get(p, "#333333") for p in pipes}

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
# 3. Replication Artifact Harvesting
# ==========================================

def organize_replications_from_work(work_dir, rep_outdir, mode="copy", threads=16):
    """
    Harvests all tree files, FASTA, MSAs, and matrices from work/ into results/replications/
    """
    try:
        from bin.organize_replications import organize_replications
        organize_replications(work_dir, rep_outdir, mode=mode, threads=threads)
    except ImportError:
        # Fallback if imported from different path
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        from organize_replications import organize_replications
        organize_replications(work_dir, rep_outdir, mode=mode, threads=threads)

# ==========================================
# Main CLI Entry Point
# ==========================================

def main():
    parser = argparse.ArgumentParser(description="Integrated Benchmark Reporting, Sequence Length Analysis, and Replication Organizer")
    parser.add_argument("--csv", help="Path to aggregated benchmark_summary.csv")
    parser.add_argument("--workdir", help="Path to Nextflow work/ directory (optional: enables auto-aggregation & organization)")
    parser.add_argument("--repdir", help="Path to results/replications/ directory (optional: for sequence length analysis)")
    parser.add_argument("--outdir", default="results", help="Output directory for reports and figures (default: results)")
    parser.add_argument("--organize_replications", action="store_true", help="Harvest and organize all replication files from work/ into results/replications/")
    parser.add_argument("--mode", choices=["copy", "hardlink"], default="copy", help="Transfer mode for replication files (default: copy)")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # 1. Phylogenetic Benchmark Analysis
    csv_file = args.csv
    if not csv_file and args.workdir:
        # Auto-discover CSV
        cand = os.path.join(args.outdir, "benchmark_summary.csv")
        if os.path.exists(cand):
            csv_file = cand

    if csv_file and os.path.exists(csv_file) and os.path.getsize(csv_file) > 0:
        df = pd.read_csv(csv_file)
        if not df.empty and "nrf_distance" in df.columns:
            print(f"Generating phylogenetic benchmark reports from '{csv_file}' ({len(df)} records)...")
            generate_regime_map(df, args.outdir)
            generate_method_comparisons(df, args.outdir)
            generate_boxplots(df, args.outdir)
            generate_summary_table(df, args.outdir)

    # 2. Sequence Length Analysis
    scan_source = args.repdir or (os.path.join(args.outdir, "replications") if os.path.exists(os.path.join(args.outdir, "replications")) else args.workdir)
    if scan_source and os.path.exists(scan_source):
        analyze_sequence_lengths(scan_source, args.outdir)

    # 3. Replication Artifact Harvesting
    if args.organize_replications and args.workdir and os.path.exists(args.workdir):
        rep_dst = args.repdir or os.path.join(args.outdir, "replications")
        organize_replications_from_work(args.workdir, rep_dst, mode=args.mode)

    print(f"\nAll post-processing reports and analyses completed in '{args.outdir}'!")

if __name__ == "__main__":
    main()
