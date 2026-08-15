#!/usr/bin/env python3
"""
Sequence Length and Residue Retention Analysis Script
Location: bin/analyze_sequence_lengths.py

Measures the actual remaining amino acid sequence lengths and gap proportions
for simulated sequences (especially L=100 across D=0.1, 0.5, 1.0, 2.0, 3.0).

Supports scanning either results/replications/ directory or Nextflow work/ directory.
Accepts input data path via command line argument (--data / --input_dir / -i).
"""

import os
import re
import sys
import glob
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO

def scan_replications_dir(rep_root):
    """Scans results/replications/D{D}_L{L}_rep{rep}/seqs.fasta"""
    records = []
    condition_dirs = glob.glob(os.path.join(rep_root, "D*_L*_rep*"))
    print(f"Scanning {len(condition_dirs)} replication directories in '{rep_root}'...")

    for cdir in condition_dirs:
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

    return pd.DataFrame(records)

def scan_work_dir(work_root):
    """Scans Nextflow work/ directories for seqs.fasta and parses .command.sh"""
    records = []
    fasta_files = glob.glob(os.path.join(work_root, "**", "seqs.fasta"), recursive=True)
    print(f"Scanning {len(fasta_files)} seqs.fasta files in '{work_root}'...")

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

    return pd.DataFrame(records)

def plot_length_distributions(df, outdir, target_length=100):
    """Plots actual sequence length distribution for target initial length (default: L=100)."""
    sub_df = df[df["initial_length"] == target_length]
    if sub_df.empty:
        print(f"Warning: No data found for initial_length = {target_length}")
        return

    plt.figure(figsize=(10, 6))
    sns.set_theme(style="whitegrid", font_scale=1.1)

    palette = sns.color_palette("viridis", n_colors=len(sub_df["distance"].unique()))

    ax = sns.violinplot(
        data=sub_df,
        x="distance",
        y="actual_length",
        palette=palette,
        inner="quartile",
        cut=0
    )

    # Reference line for initial length
    plt.axhline(target_length, color="red", linestyle="--", linewidth=1.5, label=f"Initial Target Length (L={target_length})")

    # Annotate mean, std, min-max above each violin
    distances = sorted(sub_df["distance"].unique())
    for i, d in enumerate(distances):
        d_vals = sub_df[sub_df["distance"] == d]["actual_length"]
        mean_v = d_vals.mean()
        std_v  = d_vals.std()
        min_v  = d_vals.min()
        max_v  = d_vals.max()
        ax.text(
            i, 
            max_v + 3, 
            f"μ={mean_v:.1f}\nσ={std_v:.1f}\n[{min_v}-{max_v}]", 
            ha="center", 
            va="bottom", 
            fontsize=9.5, 
            fontweight="bold",
            bbox=dict(boxstyle="round,pad=0.2", facecolor="white", alpha=0.8, edgecolor="#ccc")
        )

    plt.title(f"Remaining Amino Acid Sequence Lengths at Initial Target Length L = {target_length} aa\n(Indel Rate = 0.05, 100 Replicates x 16 Taxa)", fontsize=13, fontweight="bold", pad=15)
    plt.xlabel("Evolutionary Distance D (substitutions/site)", fontsize=12, fontweight="bold")
    plt.ylabel("Actual Sequence Length (amino acids)", fontsize=12, fontweight="bold")
    plt.ylim(0, max(sub_df["actual_length"].max() + 25, target_length + 20))
    plt.legend(loc="lower left", frameon=True)
    plt.tight_layout()

    out_png = os.path.join(outdir, f"sequence_lengths_L{target_length}.png")
    out_pdf = os.path.join(outdir, f"sequence_lengths_L{target_length}.pdf")
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.savefig(out_pdf, bbox_inches="tight")
    plt.close()
    print(f"Saved length distribution plot: {out_png}")
    print(f"Saved length distribution plot (PDF): {out_pdf}")

def plot_all_lengths_summary(df, outdir):
    """Plots faceted sequence length summary for all lengths (100, 300, 500, 1000, 1500)."""
    lengths = sorted(df["initial_length"].unique())
    fig, axes = plt.subplots(1, len(lengths), figsize=(22, 5.5), sharey=False)

    for i, l in enumerate(lengths):
        ax = axes[i]
        sub_df = df[df["initial_length"] == l]
        sns.boxplot(
            data=sub_df,
            x="distance",
            y="actual_length",
            palette="Blues",
            ax=ax
        )
        ax.axhline(l, color="red", linestyle="--", linewidth=1.2, label=f"Target L={l}")
        ax.set_title(f"Initial L = {l} aa", fontsize=12, fontweight="bold")
        ax.set_xlabel("Distance D", fontsize=11, fontweight="bold")
        if i == 0:
            ax.set_ylabel("Actual Sequence Length (aa)", fontsize=11, fontweight="bold")
        else:
            ax.set_ylabel("")
        ax.grid(True, linestyle="--", alpha=0.6)

    fig.suptitle("Actual Sequence Length vs Evolutionary Distance (D) across All Target Lengths (L)", fontsize=14, fontweight="bold", y=1.02)
    plt.tight_layout()

    out_png = os.path.join(outdir, "sequence_lengths_all_conditions.png")
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved all-lengths summary plot: {out_png}")

def generate_summary_table(df, outdir):
    """Generates summary statistics table of sequence lengths."""
    grouped = df.groupby(["initial_length", "distance"])["actual_length"]
    summary_df = grouped.agg(
        taxa_count="count",
        mean_length="mean",
        std_length="std",
        min_length="min",
        median_length="median",
        max_length="max",
        q25=lambda x: np.percentile(x, 25),
        q75=lambda x: np.percentile(x, 75)
    ).reset_index()

    summary_df["retention_pct"] = (summary_df["mean_length"] / summary_df["initial_length"]) * 100.0
    summary_df["mean_length"] = summary_df["mean_length"].round(2)
    summary_df["std_length"]  = summary_df["std_length"].round(2)
    summary_df["retention_pct"] = summary_df["retention_pct"].round(1)

    csv_path = os.path.join(outdir, "sequence_length_statistics.csv")
    summary_df.to_csv(csv_path, index=False)
    print(f"Saved sequence length statistics table: {csv_path}")

    # Print L=100 highlight table
    l100_table = summary_df[summary_df["initial_length"] == 100]
    print("\n" + "=" * 65)
    print("=== Sequence Length Summary for Initial L = 100 aa ===")
    print(l100_table[["distance", "mean_length", "std_length", "min_length", "median_length", "max_length", "retention_pct"]].to_string(index=False))
    print("=" * 65)

def main():
    parser = argparse.ArgumentParser(
        description="Analyze remaining amino acid sequence lengths across evolutionary distances (D) and target lengths (L)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Scan organized replications directory:
  python bin/analyze_sequence_lengths.py --data results/replications --outdir results/length_analysis

  # Scan raw Nextflow work directory directly:
  python bin/analyze_sequence_lengths.py --data work --outdir results/length_analysis
        """
    )
    parser.add_argument("--data", "-i", "--input_dir", dest="input_dir", required=True, help="Path to data directory (e.g. results/replications or work)")
    parser.add_argument("--outdir", "-o", default="results/length_analysis", help="Output directory for plots and tables (default: results/length_analysis)")
    args = parser.parse_args()

    input_dir = os.path.abspath(args.input_dir)
    outdir = os.path.abspath(args.outdir)

    print(f"=== Sequence Length Analysis ===")
    print(f"Input data directory : {input_dir}")
    print(f"Output directory     : {outdir}")
    print("=" * 38)

    if not os.path.exists(input_dir):
        print(f"Error: Specified data directory does not exist: {input_dir}")
        sys.exit(1)

    os.makedirs(outdir, exist_ok=True)

    # Check if input is a replications directory or a work directory
    sample_subdirs = glob.glob(os.path.join(input_dir, "D*_L*_rep*"))
    if sample_subdirs:
        print("Detected format: results/replications/ structure")
        df = scan_replications_dir(input_dir)
    else:
        print("Detected format: Nextflow work/ directory structure")
        df = scan_work_dir(input_dir)

    if df.empty:
        print(f"Error: No valid sequence files (seqs.fasta) found under '{input_dir}'.")
        sys.exit(1)

    total_taxa = len(df)
    unique_reps = len(df["replicate"].unique())
    print(f"Successfully loaded {total_taxa} sequence length records across {unique_reps} replicates.")

    # 1. Plot L=100 detailed distribution (Violin plot with μ, σ, min, max)
    plot_length_distributions(df, outdir, target_length=100)

    # 2. Plot all lengths (100, 300, 500, 1000, 1500)
    plot_all_lengths_summary(df, outdir)

    # 3. Generate summary table
    generate_summary_table(df, outdir)

    print("\n" + "=" * 38)
    print(f"=== Analysis Completed! ===")
    print(f"All figures and summary tables saved to: {outdir}")
    print("=" * 38)

if __name__ == "__main__":
    main()
