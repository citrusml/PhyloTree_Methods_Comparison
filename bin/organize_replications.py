#!/usr/bin/env python3
"""
Organize Replications Script
Scans the Nextflow work/ directory, identifies all intermediate outputs
(true trees, simulated fasta, PWA trees/matrices, MSA trees/alignments/matrices, ML trees/metadata),
and organizes them cleanly into results/replications/D{dist}_L{len}_rep{rep}/.

Supports multi-threaded copying, hard-linking (instant & zero disk overhead), or symlinking.
"""

import os
import re
import sys
import glob
import shutil
import argparse
from concurrent.futures import ThreadPoolExecutor
from collections import defaultdict

def parse_rep_key(csv_filename):
    """
    Extracts (dist, len, rep, pipeline) from CSV filenames like:
    pwa_nj_D0.1_L100_rep1.csv, msa_nj_D1.0_L500_rep12.csv, msa_ml_D3.0_L1500_rep100.csv
    """
    pattern = r"^(pwa_nj|msa_nj|msa_ml)_D([\d\.]+)_L(\d+)_rep(\d+)\.csv$"
    match = re.match(pattern, csv_filename)
    if match:
        pipeline, dist_str, len_str, rep_str = match.groups()
        return pipeline, dist_str, len_str, rep_str
    return None

def parse_simulate_data_dir(dirpath):
    """
    Parses .command.sh in SIMULATE_DATA work directories to find (dist, len, rep).
    """
    cmd_sh = os.path.join(dirpath, ".command.sh")
    if not os.path.exists(cmd_sh):
        return None

    try:
        with open(cmd_sh, "r") as f:
            content = f.read()

        dist_m = re.search(r"--distance\s+([\d\.]+)", content)
        len_m  = re.search(r"--length\s+(\d+)", content)
        rep_m  = re.search(r"--seed\s+(\d+)", content)

        if dist_m and len_m and rep_m:
            return dist_m.group(1), len_m.group(1), rep_m.group(1)
    except Exception:
        pass
    return None

def transfer_file(src, dst, mode="copy"):
    """Transfers file via copy, hardlink, or symlink."""
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    if os.path.exists(dst):
        try:
            os.remove(dst)
        except OSError:
            pass

    if mode == "hardlink":
        try:
            os.link(src, dst)
            return True
        except OSError:
            mode = "copy"  # Fallback to copy if cross-device

    if mode == "symlink":
        try:
            os.symlink(os.path.abspath(src), dst)
            return True
        except OSError:
            mode = "copy"

    if mode == "copy":
        shutil.copy2(src, dst)
        return True

    return False

def main():
    parser = argparse.ArgumentParser(description="Harvest and organize all trees and MSAs from work/ to results/replications/")
    parser.add_argument("--workdir", default="work", help="Path to Nextflow work directory (default: work)")
    parser.add_argument("--outdir", default="results/replications", help="Output directory for organized replications (default: results/replications)")
    parser.add_argument("--mode", choices=["copy", "hardlink", "symlink"], default="copy", help="File transfer mode: copy, hardlink (fastest & zero disk usage), or symlink (default: copy)")
    parser.add_argument("--threads", type=int, default=16, help="Number of worker threads (default: 16)")
    args = parser.parse_args()

    print(f"=== PhyloMethod Results Organizer ===")
    print(f"Scanning work directory : {args.workdir}")
    print(f"Destination directory   : {args.outdir}")
    print(f"Transfer mode           : {args.mode}")
    print(f"Worker threads          : {args.threads}")
    print("=" * 38)

    if not os.path.exists(args.workdir):
        print(f"Error: Work directory '{args.workdir}' does not exist.")
        sys.exit(1)

    os.makedirs(args.outdir, exist_ok=True)

    # 1. Discover all CSV files to map work directories
    print("Scanning work directories for task outputs...")
    csv_files = glob.glob(os.path.join(args.workdir, "**", "*.csv"), recursive=True)
    print(f"Found {len(csv_files)} task CSV files.")

    # Storage map: rep_key -> dict of (target_filename -> src_path)
    replication_files = defaultdict(dict)

    for csv_path in csv_files:
        fname = os.path.basename(csv_path)
        parsed = parse_rep_key(fname)
        if not parsed:
            continue

        pipeline, dist_str, len_str, rep_str = parsed
        rep_key = f"D{dist_str}_L{len_str}_rep{rep_str}"
        task_dir = os.path.dirname(csv_path)

        if pipeline == "pwa_nj":
            replication_files[rep_key]["pwa_nj.csv"] = csv_path
            pwa_nwk = os.path.join(task_dir, "pwa_nj.nwk")
            if os.path.exists(pwa_nwk):
                replication_files[rep_key]["pwa_nj.nwk"] = pwa_nwk
            pwa_mat = os.path.join(task_dir, "pwa_matrix.phylip")
            if os.path.exists(pwa_mat):
                replication_files[rep_key]["pwa_matrix.phylip"] = pwa_mat

        elif pipeline == "msa_nj":
            replication_files[rep_key]["msa_nj.csv"] = csv_path
            msa_nwk = os.path.join(task_dir, "msa_nj.nwk")
            if os.path.exists(msa_nwk):
                replication_files[rep_key]["msa_nj.nwk"] = msa_nwk
            msa_fa = os.path.join(task_dir, "msa.fasta")
            if os.path.exists(msa_fa):
                replication_files[rep_key]["msa.fasta"] = msa_fa
            msa_mat = os.path.join(task_dir, "msa_matrix.phylip")
            if os.path.exists(msa_mat):
                replication_files[rep_key]["msa_matrix.phylip"] = msa_mat

        elif pipeline == "msa_ml":
            replication_files[rep_key]["msa_ml.csv"] = csv_path
            ml_nwk = os.path.join(task_dir, "msa_ml.nwk")
            if os.path.exists(ml_nwk):
                replication_files[rep_key]["msa_ml.nwk"] = ml_nwk
            ml_json = os.path.join(task_dir, "msa_ml_meta.json")
            if os.path.exists(ml_json):
                replication_files[rep_key]["msa_ml_meta.json"] = ml_json

    # 2. Discover SIMULATE_DATA directories (true_tree.nwk, seqs.fasta)
    print("Locating simulated true trees and sequences...")
    true_tree_files = glob.glob(os.path.join(args.workdir, "**", "true_tree.nwk"), recursive=True)
    print(f"Found {len(true_tree_files)} simulation directories.")

    for tt_path in true_tree_files:
        task_dir = os.path.dirname(tt_path)
        parsed = parse_simulate_data_dir(task_dir)
        if parsed:
            dist_str, len_str, rep_str = parsed
            rep_key = f"D{dist_str}_L{len_str}_rep{rep_str}"
            replication_files[rep_key]["true_tree.nwk"] = tt_path
            seq_fa = os.path.join(task_dir, "seqs.fasta")
            if os.path.exists(seq_fa):
                replication_files[rep_key]["seqs.fasta"] = seq_fa

    total_conditions = len(replication_files)
    print(f"\nTotal unique conditions identified: {total_conditions}")

    # Build flat task list for parallel transfer
    tasks = []
    for rep_key, file_dict in replication_files.items():
        dest_dir = os.path.join(args.outdir, rep_key)
        for target_name, src_file in file_dict.items():
            dst_file = os.path.join(dest_dir, target_name)
            tasks.append((src_file, dst_file, args.mode))

    print(f"Total files to transfer: {len(tasks)}")
    print(f"Starting parallel transfer ({args.mode})...")

    completed = 0
    with ThreadPoolExecutor(max_workers=args.threads) as executor:
        futures = [executor.submit(transfer_file, src, dst, mode) for src, dst, mode in tasks]
        for f in futures:
            if f.result():
                completed += 1
            if completed % 2000 == 0 or completed == len(tasks):
                print(f"  Progress: {completed} / {len(tasks)} files transferred ({completed / len(tasks) * 100:.1f}%)")

    # 3. Verification and summary
    print("\n" + "=" * 38)
    print("=== Replications Organization Summary ===")
    print(f"Total replications processed : {total_conditions}")
    print(f"Total files organized        : {completed} / {len(tasks)}")
    print(f"Output directory             : {os.path.abspath(args.outdir)}")
    print("=" * 38)
    print("\nSample organized replication content:")
    sample_key = next(iter(replication_files.keys())) if replication_files else None
    if sample_key:
        sample_path = os.path.join(args.outdir, sample_key)
        print(f"  {sample_path}/")
        for fn in sorted(os.listdir(sample_path)):
            print(f"    ├── {fn}")
    print("\nOrganization complete! You can now safely archive or inspect results.")

if __name__ == "__main__":
    main()
