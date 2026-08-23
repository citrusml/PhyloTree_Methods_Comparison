#!/usr/bin/env python3
"""
Replication Artifact Harvesting & Organization Script
Location: bin/organize_replications.py

Harvests all tree files, unaligned FASTA, MSAs, distance matrices, metadata, and evaluation CSVs
from Nextflow work/ directory into structured replication directories:
  results/replications/D{dist}_L{len}_rep{rep}/  (for Experiment 1)
  or
  results/results_taxon/replications/N{taxa}_D{dist}_L{len}_rep{rep}/  (for Experiment 2)
"""

import os
import re
import sys
import glob
import shutil
import argparse
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor

TARGET_FILES = [
    # Ground Truth Data
    "true_tree.nwk",
    "seqs.fasta",
    "true_msa.fasta",
    "true_matrix.phylip",

    # PWA Pipeline
    "pwa_nj.nwk",
    "pwa_matrix.phylip",

    # Shared MAFFT & MSA+NJ / MSA+ML Pipelines
    "msa.fasta",
    "msa_nj.nwk",
    "msa_matrix.phylip",
    "msa_ml.nwk",
    "msa_ml_meta.json",

    # True MSA Pipelines (Experiment 1 optional)
    "true_msa_nj.nwk",
    "true_msa_matrix.phylip",
    "true_msa_ml.nwk",
    "true_msa_ml_meta.json",

    # True Patristic Distance Pipeline (Experiment 1 optional)
    "true_dist_nj.nwk"
]

def extract_rep_key_from_task_dir(task_dir):
    """
    Extracts replication key (e.g., 'D0.1_L100_rep1' or 'N8_D0.1_L100_rep1')
    from task directory metadata (.command.sh, .command.run, or evaluation CSVs).
    """
    # 1. Try extracting from evaluation CSV filenames in task directory
    csv_files = glob.glob(os.path.join(task_dir, "*.csv"))
    for cf in csv_files:
        bname = os.path.basename(cf)
        # Taxon experiment format: ..._N8_D0.1_L100_rep1.csv
        m_taxa = re.search(r"_N(\d+)_D([\d\.]+)_L(\d+)_rep(\d+)\.csv$", bname)
        if m_taxa:
            return f"N{m_taxa.group(1)}_D{m_taxa.group(2)}_L{m_taxa.group(3)}_rep{m_taxa.group(4)}"
        # Standard experiment format: ..._D0.1_L100_rep1.csv
        m_std = re.search(r"_D([\d\.]+)_L(\d+)_rep(\d+)\.csv$", bname)
        if m_std:
            return f"D{m_std.group(1)}_L{m_std.group(2)}_rep{m_std.group(3)}"

    # 2. Try extracting from .command.run (contains Nextflow task tag)
    cmd_run = os.path.join(task_dir, ".command.run")
    if os.path.exists(cmd_run):
        try:
            with open(cmd_run, "r", errors="ignore") as rf:
                content = rf.read()
            # Tag match: D=0.1_L=100_rep=1 or N=8_D=0.1_L=100_rep=1
            m_tag = re.search(r"(?:NXF_TASK_TAG|# NEXTFLOW TASK:[^(\n]*\()\s*['\"]?(?:N=)?(\d+)?_?D=([\d\.]+)_L=(\d+)_rep=(\d+)['\"]?", content)
            if m_tag:
                taxa, dist, length, rep = m_tag.groups()
                if taxa:
                    return f"N{taxa}_D{dist}_L{length}_rep{rep}"
                return f"D{dist}_L{length}_rep{rep}"
        except Exception:
            pass

    # 3. Try extracting from .command.sh (contains CLI arguments)
    cmd_sh = os.path.join(task_dir, ".command.sh")
    if os.path.exists(cmd_sh):
        try:
            with open(cmd_sh, "r", errors="ignore") as cf:
                txt = cf.read()
            taxa_m = re.search(r"--num_taxa\s+(\d+)", txt)
            dist_m = re.search(r"--distance\s+([\d\.]+)", txt)
            len_m  = re.search(r"--length\s+(\d+)", txt)
            rep_m  = re.search(r"--(?:seed|replicate)\s+(\d+)", txt)
            if dist_m and len_m and rep_m:
                dist = dist_m.group(1)
                length = len_m.group(1)
                rep = rep_m.group(1)
                if taxa_m:
                    return f"N{taxa_m.group(1)}_D{dist}_L{length}_rep{rep}"
                return f"D{dist}_L{length}_rep{rep}"
        except Exception:
            pass

    return None

def harvest_work_artifacts(work_dir):
    """
    Scans Nextflow work/ directory and indexes all artifact files by replication key.
    """
    print(f"Scanning Nextflow work directory: '{work_dir}'...")
    tasks_by_rep = defaultdict(dict)

    # Search for all target file names in work/
    for target_name in TARGET_FILES:
        found_files = glob.glob(os.path.join(work_dir, "**", target_name), recursive=True)
        for fpath in found_files:
            task_dir = os.path.dirname(fpath)
            # Check exit code if available (skip failed tasks)
            exit_file = os.path.join(task_dir, ".exitcode")
            if os.path.exists(exit_file):
                try:
                    with open(exit_file, "r") as ef:
                        code = ef.read().strip()
                    if code != "0":
                        continue
                except Exception:
                    pass

            rep_key = extract_rep_key_from_task_dir(task_dir)
            if rep_key:
                tasks_by_rep[rep_key][target_name] = fpath

    # Also collect individual evaluation CSVs into their respective replication folders
    csv_files = glob.glob(os.path.join(work_dir, "**", "*.csv"), recursive=True)
    for cf in csv_files:
        bname = os.path.basename(cf)
        if bname.startswith(".") or "benchmark" in bname:
            continue
        task_dir = os.path.dirname(cf)
        rep_key = extract_rep_key_from_task_dir(task_dir)
        if rep_key:
            tasks_by_rep[rep_key][bname] = cf

    return tasks_by_rep

def organize_replications(work_dir, outdir, mode="hardlink", threads=16):
    """
    Organizes replication files into results/replications/ using hardlinks or copies.
    """
    tasks_by_rep = harvest_work_artifacts(work_dir)
    if not tasks_by_rep:
        print(f"Warning: No valid replication artifacts found in '{work_dir}'.", file=sys.stderr)
        return

    os.makedirs(outdir, exist_ok=True)
    transfer_tasks = []

    for rep_key, file_dict in sorted(tasks_by_rep.items()):
        rep_dir = os.path.join(outdir, rep_key)
        for target_name, src_file in file_dict.items():
            dst_file = os.path.join(rep_dir, target_name)
            transfer_tasks.append((src_file, dst_file))

    print(f"Organizing {len(transfer_tasks)} files across {len(tasks_by_rep)} replications into '{outdir}' (mode: {mode})...")

    def _transfer(pair):
        src, dst = pair
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        if os.path.exists(dst):
            try:
                os.remove(dst)
            except OSError:
                pass
        if mode == "hardlink":
            try:
                os.link(src, dst)
                return
            except OSError:
                # Fallback to copy if cross-device link fails
                pass
        shutil.copy2(src, dst)

    with ThreadPoolExecutor(max_workers=threads) as executor:
        list(executor.map(_transfer, transfer_tasks))

    print(f"Successfully organized {len(transfer_tasks)} replication artifacts into '{outdir}'!")

def main():
    parser = argparse.ArgumentParser(
        description="Harvest and organize all replication files from work/ into results/replications/"
    )
    parser.add_argument(
        "--workdir",
        default="work",
        help="Path to Nextflow work/ directory (default: work)"
    )
    parser.add_argument(
        "--outdir",
        default="results/replications",
        help="Destination directory for replication artifacts (default: results/replications)"
    )
    parser.add_argument(
        "--mode",
        choices=["hardlink", "copy"],
        default="hardlink",
        help="Transfer mode: 'hardlink' (instantaneous, saves disk space) or 'copy' (default: hardlink)"
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=16,
        help="Number of concurrent worker threads (default: 16)"
    )
    args = parser.parse_args()

    if not os.path.exists(args.workdir):
        # Auto-detect if work/ exists in parent or current directory
        if os.path.exists("../work"):
            args.workdir = "../work"
        else:
            print(f"Error: Work directory '{args.workdir}' not found.", file=sys.stderr)
            sys.exit(1)

    organize_replications(
        work_dir=args.workdir,
        outdir=args.outdir,
        mode=args.mode,
        threads=args.threads
    )

if __name__ == "__main__":
    main()
