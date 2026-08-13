#!/usr/bin/env python3
"""
Run Local Phase 0 Benchmark Test
Runs a quick local Phase 0 validation test over mini parameter grid
(D in [0.1, 1.0], L in [100, 500], replicates=3)
verifying end-to-end functionality of PWA+NJ, MSA+NJ, and MSA+ML pipelines.
"""

import sys
import os
import subprocess
import pandas as pd

def run_cmd(cmd):
    print(f"Running: {' '.join(cmd)}")
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        print(f"Error executing command: {res.stderr}")
    return res.returncode == 0

def main():
    print("=== Phase 0 Local Pipeline Benchmark Test ===")
    outdir = "results_phase0"
    os.makedirs(outdir, exist_ok=True)
    csv_path = os.path.join(outdir, "benchmark_summary.csv")
    if os.path.exists(csv_path):
        os.remove(csv_path)

    distances = [0.1, 1.0]
    lengths = [100, 500]
    replicates = [1, 2, 3]

    project_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    bin_dir = os.path.join(project_dir, "bin")

    for dist in distances:
        for length in lengths:
            for rep in replicates:
                tag = f"D{dist}_L{length}_rep{rep}"
                true_tree = os.path.join(outdir, f"{tag}_true.nwk")
                fasta = os.path.join(outdir, f"{tag}_seqs.fasta")

                # 1. Simulate Data
                sim_cmd = [
                    sys.executable, os.path.join(bin_dir, "simulate_data.py"),
                    "--distance", str(dist),
                    "--length", str(length),
                    "--num_taxa", "16",
                    "--outtree", true_tree,
                    "--outfasta", fasta,
                    "--seed", str(rep)
                ]
                if not run_cmd(sim_cmd):
                    continue

                # 2. PWA + NJ
                pwa_tree = os.path.join(outdir, f"{tag}_pwa_nj.nwk")
                run_cmd([
                    sys.executable, os.path.join(bin_dir, "run_pwa_nj.py"),
                    "--fasta", fasta,
                    "--outtree", pwa_tree
                ])
                run_cmd([
                    sys.executable, os.path.join(bin_dir, "evaluate_trees.py"),
                    "--truetree", true_tree,
                    "--esttree", pwa_tree,
                    "--pipeline", "PWA+NJ",
                    "--distance", str(dist),
                    "--length", str(length),
                    "--replicate", str(rep),
                    "--outcsv", csv_path
                ])

                # 3. MSA + NJ
                msa_nj_tree = os.path.join(outdir, f"{tag}_msa_nj.nwk")
                run_cmd([
                    sys.executable, os.path.join(bin_dir, "run_msa_nj.py"),
                    "--fasta", fasta,
                    "--outtree", msa_nj_tree
                ])
                run_cmd([
                    sys.executable, os.path.join(bin_dir, "evaluate_trees.py"),
                    "--truetree", true_tree,
                    "--esttree", msa_nj_tree,
                    "--pipeline", "MSA+NJ",
                    "--distance", str(dist),
                    "--length", str(length),
                    "--replicate", str(rep),
                    "--outcsv", csv_path
                ])

                # 4. MSA + ML
                msa_ml_tree = os.path.join(outdir, f"{tag}_msa_ml.nwk")
                json_meta = os.path.join(outdir, f"{tag}_msa_ml.json")
                run_cmd([
                    sys.executable, os.path.join(bin_dir, "run_msa_ml.py"),
                    "--fasta", fasta,
                    "--outtree", msa_ml_tree,
                    "--outjson", json_meta
                ])
                run_cmd([
                    sys.executable, os.path.join(bin_dir, "evaluate_trees.py"),
                    "--truetree", true_tree,
                    "--esttree", msa_ml_tree,
                    "--pipeline", "MSA+ML",
                    "--distance", str(dist),
                    "--length", str(length),
                    "--replicate", str(rep),
                    "--json", json_meta,
                    "--outcsv", csv_path
                ])

    # 5. Plot Regime Map
    plot_cmd = [
        sys.executable, os.path.join(bin_dir, "plot_regime_map.py"),
        "--csv", csv_path,
        "--outdir", outdir
    ]
    run_cmd(plot_cmd)
    print("=== Phase 0 Local Test Completed Successfully! ===")

if __name__ == "__main__":
    main()
