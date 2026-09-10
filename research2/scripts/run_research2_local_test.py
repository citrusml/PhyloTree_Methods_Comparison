#!/usr/bin/env python3
"""
run_research2_local_test.py

Runs a lightweight local test (4 conditions x 1 replicate x distance 0.5 and 2.0 x length 500)
without requiring Nextflow on local machine. Verifies data generation, tree estimation,
MSA property calculation, and plotting end-to-end.
"""

import os
import sys
import shutil
import subprocess
from pathlib import Path
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
BIN_DIR = PROJECT_ROOT / "bin"
RESEARCH2_DIR = PROJECT_ROOT / "research2"
TEST_WORK_DIR = RESEARCH2_DIR / "test_run"
RESULTS_DIR = RESEARCH2_DIR / "results"
PLOTS_DIR = RESEARCH2_DIR / "plots"

TEST_WORK_DIR.mkdir(parents=True, exist_ok=True)
RESULTS_DIR.mkdir(parents=True, exist_ok=True)
PLOTS_DIR.mkdir(parents=True, exist_ok=True)

CONDITIONS = [
    ("indelible_rate0.10", "indelible", 0.10, 0.10),
    ("indelible_rate0.05", "indelible", 0.05, 0.05),
    ("alisim_rate0.10",    "alisim",    0.10, 0.10),
    ("alisim_rate0.05",    "alisim",    0.05, 0.05),
]

DISTANCES = [0.5, 2.0]
LENGTHS = [500]
REP = 1
TAXA = 20
MODEL = "LG+G4"
ALPHA = 1.0
INDEL_SIZE = "POW 1.7 50"

print("=" * 60)
print("Starting research2 Local End-to-End Test")
print("=" * 60)

csv_outputs = []

for cond_name, simulator, ins_rate, del_rate in CONDITIONS:
    for length in LENGTHS:
        for dist in DISTANCES:
            run_key = f"{cond_name}_D{dist}_L{length}_rep{REP}"
            run_dir = TEST_WORK_DIR / run_key
            run_dir.mkdir(parents=True, exist_ok=True)
            print(f"\n--- Running: {run_key} ---")

            # 1. Generate Tree
            tree_file = run_dir / f"true_tree_{REP}.nwk"
            cmd_tree = [
                sys.executable, str(BIN_DIR / "generate_tree.py"),
                "--taxa", str(TAXA),
                "--scale", str(dist),
                "--seed", str(100 + int(dist * 10)),
                "--model", "paper_yule",
                "--rate_sd", "0.0",
                "--lba_ratio", "1.0",
                "--outtree", str(tree_file)
            ]
            subprocess.run(cmd_tree, check=True, stdout=subprocess.DEVNULL)

            # 2. Simulate Data
            seqs_file = run_dir / f"seqs_{REP}.fasta"
            true_msa_file = run_dir / f"true_msa_{REP}.fasta"

            if simulator == "indelible":
                with open(tree_file) as f:
                    tree_str = f.read().strip()
                if not tree_str.endswith(";"):
                    tree_str += ";"
                ctrl_content = f"""[TYPE] AMINOACID 1
[SETTINGS]
  [randomseed] {100 + int(dist * 10)}
  [output] FASTA
[MODEL] mymodel
  [submodel] LG
  [rates] 0 {ALPHA} 4
  [insertmodel] {INDEL_SIZE}
  [deletemodel] {INDEL_SIZE}
  [insertrate] {ins_rate}
  [deleterate] {del_rate}
[TREE] mytree {tree_str}
[PARTITIONS] mypart [mytree mymodel {length}]
[EVOLVE] mypart 1 sim_{REP}
"""
                with open(run_dir / "control.txt", "w") as f:
                    f.write(ctrl_content)
                subprocess.run(
                    [str(BIN_DIR / "indelible"), "control.txt"],
                    cwd=run_dir, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
                )
                shutil.move(run_dir / f"sim_{REP}.fas", seqs_file)
                shutil.move(run_dir / f"sim_{REP}_TRUE.fas", true_msa_file)
            else:
                iqtree_bin = BIN_DIR / "iqtree2"
                cmd_ali = [
                    str(iqtree_bin), "--alisim", str(run_dir / f"sim_{REP}"),
                    "-m", MODEL,
                    "--length", str(length),
                    "-t", str(tree_file),
                    "--indel", f"{ins_rate},{del_rate}",
                    "--indel-size", "POW{1.7/50},POW{1.7/50}",
                    "-af", "fasta",
                    "-seed", str(100 + int(dist * 10)),
                    "--redo"
                ]
                subprocess.run(cmd_ali, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                shutil.move(run_dir / f"sim_{REP}.unaligned.fa", seqs_file)
                shutil.move(run_dir / f"sim_{REP}.fa", true_msa_file)

            # 3. MAFFT
            mafft_msa = run_dir / f"mafft_aln_{REP}.fasta"
            mafft_bin = shutil.which("mafft") or "/opt/homebrew/bin/mafft"
            subprocess.run(
                f"{mafft_bin} --auto --quiet {seqs_file} > {mafft_msa}",
                shell=True, check=True
            )

            # 4. PWA+NJ
            pwa_nj_tree = run_dir / f"tree_pwa_nj_{REP}.nwk"
            cmd_pwa = [
                sys.executable, str(BIN_DIR / "run_pwa_nj.py"),
                "--fasta", str(seqs_file),
                "--gap_open", "20.0",
                "--gap_extend", "2.0",
                "--dist_model", "poisson",
                "--tool", "rapidnj",
                "--outtree", str(pwa_nj_tree)
            ]
            subprocess.run(cmd_pwa, check=True, stdout=subprocess.DEVNULL)

            # 5. MSA+NJ
            msa_nj_tree = run_dir / f"tree_msa_nj_{REP}.nwk"
            cmd_msa_nj = [
                sys.executable, str(BIN_DIR / "run_msa_nj.py"),
                "--msa", str(mafft_msa),
                "--dist_model", "poisson",
                "--tool", "rapidnj",
                "--outtree", str(msa_nj_tree)
            ]
            subprocess.run(cmd_msa_nj, check=True, stdout=subprocess.DEVNULL)

            # 6. MSA+ML (IQ-TREE 2 fast mode for testing)
            msa_ml_tree = run_dir / f"tree_msa_ml_{REP}.nwk"
            cmd_ml = [
                str(BIN_DIR / "iqtree2"),
                "-s", str(mafft_msa),
                "-m", MODEL,
                "--prefix", str(run_dir / "iqtree_msa_ml"),
                "-T", "2",
                "--quiet", "-redo"
            ]
            subprocess.run(cmd_ml, check=True)
            shutil.move(run_dir / "iqtree_msa_ml.treefile", msa_ml_tree)

            # 7. TRUE_MSA+NJ
            true_nj_tree = run_dir / f"tree_true_msa_nj_{REP}.nwk"
            cmd_t_nj = [
                sys.executable, str(BIN_DIR / "run_msa_nj.py"),
                "--msa", str(true_msa_file),
                "--dist_model", "poisson",
                "--tool", "rapidnj",
                "--outtree", str(true_nj_tree)
            ]
            subprocess.run(cmd_t_nj, check=True, stdout=subprocess.DEVNULL)

            # 8. TRUE_MSA+ML
            true_ml_tree = run_dir / f"tree_true_msa_ml_{REP}.nwk"
            cmd_t_ml = [
                str(BIN_DIR / "iqtree2"),
                "-s", str(true_msa_file),
                "-m", MODEL,
                "--prefix", str(run_dir / "iqtree_true_ml"),
                "-T", "2",
                "--quiet", "-redo"
            ]
            subprocess.run(cmd_t_ml, check=True)
            shutil.move(run_dir / "iqtree_true_ml.treefile", true_ml_tree)

            # 9. SSS Calculation
            sss_csv = run_dir / f"sss_{REP}.csv"
            cmd_sss = [
                sys.executable, str(BIN_DIR / "calculate_sss.py"),
                "--fasta", str(seqs_file),
                "--distance", str(dist),
                "--length", str(length),
                "--rep_start", str(REP),
                "--rep_end", str(REP),
                "--outcsv", str(sss_csv)
            ]
            subprocess.run(cmd_sss, check=True, stdout=subprocess.DEVNULL)

            # 10. Calculate Comprehensive Properties
            out_props_csv = run_dir / f"props_{run_key}.csv"
            cmd_calc = [
                sys.executable, str(RESEARCH2_DIR / "scripts" / "calculate_msa_properties.py"),
                "--cond_name", cond_name,
                "--simulator", simulator,
                "--ins_rate", str(ins_rate),
                "--del_rate", str(del_rate),
                "--dist", str(dist),
                "--length", str(length),
                "--rep", str(REP),
                "--true_tree", str(tree_file),
                "--true_msa", str(true_msa_file),
                "--unaligned", str(seqs_file),
                "--mafft_msa", str(mafft_msa),
                "--tree_pwa_nj", str(pwa_nj_tree),
                "--tree_msa_nj", str(msa_nj_tree),
                "--tree_msa_ml", str(msa_ml_tree),
                "--tree_true_msa_nj", str(true_nj_tree),
                "--tree_true_msa_ml", str(true_ml_tree),
                "--sss_csv", str(sss_csv),
                "--outfile", str(out_props_csv)
            ]
            subprocess.run(cmd_calc, check=True)
            csv_outputs.append(out_props_csv)

# 11. Combine CSVs and Plot
print("\n" + "=" * 60)
print("Aggregating Results and Generating Plots")
print("=" * 60)
dfs = [pd.read_csv(f) for f in csv_outputs]
all_df = pd.concat(dfs, ignore_index=True)
combined_csv = RESULTS_DIR / "all_properties_research2_test.csv"
all_df.to_csv(combined_csv, index=False)
print(f"Saved combined CSV: {combined_csv} ({len(all_df)} rows)")

# Run Plotter
cmd_plot = [
    sys.executable, str(RESEARCH2_DIR / "scripts" / "plot_research2_comparison.py"),
    "--csv", str(combined_csv),
    "--outdir", str(PLOTS_DIR)
]
subprocess.run(cmd_plot, check=True)

print("\nSUCCESS! Test completed successfully.")
print(f"Check results at: {RESULTS_DIR}")
print(f"Check plots at: {PLOTS_DIR}")
