#!/usr/bin/env python3
"""
run_dataset_methods.py
======================
Runs the 4 phylogenetic reconstruction methods on the 400 real biological HOG datasets
(Amniota, Eukaryota, Bacteria, Actinobacteria) located in dataset/data/:
  1. MSA + ML (MAFFT + IQ-TREE 2 with ModelFinder -m MFP)
  2. MSA + NJ (MAFFT + Poisson distance + RapidNJ)
  3. PSA + NJ (Needleman-Wunsch BLOSUM62 gap_open=20, gap_extend=2 + Poisson distance + RapidNJ)
  4. GS       (Graph Splitting: gs2 -s -l -t 1 -m 7.5)

The execution options are strictly identical to Experiment 15 (Paper Tree Benchmark).
Outputs:
  - Trees: dataset/results/trees/<method>/<group>_<hog_id>.nwk
  - Execution summary: dataset/results/dataset_trees_summary.csv
  - Method concordance: dataset/results/method_concordance_summary.csv
"""

import os
import sys
import time
import shutil
import argparse
import subprocess
from pathlib import Path
from typing import Dict, List, Tuple, Any, Optional
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import pandas as pd
from Bio import SeqIO
import dendropy
from dendropy.calculate import treecompare

# Add project root to sys.path
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))
sys.path.insert(0, str(PROJECT_ROOT / "bin"))

try:
    from run_pwa_nj import calculate_distance, needleman_wunsch, run_nj_tool
except ImportError:
    # Direct import fallback if needed
    import importlib.util
    spec = importlib.util.spec_from_file_location("run_pwa_nj", str(PROJECT_ROOT / "bin" / "run_pwa_nj.py"))
    run_pwa_nj_mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(run_pwa_nj_mod)
    calculate_distance = run_pwa_nj_mod.calculate_distance
    needleman_wunsch = run_pwa_nj_mod.needleman_wunsch
    run_nj_tool = run_pwa_nj_mod.run_nj_tool

# Map short dataset keys to directory names
DATASET_DIRS = {
    "amn": "amniota",
    "euk": "eukaryota",
    "bac": "bacteria",
    "act": "actinomycete"
}

DATASET_NAMES = {
    "amn": "Amniota",
    "euk": "Eukaryota",
    "bac": "Bacteria",
    "act": "Actinobacteria"
}


def find_iqtree_binary() -> str:
    """Finds available IQ-TREE executable."""
    py_dir = os.path.dirname(sys.executable)
    for bin_name in ["iqtree3", "iqtree2", "iqtree"]:
        full_path = os.path.join(py_dir, bin_name)
        if os.path.isfile(full_path) and os.access(full_path, os.X_OK):
            return full_path
        which_path = shutil.which(bin_name)
        if which_path:
            return which_path
    raise RuntimeError("IQ-TREE executable ('iqtree3', 'iqtree2', or 'iqtree') not found.")


def find_gs2_binary() -> str:
    """Finds gs2 executable."""
    candidate_paths = [
        PROJECT_ROOT / "bin" / "gs2",
        Path(shutil.which("gs2") or "")
    ]
    for p in candidate_paths:
        if p and p.is_file() and os.access(p, os.X_OK):
            return str(p)
    raise RuntimeError("gs2 binary not found. Please build it via bin/install_gs.sh.")


def sanitize_fasta_headers(input_fasta: Path, output_fasta: Path) -> List[str]:
    """
    Sanitizes FASTA headers so that each sequence has only its clean ID (e.g. 'MOUSE20754')
    without long descriptions or special characters that confuse Newick parsers.
    Returns the list of sequence IDs.
    """
    records = list(SeqIO.parse(input_fasta, "fasta"))
    if not records:
        raise ValueError(f"No sequences found in {input_fasta}")

    clean_ids = []
    with open(output_fasta, "w") as out_f:
        for rec in records:
            clean_id = rec.id.strip().split()[0]
            clean_ids.append(clean_id)
            out_f.write(f">{clean_id}\n{str(rec.seq)}\n")
    return clean_ids


def run_mafft(in_fasta: Path, out_msa: Path) -> None:
    """Runs MAFFT with --threadit 0 --auto."""
    cmd = ["mafft", "--threadit", "0", "--auto", str(in_fasta)]
    with open(out_msa, "w") as out_f:
        res = subprocess.run(cmd, stdout=out_f, stderr=subprocess.PIPE, text=True)
    if res.returncode != 0 or not out_msa.exists() or out_msa.stat().st_size == 0:
        raise RuntimeError(f"MAFFT failed (exit code {res.returncode}):\n{res.stderr}")


def run_method_msa_ml(msa_file: Path, out_tree: Path, work_dir: Path, iqtree_bin: str) -> str:
    """
    Executes MSA + ML (IQ-TREE 2 with ModelFinder -m MFP, no bootstrap, -T 1).
    Returns the best-fit model according to BIC.
    """
    prefix = work_dir / "iqtree_run"
    cmd = [
        iqtree_bin,
        "-s", str(msa_file),
        "-m", "MFP",
        "--prefix", str(prefix),
        "--redo",
        "-T", "1"
    ]
    res = subprocess.run(cmd, capture_output=True, text=True)
    tree_file = Path(str(prefix) + ".treefile")
    iq_file = Path(str(prefix) + ".iqtree")

    if not tree_file.exists() or tree_file.stat().st_size == 0:
        raise RuntimeError(f"IQ-TREE failed:\nSTDOUT: {res.stdout}\nSTDERR: {res.stderr}")

    shutil.copyfile(tree_file, out_tree)

    best_model = "Unknown"
    if iq_file.exists():
        with open(iq_file, "r") as f:
            for line in f:
                if "Best-fit model according to BIC:" in line:
                    best_model = line.split(":")[-1].strip()
                    break
    return best_model


def run_method_msa_nj(msa_file: Path, out_tree: Path, work_dir: Path) -> None:
    """
    Executes MSA + NJ (MAFFT alignment + Poisson distance + RapidNJ).
    Matches Experiment 15 options.
    """
    records = list(SeqIO.parse(msa_file, "fasta"))
    names = [rec.id for rec in records]
    seqs = [str(rec.seq) for rec in records]
    N = len(names)

    dist_matrix = [[0.0] * N for _ in range(N)]
    for i in range(N):
        for j in range(i + 1, N):
            d = calculate_distance(seqs[i], seqs[j], dist_model="poisson", alpha=1.0)
            dist_matrix[i][j] = d
            dist_matrix[j][i] = d

    matrix_file = work_dir / "msa_matrix.phylip"
    matrix_str = f"   {N}\n"
    for i in range(N):
        row = f"{names[i]:<10}" + "".join(f"  {dist_matrix[i][j]:.6f}" for j in range(N)) + "\n"
        matrix_str += row

    with open(matrix_file, "w") as f:
        f.write(matrix_str)

    run_nj_tool(str(matrix_file), str(out_tree), tool="rapidnj")


def run_method_psa_nj(clean_fasta: Path, out_tree: Path, work_dir: Path) -> None:
    """
    Executes PSA + NJ (Pairwise Needleman-Wunsch BLOSUM62 + Poisson distance + RapidNJ).
    Options strictly follow Experiment 15: Gap Open = 20.0, Gap Extension = 2.0.
    """
    records = list(SeqIO.parse(clean_fasta, "fasta"))
    names = [rec.id for rec in records]
    seqs = [str(rec.seq) for rec in records]
    N = len(names)

    dist_matrix = [[0.0] * N for _ in range(N)]
    for i in range(N):
        for j in range(i + 1, N):
            al1, al2 = needleman_wunsch(seqs[i], seqs[j], gap_open=20.0, gap_extend=2.0)
            d = calculate_distance(al1, al2, dist_model="poisson", alpha=1.0)
            dist_matrix[i][j] = d
            dist_matrix[j][i] = d

    matrix_file = work_dir / "psa_matrix.phylip"
    matrix_str = f"   {N}\n"
    for i in range(N):
        row = f"{names[i]:<10}" + "".join(f"  {dist_matrix[i][j]:.6f}" for j in range(N)) + "\n"
        matrix_str += row

    with open(matrix_file, "w") as f:
        f.write(matrix_str)

    run_nj_tool(str(matrix_file), str(out_tree), tool="rapidnj")


def run_method_gs(clean_fasta: Path, out_tree: Path, work_dir: Path, gs2_bin: str) -> None:
    """
    Executes Graph Splitting: gs2 -s -l -t 1 -m 7.5 <fasta> | tr -d '"' > <outtree>.
    Options strictly follow Experiment 15: sensitivity = 7.5.
    Runs inside work_dir to isolate temporary files across parallel workers.
    """
    cmd = f"{gs2_bin} -s -l -t 1 -m 7.5 '{clean_fasta.resolve()}' | tr -d '\"' > '{out_tree.resolve()}'"
    res = subprocess.run(cmd, shell=True, cwd=str(work_dir), capture_output=True, text=True)
    if res.returncode != 0 or not out_tree.exists() or out_tree.stat().st_size == 0:
        raise RuntimeError(f"GS execution failed (exit code {res.returncode}):\n{res.stderr}\n{res.stdout}")


def compute_pairwise_rf_metrics(trees_dict: Dict[str, Path]) -> Dict[str, float]:
    """
    Calculates pairwise Robinson-Foulds (RF and normalized nRF) topological distances
    between the inferred trees for all method pairs using DendroPy.
    """
    methods = ["msa_ml", "msa_nj", "psa_nj", "gs"]
    metrics = {}

    # Initialise default values
    for i in range(len(methods)):
        for j in range(i + 1, len(methods)):
            m1, m2 = methods[i], methods[j]
            metrics[f"rf_{m1}_vs_{m2}"] = np.nan
            metrics[f"nrf_{m1}_vs_{m2}"] = np.nan

    loaded_trees = {}
    tns = dendropy.TaxonNamespace()

    for m in methods:
        p = trees_dict.get(m)
        if p and p.exists() and p.stat().st_size > 0:
            try:
                tree = dendropy.Tree.get(path=str(p), schema="newick", taxon_namespace=tns, preserve_underscores=True)
                tree.is_rooted = False
                tree.deroot()
                tree.encode_bipartitions()
                loaded_trees[m] = tree
            except Exception as e:
                # Tree parsing error
                pass

    num_taxa = len(tns)
    max_rf = 2 * (num_taxa - 3) if num_taxa > 3 else 1

    for i in range(len(methods)):
        for j in range(i + 1, len(methods)):
            m1, m2 = methods[i], methods[j]
            if m1 in loaded_trees and m2 in loaded_trees:
                t1, t2 = loaded_trees[m1], loaded_trees[m2]
                rf = float(treecompare.symmetric_difference(t1, t2))
                nrf = float(rf / max_rf if max_rf > 0 else 0.0)
                metrics[f"rf_{m1}_vs_{m2}"] = rf
                metrics[f"nrf_{m1}_vs_{m2}"] = nrf

    return metrics


def process_single_hog(task_args: Tuple[str, str, Path, Path, Path, str, str, bool]) -> Dict[str, Any]:
    """
    Processes one HOG:
      - Sanitizes FASTA headers
      - Runs MSA+ML, MSA+NJ, PSA+NJ, GS
      - Calculates method concordance (RF and nRF)
      - Returns execution summary
    """
    dataset_key, hog_id, fasta_file, out_trees_dir, scratch_dir, iqtree_bin, gs2_bin, resume = task_args

    hog_scratch = scratch_dir / f"{dataset_key}_{hog_id}"
    hog_scratch.mkdir(parents=True, exist_ok=True)

    tree_files = {
        "msa_ml": out_trees_dir / "msa_ml" / f"{dataset_key}_{hog_id}.nwk",
        "msa_nj": out_trees_dir / "msa_nj" / f"{dataset_key}_{hog_id}.nwk",
        "psa_nj": out_trees_dir / "psa_nj" / f"{dataset_key}_{hog_id}.nwk",
        "gs": out_trees_dir / "gs" / f"{dataset_key}_{hog_id}.nwk"
    }

    status = {
        "msa_ml": "PENDING",
        "msa_nj": "PENDING",
        "psa_nj": "PENDING",
        "gs": "PENDING"
    }
    elapsed_times = {
        "msa_ml": 0.0,
        "msa_nj": 0.0,
        "psa_nj": 0.0,
        "gs": 0.0
    }
    best_model = "None"

    t_start = time.time()

    try:
        clean_fasta = hog_scratch / f"{hog_id}_clean.fa"
        clean_taxa = sanitize_fasta_headers(fasta_file, clean_fasta)
        num_taxa = len(clean_taxa)

        # 1. MSA generation (shared by MSA+ML and MSA+NJ)
        clean_msa = hog_scratch / f"{hog_id}_clean.msa"
        msa_done = False
        if not (resume and tree_files["msa_ml"].exists() and tree_files["msa_nj"].exists()):
            run_mafft(clean_fasta, clean_msa)
            msa_done = True

        # --- Method 1: MSA + ML ---
        if resume and tree_files["msa_ml"].exists() and tree_files["msa_ml"].stat().st_size > 0:
            status["msa_ml"] = "SKIPPED_EXISTING"
        else:
            t0 = time.time()
            try:
                best_model = run_method_msa_ml(clean_msa, tree_files["msa_ml"], hog_scratch, iqtree_bin)
                status["msa_ml"] = "SUCCESS"
            except Exception as e:
                status["msa_ml"] = f"ERROR: {e}"
            elapsed_times["msa_ml"] = round(time.time() - t0, 3)

        # --- Method 2: MSA + NJ ---
        if resume and tree_files["msa_nj"].exists() and tree_files["msa_nj"].stat().st_size > 0:
            status["msa_nj"] = "SKIPPED_EXISTING"
        else:
            t0 = time.time()
            try:
                if not clean_msa.exists():
                    run_mafft(clean_fasta, clean_msa)
                run_method_msa_nj(clean_msa, tree_files["msa_nj"], hog_scratch)
                status["msa_nj"] = "SUCCESS"
            except Exception as e:
                status["msa_nj"] = f"ERROR: {e}"
            elapsed_times["msa_nj"] = round(time.time() - t0, 3)

        # --- Method 3: PSA + NJ ---
        if resume and tree_files["psa_nj"].exists() and tree_files["psa_nj"].stat().st_size > 0:
            status["psa_nj"] = "SKIPPED_EXISTING"
        else:
            t0 = time.time()
            try:
                run_method_psa_nj(clean_fasta, tree_files["psa_nj"], hog_scratch)
                status["psa_nj"] = "SUCCESS"
            except Exception as e:
                status["psa_nj"] = f"ERROR: {e}"
            elapsed_times["psa_nj"] = round(time.time() - t0, 3)

        # --- Method 4: GS ---
        if resume and tree_files["gs"].exists() and tree_files["gs"].stat().st_size > 0:
            status["gs"] = "SKIPPED_EXISTING"
        else:
            t0 = time.time()
            try:
                run_method_gs(clean_fasta, tree_files["gs"], hog_scratch, gs2_bin)
                status["gs"] = "SUCCESS"
            except Exception as e:
                status["gs"] = f"ERROR: {e}"
            elapsed_times["gs"] = round(time.time() - t0, 3)

        # Compute topological concordance (RF / nRF)
        rf_metrics = compute_pairwise_rf_metrics(tree_files)

    finally:
        # Clean up temporary files in scratch
        if hog_scratch.exists():
            try:
                shutil.rmtree(hog_scratch)
            except OSError:
                pass

    total_time = round(time.time() - t_start, 3)

    summary_row = {
        "dataset": dataset_key,
        "dataset_name": DATASET_NAMES[dataset_key],
        "hog_id": hog_id,
        "num_taxa": num_taxa,
        "msa_ml_status": status["msa_ml"],
        "msa_ml_time_s": elapsed_times["msa_ml"],
        "msa_ml_model": best_model,
        "msa_nj_status": status["msa_nj"],
        "msa_nj_time_s": elapsed_times["msa_nj"],
        "psa_nj_status": status["psa_nj"],
        "psa_nj_time_s": elapsed_times["psa_nj"],
        "gs_status": status["gs"],
        "gs_time_s": elapsed_times["gs"],
        "total_time_s": total_time
    }

    concordance_row = {
        "dataset": dataset_key,
        "dataset_name": DATASET_NAMES[dataset_key],
        "hog_id": hog_id,
        "num_taxa": num_taxa,
        **rf_metrics
    }

    return {
        "summary": summary_row,
        "concordance": concordance_row
    }


def main():
    parser = argparse.ArgumentParser(description="Run 4 Phylogenetic Methods on dataset/data HOGs")
    parser.add_argument("--data_dir", type=Path, default=PROJECT_ROOT / "dataset" / "data", help="Path to dataset/data directory")
    parser.add_argument("--results_dir", type=Path, default=PROJECT_ROOT / "dataset" / "results", help="Path to dataset/results directory")
    parser.add_argument("--scratch_dir", type=Path, default=PROJECT_ROOT / "scratch" / "tmp_dataset_methods", help="Temporary working directory")
    parser.add_argument("--workers", type=int, default=6, help="Number of parallel worker processes (default: 6)")
    parser.add_argument("--groups", nargs="+", default=["amn", "euk", "bac", "act"], help="Datasets to run (default: amn euk bac act)")
    parser.add_argument("--single_hog", type=str, default=None, help="Run only a single HOG ID for testing")
    parser.add_argument("--force", action="store_true", help="Force recomputation of existing trees")
    args = parser.parse_args()

    iqtree_bin = find_iqtree_binary()
    gs2_bin = find_gs2_binary()
    resume = not args.force

    out_trees_dir = args.results_dir / "trees"
    for m in ["msa_ml", "msa_nj", "psa_nj", "gs"]:
        (out_trees_dir / m).mkdir(parents=True, exist_ok=True)
    args.scratch_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 80)
    print("  Real Dataset Phylogenetic Methods Benchmark (Experiment 15 Options)")
    print(f"  Data Directory    : {args.data_dir}")
    print(f"  Results Directory : {args.results_dir}")
    print(f"  Tree Output Dir   : {out_trees_dir}")
    print(f"  Parallel Workers  : {args.workers}")
    print(f"  Target Groups     : {args.groups}")
    print(f"  IQ-TREE Binary    : {iqtree_bin}")
    print(f"  GS2 Binary        : {gs2_bin}")
    print(f"  Resume Mode       : {resume}")
    print("=" * 80 + "\n")

    # Collect tasks
    tasks = []
    for d_key in args.groups:
        d_dir_name = DATASET_DIRS[d_key]
        group_dir = args.data_dir / d_dir_name
        if not group_dir.exists():
            print(f"Warning: Directory {group_dir} not found. Skipping {d_key}.")
            continue

        fasta_files = sorted(list(group_dir.glob("*.fa")) + list(group_dir.glob("*.fasta")))
        for fpath in fasta_files:
            raw_id = fpath.stem
            hog_id = raw_id.replace("HOG_", "").replace("HOG:", "")
            if args.single_hog and hog_id != args.single_hog and raw_id != args.single_hog:
                continue
            tasks.append((d_key, hog_id, fpath, out_trees_dir, args.scratch_dir, iqtree_bin, gs2_bin, resume))

    print(f"Total HOG tasks queued: {len(tasks)}")
    if not tasks:
        print("No tasks to run. Exiting.")
        return

    all_summaries = []
    all_concordance = []

    t_global_start = time.time()
    completed_count = 0

    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        future_to_task = {executor.submit(process_single_hog, t): t for t in tasks}
        for future in as_completed(future_to_task):
            t_info = future_to_task[future]
            try:
                res = future.result()
                all_summaries.append(res["summary"])
                all_concordance.append(res["concordance"])
                completed_count += 1
                if completed_count % 10 == 0 or completed_count == len(tasks):
                    elapsed = time.time() - t_global_start
                    pct = (completed_count / len(tasks)) * 100
                    print(f"  [{completed_count:3d}/{len(tasks):3d}] ({pct:5.1f}%) Completed | Elapsed: {elapsed:6.1f}s | Latest: {t_info[0]}_{t_info[1]}")
            except Exception as e:
                print(f"Error processing task {t_info[0]}_{t_info[1]}: {e}", file=sys.stderr)

    # Save summary CSVs
    summary_df = pd.DataFrame(all_summaries).sort_values(by=["dataset", "hog_id"])
    concordance_df = pd.DataFrame(all_concordance).sort_values(by=["dataset", "hog_id"])

    summary_csv = args.results_dir / "dataset_trees_summary.csv"
    concordance_csv = args.results_dir / "method_concordance_summary.csv"

    summary_df.to_csv(summary_csv, index=False)
    concordance_df.to_csv(concordance_csv, index=False)

    print("\n" + "=" * 80)
    print("  Benchmark Completed Successfully!")
    print(f"  Total Elapsed Time: {time.time() - t_global_start:.2f}s")
    print(f"  Saved Tree Summary      : {summary_csv} ({len(summary_df)} rows)")
    print(f"  Saved Method Concordance: {concordance_csv} ({len(concordance_df)} rows)")
    print("=" * 80)


if __name__ == "__main__":
    main()
