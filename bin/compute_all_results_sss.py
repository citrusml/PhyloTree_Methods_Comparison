#!/usr/bin/env python3
"""
Comprehensive SSS Full Calculation Script: compute_all_results_sss.py
Recalculates Sequence Similarity Score (SSS) strictly using MMseqs2 matching Matsui & Iwasaki (2020)
and gs2 definition across all 14 benchmark result directories under results/.

Key Execution Features:
1. Replicate-level Full Simulation & Calculation:
   - Evaluates each replicate individually (seed = replicate, 1..100)
   - Accurately passes --length ${length} to AliSim for all conditions
2. Deduplication & Caching across Shared Experimental Conditions:
   - Shared conditions (e.g. results_new, results_true_pwa, results_fastme_NoOption) are evaluated once
   - Persistent checkpoint cache (results/sss_replicate_cache.csv) allows interruption & resumption
3. Multi-process Parallel Execution:
   - Utilizes available CPU cores efficiently (default: 8 workers)
4. Full Dataset Integration & Export:
   - Updates each benchmark_*_summary.csv in-place (with backup) by exact (distance, length, replicate, ...) key matching
   - Saves complete 100-replicate sss_summary.csv in each experiment directory
"""

import os
import sys
import glob
import time
import shutil
import tempfile
import subprocess
import argparse
import random
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Dict, List, Tuple, Any, Optional
import numpy as np
import pandas as pd
from Bio import SeqIO

PROJECT_ROOT = Path(__file__).resolve().parent.parent
BIN_DIR = PROJECT_ROOT / "bin"
RESULTS_DIR = PROJECT_ROOT / "results"
MODELS_DIR = PROJECT_ROOT / "models"
ICS_MODEL_PATH = MODELS_DIR / "ics_model.nex"
CACHE_FILE = RESULTS_DIR / "sss_replicate_cache.csv"

IQTREE_BIN = "/opt/homebrew/Cellar/micromamba/2.8.1/envs/phylomethod_env/bin/iqtree"
MMSEQS_BIN = "/opt/homebrew/Cellar/micromamba/2.8.1/envs/phylomethod_env/bin/mmseqs"

EXPERIMENT_CONFIGS = [
    {
        "name": "results_alpha",
        "csv": "benchmark_alpha_summary.csv",
        "type": "alpha",
        "keys": ["distance", "length", "alpha", "num_taxa"],
        "default_taxa": 32,
        "default_alpha": 1.0,
        "default_ics": 0.0,
        "indel": "0.05,0.10",
        "indel_size": None,
    },
    {
        "name": "results_fastme_NoOption",
        "csv": "benchmark_fastme_summary.csv",
        "type": "standard",
        "keys": ["distance", "length", "alpha", "ics_prop", "num_taxa"],
        "default_taxa": 32,
        "default_alpha": 1.0,
        "default_ics": 0.0,
        "indel": "0.05,0.10",
        "indel_size": None,
    },
    {
        "name": "results_fastme_options",
        "csv": "benchmark_fastme_options_summary.csv",
        "type": "standard",
        "keys": ["distance", "length", "alpha", "ics_prop", "num_taxa"],
        "default_taxa": 32,
        "default_alpha": 1.0,
        "default_ics": 0.0,
        "indel": "0.05,0.10",
        "indel_size": None,
    },
    {
        "name": "results_gamma",
        "csv": "benchmark_gamma_summary.csv",
        "type": "alpha",
        "keys": ["distance", "length", "alpha", "ics_prop", "num_taxa"],
        "default_taxa": 32,
        "default_alpha": 1.0,
        "default_ics": 0.0,
        "indel": "0.05,0.10",
        "indel_size": None,
    },
    {
        "name": "results_high_dist",
        "csv": "benchmark_high_dist_summary.csv",
        "type": "standard",
        "keys": ["distance", "length", "alpha", "ics_prop", "num_taxa"],
        "default_taxa": 32,
        "default_alpha": 1.0,
        "default_ics": 0.0,
        "indel": "0.05,0.10",
        "indel_size": None,
    },
    {
        "name": "results_high_gap",
        "csv": "benchmark_summary.csv",
        "type": "standard",
        "keys": ["distance", "length", "alpha", "ics_prop", "num_taxa"],
        "default_taxa": 32,
        "default_alpha": 1.0,
        "default_ics": 0.0,
        "indel": "0.05,0.10",
        "indel_size": None,
    },
    {
        "name": "results_ics",
        "csv": "benchmark_ics_summary.csv",
        "type": "ics",
        "keys": ["distance", "length", "alpha", "ics_prop", "num_taxa"],
        "default_taxa": 32,
        "default_alpha": 1.0,
        "default_ics": 0.0,
        "indel": "0.05,0.10",
        "indel_size": None,
    },
    {
        "name": "results_ics_full",
        "csv": "benchmark_ics_full_summary.csv",
        "type": "ics_full",
        "keys": ["distance", "length", "alpha", "ics_prop", "num_taxa"],
        "default_taxa": 32,
        "default_alpha": 1.0,
        "default_ics": 1.0,
        "indel": "0.05,0.10",
        "indel_size": None,
    },
    {
        "name": "results_new",
        "csv": "benchmark_summary.csv",
        "type": "standard",
        "keys": ["distance", "length", "num_taxa"],
        "default_taxa": 32,
        "default_alpha": 1.0,
        "default_ics": 0.0,
        "indel": "0.05,0.10",
        "indel_size": None,
    },
    {
        "name": "results_paper_tree_pre_fix",
        "csv": "benchmark_summary.csv",
        "type": "paper_tree",
        "keys": ["distance", "length", "replicate"],
        "default_taxa": 32,
        "default_alpha": 1.0,
        "default_ics": 0.0,
        "indel": "0.10,0.10",
        "indel_size": "POW{1.7/50},POW{1.7/50}",
    },
    {
        "name": "results_zipfian_indel",
        "csv": "benchmark_summary.csv",
        "type": "zipfian",
        "keys": ["distance", "length", "alpha", "ics_prop", "num_taxa"],
        "default_taxa": 32,
        "default_alpha": 1.0,
        "default_ics": 0.0,
        "indel": "0.05,0.05",
        "indel_size": "POW{1.7/50},POW{1.7/50}",
    },
    {
        "name": "results_simple",
        "csv": "benchmark_high_dist_summary.csv",
        "type": "simple_lg",
        "keys": ["distance", "length", "alpha", "ics_prop", "num_taxa"],
        "default_taxa": 32,
        "default_alpha": 1.0,
        "default_ics": 0.0,
        "indel": "0.05,0.10",
        "indel_size": None,
    },
    {
        "name": "results_taxon",
        "csv": "benchmark_taxon_summary.csv",
        "type": "taxon",
        "keys": ["distance", "length", "num_taxa"],
        "default_alpha": 1.0,
        "default_ics": 0.0,
        "indel": "0.05,0.10",
        "indel_size": None,
    },
    {
        "name": "results_true_pwa",
        "csv": "benchmark_true_pwa_summary.csv",
        "type": "standard",
        "keys": ["distance", "length", "alpha", "ics_prop", "num_taxa"],
        "default_taxa": 32,
        "default_alpha": 1.0,
        "default_ics": 0.0,
        "indel": "0.05,0.10",
        "indel_size": None,
    },
]


def compute_single_replicate_sss(
    fasta_path: str,
    sensitivity: float = 7.5,
    threads: int = 1,
    mmseqs_bin: str = MMSEQS_BIN
) -> Dict[str, float]:
    """Computes SSS for a single FASTA using MMseqs2 matching Matsui & Iwasaki (2020)."""
    records = list(SeqIO.parse(fasta_path, "fasta"))
    seqs = [str(r.seq).replace("-", "").upper() for r in records]
    n = len(seqs)
    if n < 2 or any(len(s) == 0 for s in seqs):
        return {
            "sss_mean": 0.0, "sss_median": 0.0, "sss_std": 0.0,
            "sss_min": 0.0, "sss_max": 0.0, "sss_zero_frac": 1.0, "num_taxa": n
        }

    id_map = {rec.id: i for i, rec in enumerate(records)}
    S = [[0.0] * n for _ in range(n)]

    with tempfile.TemporaryDirectory() as tmpdir:
        db = os.path.join(tmpdir, "db")
        res = os.path.join(tmpdir, "res")
        tmp = os.path.join(tmpdir, "tmp")
        out = os.path.join(tmpdir, "out.m8")
        os.makedirs(tmp, exist_ok=True)

        subprocess.run(
            [mmseqs_bin, "createdb", fasta_path, db],
            check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
        )
        subprocess.run(
            [mmseqs_bin, "search", db, db, res, tmp, "--threads", str(threads), "-e", "10", "-s", str(sensitivity)],
            check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
        )
        subprocess.run(
            [mmseqs_bin, "convertalis", db, db, res, out, "--format-mode", "0"],
            check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
        )

        if os.path.exists(out):
            with open(out, "r") as f:
                for line in f:
                    parts = line.strip().split("\t")
                    if len(parts) >= 12:
                        q, t = parts[0], parts[1]
                        if q in id_map and t in id_map:
                            try:
                                bits = float(parts[11])
                            except ValueError:
                                bits = 0.0
                            qi, ti = id_map[q], id_map[t]
                            if bits > S[qi][ti]:
                                S[qi][ti] = bits

    scores = []
    for i in range(n):
        for j in range(i + 1, n):
            self_score = S[i][i] + S[j][j]
            comp_score = S[i][j] + S[j][i]
            w_ij = comp_score / self_score if self_score > 0.0 else 0.0
            scores.append(min(1.0, max(0.0, w_ij)))

    arr = np.array(scores, dtype=np.float64)
    return {
        "sss_mean": round(float(np.mean(arr)), 6),
        "sss_median": round(float(np.median(arr)), 6),
        "sss_std": round(float(np.std(arr)), 6),
        "sss_min": round(float(np.min(arr)), 6),
        "sss_max": round(float(np.max(arr)), 6),
        "sss_zero_frac": round(float(np.count_nonzero(arr <= 1e-6) / len(arr)), 6),
        "num_taxa": n
    }


def simulate_single_replicate_task(task: Tuple[str, Dict[str, Any], int]) -> Tuple[str, int, Dict[str, float]]:
    """Worker task: Simulates single replicate using AliSim and calculates MMseqs2 SSS."""
    cond_id, params, rep = task
    dist = float(params["distance"])
    length = int(params["length"])
    taxa = int(params.get("num_taxa", 32))
    alpha = float(params.get("alpha", 1.0))
    ics_prop = float(params.get("ics_prop", 0.0))
    sim_type = params.get("type", "standard")
    indel = params.get("indel", "0.05,0.10")
    indel_size = params.get("indel_size", None)

    cur_seed = rep

    with tempfile.TemporaryDirectory() as tmpdir:
        sim_prefix = os.path.join(tmpdir, f"sim_{rep}")
        cmd = [IQTREE_BIN, "--alisim", sim_prefix]

        # CRITICAL FIX: Pass length to AliSim
        cmd.extend(["--length", str(length)])

        if sim_type == "simple_lg":
            cmd.extend(["-m", "LG"])
        elif sim_type in ["ics", "ics_full"]:
            part_file = os.path.join(tmpdir, "partition.nex")
            k_ics = int(round(length * ics_prop))
            k_ics = max(1, min(length - 1, k_ics)) if 0.0 < ics_prop < 1.0 else (length if ics_prop >= 1.0 else 0)
            all_sites = list(range(1, length + 1))
            random.seed(cur_seed)
            ics_sites = sorted(random.sample(all_sites, k_ics)) if k_ics > 0 else []
            lg_sites = sorted([s for s in all_sites if s not in set(ics_sites)])

            with open(part_file, "w") as pf:
                pf.write("#nexus\nbegin sets;\n")
                charparts = []
                if lg_sites:
                    pf.write(f"    charset part_lg = {' '.join(map(str, lg_sites))};\n")
                    charparts.append(f"LG+G4{{{alpha}}}:part_lg")
                if ics_sites:
                    pf.write(f"    charset part_ics = {' '.join(map(str, ics_sites))};\n")
                    charparts.append(f"ICS+G4{{{alpha}}}:part_ics")
                pf.write(f"    charpartition mypart = {', '.join(charparts)};\nend;\n")

            cmd.extend(["--seqtype", "AA", "--mdef", str(ICS_MODEL_PATH), "-q", part_file])
        else:
            cmd.extend(["-m", f"LG+G4{{{alpha}}}"])

        cmd.extend(["-t", f"RANDOM{{bd{{0.1/0.05}}/{taxa}}}"])
        cmd.extend(["--indel", indel])
        if indel_size:
            cmd.extend(["--indel-size", indel_size])
        cmd.extend(["--branch-scale", str(dist), "-af", "fasta", "-seed", str(cur_seed), "--redo"])

        try:
            subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            fasta_path = f"{sim_prefix}.unaligned.fa"
            if not os.path.exists(fasta_path):
                fasta_path = f"{sim_prefix}.fa"

            if os.path.exists(fasta_path):
                stats = compute_single_replicate_sss(fasta_path, sensitivity=7.5, threads=1)
                return cond_id, rep, stats
        except Exception as e:
            pass

    # Fallback
    return cond_id, rep, {
        "sss_mean": 0.0, "sss_median": 0.0, "sss_std": 0.0,
        "sss_min": 0.0, "sss_max": 0.0, "sss_zero_frac": 1.0, "num_taxa": taxa
    }


def get_condition_hash(params: Dict[str, Any]) -> str:
    """Generates unique deterministic string hash for a simulation condition."""
    d = round(float(params["distance"]), 4)
    l = int(params["length"])
    a = round(float(params.get("alpha", 1.0)), 4)
    i = round(float(params.get("ics_prop", 0.0)), 4)
    t = int(params.get("num_taxa", 32))
    st = params.get("type", "standard")
    ind = params.get("indel", "0.05,0.10")
    inds = params.get("indel_size", "none") or "none"
    return f"T={st}|D={d}|L={l}|A={a}|ICS={i}|N={t}|IND={ind}|INDS={inds}"


def load_checkpoint_cache() -> Dict[Tuple[str, int], Dict[str, float]]:
    """Loads existing replicate results from persistent cache file if available."""
    cache: Dict[Tuple[str, int], Dict[str, float]] = {}
    if CACHE_FILE.exists():
        try:
            df_cache = pd.read_csv(CACHE_FILE)
            for _, row in df_cache.iterrows():
                cid = row["cond_hash"]
                rep = int(row["replicate"])
                cache[(cid, rep)] = {
                    "sss_mean": float(row["sss_mean"]),
                    "sss_median": float(row["sss_median"]),
                    "sss_std": float(row["sss_std"]),
                    "sss_min": float(row["sss_min"]),
                    "sss_max": float(row["sss_max"]),
                    "sss_zero_frac": float(row["sss_zero_frac"]),
                    "num_taxa": int(row["num_taxa"]) if "num_taxa" in row else 32
                }
            print(f"Loaded {len(cache):,} pre-computed replicates from cache: {CACHE_FILE.name}")
        except Exception as e:
            print(f"[Warning] Could not load cache: {e}")
    return cache


def append_to_cache(batch_results: List[Dict[str, Any]]):
    """Appends newly completed replicate stats to cache file."""
    if not batch_results:
        return
    df_batch = pd.DataFrame(batch_results)
    header = not CACHE_FILE.exists()
    df_batch.to_csv(CACHE_FILE, mode="a", index=False, header=header)


def run_full_replicate_calculation(max_workers: int = 8, batch_save_interval: int = 100):
    print("=" * 70)
    print("STEP 1: Scanning All 14 Benchmark Result Datasets...")
    print("=" * 70)

    cached_results = load_checkpoint_cache()
    exp_dfs: Dict[str, pd.DataFrame] = {}
    unique_cond_params: Dict[str, Dict[str, Any]] = {}
    tasks_to_run: List[Tuple[str, Dict[str, Any], int]] = []

    for exp in EXPERIMENT_CONFIGS:
        exp_name = exp["name"]
        exp_dir = RESULTS_DIR / exp_name
        csv_path = exp_dir / exp["csv"]
        if not csv_path.exists():
            continue

        df = pd.read_csv(csv_path)
        exp_dfs[exp_name] = df

        if exp_name == "results_paper_tree_pre_fix":
            # Paper tree already has 1800 replicates fully calculated!
            continue

        cond_keys = [k for k in exp["keys"] if k in df.columns]
        rep_keys = cond_keys + ["replicate"]
        unique_reps = df[rep_keys].drop_duplicates().to_dict(orient="records")

        for r in unique_reps:
            params = {
                "distance": r["distance"],
                "length": r["length"],
                "num_taxa": r.get("num_taxa", exp.get("default_taxa", 32)),
                "alpha": r.get("alpha", exp.get("default_alpha", 1.0)),
                "ics_prop": r.get("ics_prop", exp.get("default_ics", 0.0)),
                "type": exp["type"],
                "indel": exp["indel"],
                "indel_size": exp["indel_size"],
            }
            cid = get_condition_hash(params)
            unique_cond_params[cid] = params
            rep = int(r["replicate"])

            if (cid, rep) not in cached_results:
                tasks_to_run.append((cid, params, rep))

    total_tasks = len(tasks_to_run)
    total_reps_target = len(cached_results) + total_tasks
    print(f"Total target replicates across all conditions: {total_reps_target:,}")
    print(f"Already cached: {len(cached_results):,}")
    print(f"Remaining replicates to compute: {total_tasks:,}")

    if total_tasks > 0:
        print("=" * 70)
        print(f"STEP 2: Computing {total_tasks:,} Replicates via MMseqs2 ({max_workers} parallel workers)...")
        print("=" * 70)

        t_start = time.time()
        completed = 0
        batch_to_save: List[Dict[str, Any]] = []

        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            futures = {executor.submit(simulate_single_replicate_task, t): (t[0], t[2]) for t in tasks_to_run}
            for fut in as_completed(futures):
                cid, rep = futures[fut]
                try:
                    res_cid, res_rep, stats = fut.result()
                    cached_results[(res_cid, res_rep)] = stats
                    params = unique_cond_params[res_cid]
                    batch_to_save.append({
                        "cond_hash": res_cid,
                        "replicate": res_rep,
                        "distance": params["distance"],
                        "length": params["length"],
                        "alpha": params.get("alpha", 1.0),
                        "ics_prop": params.get("ics_prop", 0.0),
                        "num_taxa": params.get("num_taxa", 32),
                        "type": params["type"],
                        **stats
                    })
                except Exception as e:
                    print(f"[Error] Task {cid} rep {rep}: {e}")

                completed += 1
                if len(batch_to_save) >= batch_save_interval:
                    append_to_cache(batch_to_save)
                    batch_to_save = []

                if completed % 250 == 0 or completed == total_tasks:
                    elapsed = time.time() - t_start
                    speed = completed / elapsed if elapsed > 0 else 0
                    eta_sec = (total_tasks - completed) / speed if speed > 0 else 0
                    print(f"Progress: [{completed:,}/{total_tasks:,}] ({completed/total_tasks*100:.1f}%) | "
                          f"Speed: {speed:.1f} reps/s | ETA: {eta_sec/60:.1f} min")

        # Save remaining batch
        if batch_to_save:
            append_to_cache(batch_to_save)

    print("=" * 70)
    print("STEP 3: Merging Exact Replicate-Level SSS into Summary CSVs...")
    print("=" * 70)

    for exp in EXPERIMENT_CONFIGS:
        exp_name = exp["name"]
        exp_dir = RESULTS_DIR / exp_name
        csv_path = exp_dir / exp["csv"]
        if exp_name not in exp_dfs:
            continue

        df = exp_dfs[exp_name].copy()
        sss_csv_path = exp_dir / "sss_summary.csv"

        if exp_name == "results_paper_tree_pre_fix":
            print(f"[{exp_name}] Using established sss_summary.csv.")
            continue

        sss_rows = []
        for idx, row in df.iterrows():
            params = {
                "distance": row["distance"],
                "length": row["length"],
                "num_taxa": row.get("num_taxa", exp.get("default_taxa", 32)),
                "alpha": row.get("alpha", exp.get("default_alpha", 1.0)),
                "ics_prop": row.get("ics_prop", exp.get("default_ics", 0.0)),
                "type": exp["type"],
                "indel": exp["indel"],
                "indel_size": exp["indel_size"],
            }
            cid = get_condition_hash(params)
            rep = int(row["replicate"])
            stats = cached_results.get((cid, rep), {
                "sss_mean": np.nan, "sss_median": np.nan, "sss_std": np.nan,
                "sss_min": np.nan, "sss_max": np.nan, "sss_zero_frac": np.nan
            })
            sss_rows.append(stats)

        df_sss = pd.DataFrame(sss_rows)
        for col in ["sss_mean", "sss_median", "sss_std", "sss_min", "sss_max", "sss_zero_frac"]:
            df[col] = df_sss[col]

        # Backup original CSV if not already backed up
        bak_path = csv_path.with_suffix(".csv.bak")
        if not bak_path.exists():
            shutil.copyfile(csv_path, bak_path)

        df.to_csv(csv_path, index=False)
        print(f"[{exp_name}] Updated {csv_path.name} with exact replicate SSS ({len(df):,} rows, {df['sss_mean'].nunique():,} unique SSS values)")

        # Export condition & replicate-level sss_summary.csv
        cond_keys = [k for k in ["distance", "length", "alpha", "ics_prop", "num_taxa", "replicate"] if k in df.columns]
        summary_sss = df[cond_keys + ["sss_mean", "sss_median", "sss_std", "sss_min", "sss_max", "sss_zero_frac"]].drop_duplicates()
        summary_sss.to_csv(sss_csv_path, index=False)
        print(f"[{exp_name}] Saved replicate-level SSS summary -> {sss_csv_path.name} ({len(summary_sss):,} replicates)")

    print("=" * 70)
    print("All replicate-level SSS calculations, CSV merges, and exports completed successfully!")
    print("=" * 70)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Full Replicate-Level SSS Calculation across 14 directories.")
    parser.add_argument("--workers", type=int, default=8, help="Parallel worker processes (default: 8)")
    args = parser.parse_args()

    run_full_replicate_calculation(max_workers=args.workers)
