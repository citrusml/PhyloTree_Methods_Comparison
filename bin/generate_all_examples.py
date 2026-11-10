#!/usr/bin/env python3
"""
Script to generate all example simulation conditions organized hierarchically by experiment.
"""
import os
import sys
import json
import shutil
import itertools
import subprocess
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed

import dendropy
from dendropy.calculate import treecompare
from Bio import SeqIO

PROJECT_ROOT = Path(__file__).resolve().parent.parent
ANALYSIS_DIR = PROJECT_ROOT / "analysis"
EXAMPLE_DIR = ANALYSIS_DIR / "example"

env_bin_paths = [
    "/opt/homebrew/Cellar/micromamba/2.8.1/envs/phylomethod_env/bin",
    os.path.expanduser("~/miniconda3/envs/phylomethod_env/bin"),
    os.path.expanduser("~/miniforge3/envs/phylomethod_env/bin"),
    os.path.expanduser("~/micromamba/envs/phylomethod_env/bin"),
    "/opt/homebrew/bin",
    os.path.expanduser("~/bin")
]
for p in env_bin_paths:
    if os.path.isdir(p) and p not in os.environ["PATH"]:
        os.environ["PATH"] = p + ":" + os.environ["PATH"]

def compute_unrooted_nrf(true_tree_file, est_tree_file):
    """DendroPy を用いて真の系統樹と推定系統樹間の Robinson-Foulds (RF/nRF) 距離を計算"""
    tns = dendropy.TaxonNamespace()
    t_true = dendropy.Tree.get(path=str(true_tree_file), schema="newick", taxon_namespace=tns, preserve_underscores=True)
    t_est = dendropy.Tree.get(path=str(est_tree_file), schema="newick", taxon_namespace=tns, preserve_underscores=True)
    t_true.is_rooted = False
    t_true.deroot()
    t_est.is_rooted = False
    t_est.deroot()
    t_true.encode_bipartitions()
    t_est.encode_bipartitions()
    rf = treecompare.symmetric_difference(t_true, t_est)
    num_taxa = len(tns)
    max_rf = 2 * (num_taxa - 3)
    nrf = rf / max_rf if max_rf > 0 else 0.0
    return rf, nrf

def generate_single_example(task_tuple):
    """
    単一条件のシミュレーションおよび推定パイプラインを実行
    """
    exp_id, cond_name, params, example_base_dir_str, project_root_str = task_tuple
    
    cond_dir = Path(example_base_dir_str) / exp_id / cond_name
    cond_dir.mkdir(parents=True, exist_ok=True)
    
    project_root = Path(project_root_str)
    sim_script = project_root / "bin" / "simulate_data.py"
    pwa_script = project_root / "bin" / "run_pwa_nj.py"
    msa_nj_script = project_root / "bin" / "run_msa_nj.py"
    
    out_tree = cond_dir / "true_tree.nwk"
    out_fasta = cond_dir / "unaligned.fasta"
    out_true_msa = cond_dir / "true_msa.fasta"
    out_true_mat = cond_dir / "true_matrix.phy"
    mafft_msa = cond_dir / "mafft_msa.fasta"
    pwa_tree = cond_dir / "pwa_nj_tree.nwk"
    msa_nj_tree = cond_dir / "msa_nj_tree.nwk"
    fastme_tree = cond_dir / "pwa_fastme_tree.nwk"
    meta_json = cond_dir / "summary.json"
    
    # 既存の有効な結果があれば再利用
    if meta_json.exists() and pwa_tree.exists() and msa_nj_tree.exists():
        try:
            with open(meta_json, "r") as f:
                return json.load(f)
        except Exception:
            pass

    # 1. AliSim によるシミュレーション
    cmd_sim = [
        sys.executable, str(sim_script),
        "--distance", str(params.get("distance", 1.0)),
        "--length", str(params.get("length", 500)),
        "--num_taxa", str(params.get("num_taxa", 32)),
        "--alpha", str(params.get("alpha", 1.0)),
        "--model", str(params.get("model", "LG+G4")),
        "--ics_prop", str(params.get("ics_prop", 0.0)),
        "--outtree", str(out_tree),
        "--outfasta", str(out_fasta),
        "--outtrue_msa", str(out_true_msa),
        "--outtrue_matrix", str(out_true_mat),
        "--seed", str(params.get("seed", 42))
    ]
    subprocess.run(cmd_sim, check=True, capture_output=True, env=os.environ)
    
    # 2. MAFFT による MSA 構築
    mafft_bin = shutil.which("mafft") or "/opt/homebrew/Cellar/micromamba/2.8.1/envs/phylomethod_env/bin/mafft"
    with open(mafft_msa, "w") as f_out:
        subprocess.run([mafft_bin, "--auto", "--quiet", str(out_fasta)], stdout=f_out, check=True, env=os.environ)
        
    # 3. PWA+NJ の実行 (RapidNJ)
    cmd_pwa = [
        sys.executable, str(pwa_script),
        "--fasta", str(out_fasta),
        "--dist_model", params.get("dist_model", "poisson"),
        "--alpha", str(params.get("alpha", 1.0)),
        "--tool", "rapidnj",
        "--outtree", str(pwa_tree)
    ]
    subprocess.run(cmd_pwa, check=True, capture_output=True, env=os.environ)

    # 4. MSA+NJ の実行 (RapidNJ)
    cmd_msa_nj = [
        sys.executable, str(msa_nj_script),
        "--msa", str(mafft_msa),
        "--tool", "rapidnj",
        "--outtree", str(msa_nj_tree)
    ]
    subprocess.run(cmd_msa_nj, check=True, capture_output=True, env=os.environ)

    # 5. PWA+FastME (FastME 利用可能時)
    fastme_bin = shutil.which("fastme") or "/opt/homebrew/Cellar/micromamba/2.8.1/envs/phylomethod_env/bin/fastme"
    if os.path.exists(fastme_bin):
        cmd_pwa_fastme = [
            sys.executable, str(pwa_script),
            "--fasta", str(out_fasta),
            "--dist_model", params.get("dist_model", "poisson"),
            "--alpha", str(params.get("alpha", 1.0)),
            "--tool", "fastme",
            "--outtree", str(fastme_tree)
        ]
        subprocess.run(cmd_pwa_fastme, check=True, capture_output=True, env=os.environ)

    # 6. アライメント長 & 精度評価
    records = list(SeqIO.parse(mafft_msa, "fasta"))
    msa_len = len(records[0].seq)
    
    _, pwa_nrf = compute_unrooted_nrf(out_tree, pwa_tree)
    _, msa_nj_nrf = compute_unrooted_nrf(out_tree, msa_nj_tree)

    record = {
        "experiment_id": exp_id,
        "condition_name": cond_name,
        "num_taxa": params.get("num_taxa", 32),
        "distance": params.get("distance", 1.0),
        "length": params.get("length", 500),
        "alpha": params.get("alpha", 1.0),
        "ics_prop": params.get("ics_prop", 0.0),
        "model": params.get("model", "LG+G4"),
        "dist_model": params.get("dist_model", "poisson"),
        "msa_aligned_length": msa_len,
        "pwa_nrf": round(pwa_nrf, 4),
        "msa_nj_nrf": round(msa_nj_nrf, 4),
        "dir": str(cond_dir)
    }
    
    with open(meta_json, "w") as f:
        json.dump(record, f, indent=2)
        
    return record

def build_all_tasks(example_dir=EXAMPLE_DIR, project_root=PROJECT_ROOT):
    all_experiment_definitions = {
        'exp1_default': {
            'distances': [0.1, 0.5, 1.0, 2.0, 3.0],
            'lengths': [100, 300, 500, 1000, 1500],
            'num_taxa': [32],
            'alphas': [1.0],
            'ics_props': [0.0],
            'model': 'LG+G4',
            'dist_model': 'poisson',
        },
        'exp2_taxon': {
            'distances': [0.1, 0.5, 1.0, 2.0, 3.0],
            'lengths': [100, 500, 1000],
            'num_taxa': [8, 16, 64, 128],
            'alphas': [1.0],
            'ics_props': [0.0],
            'model': 'LG+G4',
            'dist_model': 'poisson',
        },
        'exp3_alpha': {
            'distances': [0.1, 0.5, 1.0, 2.0, 3.0],
            'lengths': [100, 500, 1000],
            'num_taxa': [32],
            'alphas': [0.25, 0.5, 1.0, 2.0],
            'ics_props': [0.0],
            'model': 'LG+G4',
            'dist_model': 'poisson',
        },
        'exp5_gamma': {
            'distances': [0.1, 0.5, 1.0, 2.0, 3.0],
            'lengths': [100, 500, 1000],
            'num_taxa': [32],
            'alphas': [0.25, 0.5, 1.0, 2.0],
            'ics_props': [0.0],
            'model': 'LG+G4',
            'dist_model': 'gamma_poisson',
        },
        'exp6_ics': {
            'distances': [0.1, 0.5, 1.0, 2.0, 3.0],
            'lengths': [100, 300, 500, 1000],
            'num_taxa': [32],
            'alphas': [1.0],
            'ics_props': [0.0, 0.05, 0.1, 0.2],
            'model': 'LG+G4',
            'dist_model': 'poisson',
        },
        'exp7_fastme': {
            'distances': [0.1, 0.5, 1.0, 2.0, 3.0],
            'lengths': [100, 300, 500, 1000, 1500],
            'num_taxa': [32],
            'alphas': [1.0],
            'ics_props': [0.0],
            'model': 'LG+G4',
            'dist_model': 'poisson',
        },
        'exp8_high_dist': {
            'distances': [4.0, 5.0, 6.0],
            'lengths': [300, 500, 1000, 1500],
            'num_taxa': [32],
            'alphas': [1.0],
            'ics_props': [0.0],
            'model': 'LG+G4',
            'dist_model': 'poisson',
        },
        'exp9_simple_lg': {
            'distances': [0.1, 0.5, 1.0, 2.0, 3.0],
            'lengths': [100, 300, 500, 1000, 1500],
            'num_taxa': [32],
            'alphas': [1.0],
            'ics_props': [0.0],
            'model': 'LG',
            'dist_model': 'poisson',
        }
    }

    simulation_tasks = []

    for exp_id, cfg in all_experiment_definitions.items():
        combos = list(itertools.product(cfg['distances'], cfg['lengths'], cfg['num_taxa'], cfg['alphas'], cfg['ics_props']))
        for d, l, n, a, ics in combos:
            parts = []
            if len(cfg['num_taxa']) > 1:
                parts.append(f"N{n}")
            if len(cfg['alphas']) > 1:
                parts.append(f"a{a}")
            if len(cfg['ics_props']) > 1:
                parts.append(f"ics{ics}")
            parts.append(f"D{d}_L{l}")
            cond_name = "_".join(parts)
            
            p = {
                "distance": d,
                "length": l,
                "num_taxa": n,
                "alpha": a,
                "ics_prop": ics,
                "model": cfg["model"],
                "dist_model": cfg["dist_model"],
                "seed": 42
            }
            simulation_tasks.append((exp_id, cond_name, p, str(example_dir), str(project_root)))
            
    return simulation_tasks

def run_all(max_workers=8):
    tasks = build_all_tasks()
    print(f"🚀 Generating all {len(tasks)} example conditions across 8 experiments with {max_workers} threads...")
    
    results = []
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(generate_single_example, t) for t in tasks]
        done = 0
        for f in as_completed(futures):
            res = f.result()
            results.append(res)
            done += 1
            if done % 50 == 0 or done == len(tasks):
                print(f"  [{done:>3}/{len(tasks)}] conditions completed.")
                
    return results

if __name__ == "__main__":
    run_all()
