#!/usr/bin/env python3
"""
Sequence Similarity Benchmark Script
Simulates protein sequences across:
- Exp1 Conditions: D in [0.1, 0.5, 1.0, 2.0, 3.0] x L in [100, 300, 500, 1000, 1500] (25 conditions)
- High Distance Conditions: D in [4.0, 5.0, 6.0] x L in [300, 500, 1000, 1500] (12 conditions)
(100 replicates each = 3,700 total tasks)

Computes pairwise sequence identity under:
1. True MSA (Ground Truth Alignment - Full 496 pairs)
2. MAFFT MSA (Multiple Sequence Alignment - Full 496 pairs)
3. Pairwise Alignment (PSA - Sampled pairs)

Outputs results to analysis/similarity/similarity_summary.csv
"""

import os
import sys
import argparse
import subprocess
import shutil
import tempfile
import itertools
import random
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Align import PairwiseAligner, substitution_matrices

PROJECT_ROOT = Path(__file__).resolve().parent.parent
ANALYSIS_DIR = PROJECT_ROOT / "analysis"
SIMILARITY_DIR = ANALYSIS_DIR / "similarity"
SIMILARITY_DIR.mkdir(parents=True, exist_ok=True)

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

_BLOSUM62 = substitution_matrices.load("BLOSUM62")

def get_nw_aligner(gap_open=10.0, gap_extend=0.5):
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    aligner.substitution_matrix = _BLOSUM62
    aligner.open_gap_score = -abs(gap_open)
    aligner.extend_gap_score = -abs(gap_extend)
    return aligner

def compute_pairwise_identities_from_alignment(aligned_seqs):
    """
    Computes pairwise identities (% identical sites over aligned non-gap columns) for all 496 pairs.
    """
    N = len(aligned_seqs)
    identities = []
    gap_proportions = []
    
    seq_arrs = [np.frombuffer(s.encode('ascii'), dtype='S1') for s in aligned_seqs]
    L = len(seq_arrs[0])
    
    for i in range(N):
        gap_count = np.count_nonzero(seq_arrs[i] == b'-')
        gap_proportions.append(gap_count / L)
        
        for j in range(i + 1, N):
            s1 = seq_arrs[i]
            s2 = seq_arrs[j]
            valid_mask = (s1 != b'-') & (s2 != b'-')
            valid_len = np.count_nonzero(valid_mask)
            if valid_len > 0:
                match_count = np.count_nonzero((s1 == s2) & valid_mask)
                identities.append(match_count / valid_len)
            else:
                identities.append(0.0)
                
    return identities, gap_proportions

def compute_psa_identities(unaligned_seqs, aligner, max_pairs=15):
    """
    Computes pairwise identities from Needleman-Wunsch pairwise alignments for sampled pairs.
    """
    N = len(unaligned_seqs)
    all_pairs = list(itertools.combinations(range(N), 2))
    if len(all_pairs) > max_pairs:
        pairs = random.sample(all_pairs, max_pairs)
    else:
        pairs = all_pairs
        
    psa_identities = []
    for i, j in pairs:
        s1 = unaligned_seqs[i]
        s2 = unaligned_seqs[j]
        alns = aligner.align(s1, s2)
        top = alns[0]
        
        a1 = np.frombuffer(top[0].encode('ascii'), dtype='S1')
        a2 = np.frombuffer(top[1].encode('ascii'), dtype='S1')
        valid = (a1 != b'-') & (a2 != b'-')
        valid_len = np.count_nonzero(valid)
        if valid_len > 0:
            matches = np.count_nonzero((a1 == a2) & valid)
            psa_identities.append(matches / valid_len)
        else:
            psa_identities.append(0.0)
            
    return psa_identities

def process_single_similarity_task(task_params):
    dist = task_params["distance"]
    length = task_params["length"]
    rep = task_params["replicate"]
    num_taxa = task_params.get("num_taxa", 32)
    alpha = task_params.get("alpha", 1.0)
    model = task_params.get("model", "LG+G4")
    
    sim_script = PROJECT_ROOT / "bin" / "simulate_data.py"
    mafft_bin = shutil.which("mafft") or "/opt/homebrew/Cellar/micromamba/2.8.1/envs/phylomethod_env/bin/mafft"
    
    try:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            out_tree = tmp_path / "true_tree.nwk"
            out_fasta = tmp_path / "unaligned.fasta"
            out_true_msa = tmp_path / "true_msa.fasta"
            out_true_mat = tmp_path / "true_matrix.phy"
            out_mafft_msa = tmp_path / "mafft_msa.fasta"
            
            # 1. AliSim Simulation
            cmd_sim = [
                sys.executable, str(sim_script),
                "--num_taxa", str(num_taxa),
                "--distance", str(dist),
                "--length", str(length),
                "--model", model,
                "--alpha", str(alpha),
                "--outtree", str(out_tree),
                "--outfasta", str(out_fasta),
                "--outtrue_msa", str(out_true_msa),
                "--outtrue_matrix", str(out_true_mat),
                "--seed", str(rep)
            ]
            res = subprocess.run(cmd_sim, capture_output=True, env=os.environ, text=True)
            if res.returncode != 0:
                print(f"⚠️ Simulation retry fallback for D={dist}, L={length}, rep={rep}: {res.stderr[:100]}", flush=True)
                return None
            
            # 2. True MSA Identities (Full 496 pairs)
            true_records = list(SeqIO.parse(out_true_msa, "fasta"))
            true_seqs = [str(r.seq) for r in true_records]
            true_idents, true_gaps = compute_pairwise_identities_from_alignment(true_seqs)
            
            # 3. MAFFT Alignment & Identities (Full 496 pairs)
            with open(out_mafft_msa, "w") as f:
                subprocess.run([mafft_bin, "--auto", "--quiet", str(out_fasta)], stdout=f, check=True, env=os.environ)
                
            mafft_records = list(SeqIO.parse(out_mafft_msa, "fasta"))
            mafft_seqs = [str(r.seq) for r in mafft_records]
            mafft_idents, mafft_gaps = compute_pairwise_identities_from_alignment(mafft_seqs)
            
            # 4. Pairwise Sequence Alignment (PSA) Identities (Needleman-Wunsch with BLOSUM62)
            unaligned_seqs = [str(r.seq) for r in SeqIO.parse(out_fasta, "fasta")]
            aligner = get_nw_aligner(gap_open=10.0, gap_extend=0.5)
            psa_idents = compute_psa_identities(unaligned_seqs, aligner, max_pairs=15)

            return {
                "distance": dist,
                "length": length,
                "replicate": rep,
                "num_taxa": num_taxa,
                "alpha": alpha,
                "model": model,
                # True MSA (Ground Truth)
                "true_identity_mean": float(np.mean(true_idents)),
                "true_identity_min": float(np.min(true_idents)),
                "true_identity_max": float(np.max(true_idents)),
                "true_identity_std": float(np.std(true_idents)),
                "true_gap_mean": float(np.mean(true_gaps)),
                "true_aligned_length": len(true_seqs[0]),
                # MAFFT MSA
                "mafft_identity_mean": float(np.mean(mafft_idents)),
                "mafft_identity_min": float(np.min(mafft_idents)),
                "mafft_identity_max": float(np.max(mafft_idents)),
                "mafft_identity_std": float(np.std(mafft_idents)),
                "mafft_gap_mean": float(np.mean(mafft_gaps)),
                "mafft_aligned_length": len(mafft_seqs[0]),
                # PSA
                "psa_identity_mean": float(np.mean(psa_idents)),
                "psa_identity_min": float(np.min(psa_idents)),
                "psa_identity_max": float(np.max(psa_idents)),
                "psa_identity_std": float(np.std(psa_idents))
            }
    except Exception as e:
        print(f"Error processing D={dist}, L={length}, rep={rep}: {e}", flush=True)
        return None

def main():
    parser = argparse.ArgumentParser(description="Compute Sequence Similarity / Identity across D and L")
    parser.add_argument("--replicates", type=int, default=100, help="Number of replicates per condition (default: 100)")
    parser.add_argument("--threads", type=int, default=8, help="Number of parallel worker threads (default: 8)")
    parser.add_argument("--outcsv", default=str(SIMILARITY_DIR / "similarity_summary.csv"), help="Output CSV path")
    args = parser.parse_args()

    # Exp1 conditions (D <= 3.0 has L in [100, 300, 500, 1000, 1500])
    # Exp8 conditions (D >= 4.0 has L in [300, 500, 1000, 1500])
    conditions = []
    for d in [0.1, 0.5, 1.0, 2.0, 3.0]:
        for l in [100, 300, 500, 1000, 1500]:
            conditions.append((d, l))
            
    for d in [4.0, 5.0, 6.0]:
        for l in [300, 500, 1000, 1500]:
            conditions.append((d, l))

    tasks = []
    for d, l in conditions:
        for rep in range(1, args.replicates + 1):
            tasks.append({
                "distance": d,
                "length": l,
                "replicate": rep,
                "num_taxa": 32,
                "alpha": 1.0,
                "model": "LG+G4"
            })

    print(f"🚀 Starting sequence similarity benchmark: {len(tasks)} total tasks ({len(conditions)} conditions x {args.replicates} replicates) with {args.threads} threads...", flush=True)

    results = []
    completed = 0
    with ThreadPoolExecutor(max_workers=args.threads) as executor:
        futures = [executor.submit(process_single_similarity_task, t) for t in tasks]
        for f in as_completed(futures):
            res = f.result()
            if res is not None:
                results.append(res)
            completed += 1
            if completed % 200 == 0 or completed == len(tasks):
                print(f"  [{completed:>4}/{len(tasks)}] tasks completed ({completed/len(tasks)*100:.1f}%)", flush=True)

    df = pd.DataFrame(results)
    df = df.sort_values(by=["distance", "length", "replicate"]).reset_index(drop=True)
    df.to_csv(args.outcsv, index=False)
    print(f"\n✅ All sequence similarity computations completed! Successfully wrote {len(df)} records to -> {args.outcsv}", flush=True)

if __name__ == "__main__":
    main()
