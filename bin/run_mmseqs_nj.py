#!/usr/bin/env python3
"""
MMseqs2 Pairwise Alignment + Neighbor-Joining (MMSEQS+NJ) Pipeline Script.

Strictly implements the paper's pairwise sequence alignment (NJ_PSA) methodology
from Matsui & Iwasaki (2020, Systematic Biology):
1. Runs all-to-all MMseqs2 local alignment with sensitivity (-s 7.5, -e 10).
2. Computes symmetrical relative sequence similarity score:
   w_ij = (S_bits(i, j) + S_bits(j, i)) / (S_bits(i, i) + S_bits(j, j))
3. Transforms similarity w_ij into an additive evolutionary distance:
   - Log-transformed (default): d_ij = -ln(max(w_ij, epsilon))
   - Linear complement:         d_ij = 1.0 - w_ij
4. Exports PHYLIP distance matrix and infers Neighbor-Joining tree using RapidNJ or FastME.
"""

from __future__ import annotations

import sys
import os
import math
import time
import json
import argparse
import tempfile
import shutil
import subprocess
import numpy as np
from typing import Dict, Any, List, Optional
from Bio import SeqIO


def find_mmseqs_binary() -> str:
    """Locates MMseqs2 executable in PATH or conda environments."""
    cand = shutil.which("mmseqs")
    if cand and os.path.isfile(cand) and os.access(cand, os.X_OK):
        return cand

    py_dir = os.path.dirname(sys.executable)
    local_cand = os.path.join(py_dir, "mmseqs")
    if os.path.isfile(local_cand) and os.access(local_cand, os.X_OK):
        return local_cand

    raise RuntimeError("MMseqs2 binary 'mmseqs' was not found in PATH. Please ensure mmseqs2 is installed.")


def find_nj_binary(tool: str = "rapidnj") -> str:
    """Locates RapidNJ or FastME executable."""
    cand = shutil.which(tool)
    if cand and os.path.isfile(cand) and os.access(cand, os.X_OK):
        return cand

    py_dir = os.path.dirname(sys.executable)
    local_cand = os.path.join(py_dir, tool)
    if os.path.isfile(local_cand) and os.access(local_cand, os.X_OK):
        return local_cand

    raise RuntimeError(f"NJ tool binary '{tool}' was not found in PATH.")


def compute_mmseqs_similarity_matrix(
    fasta_path: str,
    sensitivity: float = 7.5,
    threads: int = 1,
) -> tuple[List[str], np.ndarray, Dict[str, float]]:
    """
    Computes all-to-all MMseqs2 symmetrical relative sequence similarity score w_ij.

    Returns
    -------
    taxa_names : List[str]
    w_matrix : np.ndarray of shape (N, N)
    summary_stats : Dict[str, float]
    """
    records = list(SeqIO.parse(fasta_path, "fasta"))
    if not records:
        raise ValueError(f"No FASTA records found in {fasta_path}")
    if any(len(rec.seq) == 0 for rec in records):
        raise ValueError(f"Zero-length sequence detected in {fasta_path}")

    taxa_names = [rec.id for rec in records]
    n = len(taxa_names)
    id_map = {name: i for i, name in enumerate(taxa_names)}

    mmseqs_bin = find_mmseqs_binary()
    S = np.zeros((n, n), dtype=np.float64)

    with tempfile.TemporaryDirectory(prefix="mmseqs_nj_") as tmpdir:
        db = os.path.join(tmpdir, "db")
        res = os.path.join(tmpdir, "res")
        tmp = os.path.join(tmpdir, "tmp")
        out_m8 = os.path.join(tmpdir, "out.m8")
        os.makedirs(tmp, exist_ok=True)

        # 1. createdb
        subprocess.run(
            [mmseqs_bin, "createdb", os.path.abspath(fasta_path), db],
            check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True
        )
        # 2. search
        subprocess.run(
            [mmseqs_bin, "search", db, db, res, tmp, "--threads", str(threads), "-e", "10", "-s", str(sensitivity)],
            check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True
        )
        # 3. convertalis
        subprocess.run(
            [mmseqs_bin, "convertalis", db, db, res, out_m8, "--format-mode", "0"],
            check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True
        )

        if os.path.exists(out_m8):
            with open(out_m8, "r") as f:
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
                            if bits > S[qi, ti]:
                                S[qi, ti] = bits

    # Compute symmetrical relative sequence similarity score w_ij
    # matching Matsui & Iwasaki (2020) and gs2 bl2mat:
    # w_ij = (S_bits(i, j) + S_bits(j, i)) / (S_bits(i, i) + S_bits(j, j))
    w_mat = np.zeros((n, n), dtype=np.float64)
    off_diag_scores: List[float] = []

    for i in range(n):
        w_mat[i, i] = 1.0
        for j in range(i + 1, n):
            self_score = S[i, i] + S[j, j]
            comp_score = S[i, j] + S[j, i]
            val = (comp_score / self_score) if self_score > 0.0 else 0.0
            val = min(1.0, max(0.0, float(val)))
            w_mat[i, j] = val
            w_mat[j, i] = val
            off_diag_scores.append(val)

    scores_arr = np.array(off_diag_scores, dtype=np.float64) if off_diag_scores else np.array([0.0])
    summary_stats = {
        "sss_mean": float(np.mean(scores_arr)),
        "sss_median": float(np.median(scores_arr)),
        "sss_std": float(np.std(scores_arr)),
        "sss_min": float(np.min(scores_arr)),
        "sss_max": float(np.max(scores_arr)),
        "sss_zero_frac": float(np.mean(scores_arr == 0.0)),
        "num_taxa": float(n),
    }

    return taxa_names, w_mat, summary_stats


def similarity_to_distance_matrix(
    w_mat: np.ndarray,
    dist_model: str = "log",
    epsilon: float = 1e-4,
) -> np.ndarray:
    """
    Converts sequence similarity matrix w_ij into distance matrix d_ij.

    Parameters
    ----------
    w_mat : np.ndarray
        Symmetrical relative sequence similarity matrix in range [0.0, 1.0].
    dist_model : str
        'log': d_ij = -ln(max(w_ij, epsilon))
        'linear': d_ij = 1.0 - w_ij
    epsilon : float
        Small positive floor to prevent -ln(0) divergence on completely diverged pairs.

    Returns
    -------
    d_mat : np.ndarray
    """
    n = w_mat.shape[0]
    d_mat = np.zeros((n, n), dtype=np.float64)

    for i in range(n):
        for j in range(n):
            if i == j:
                d_mat[i, j] = 0.0
            else:
                w = w_mat[i, j]
                if dist_model == "linear":
                    d_mat[i, j] = max(0.0, 1.0 - w)
                else:  # default "log"
                    eff_w = max(w, epsilon)
                    d_mat[i, j] = -math.log(eff_w)

    return d_mat


def write_phylip_matrix(taxa_names: List[str], d_mat: np.ndarray, out_file: str) -> None:
    """Writes standard PHYLIP square distance matrix."""
    n = len(taxa_names)
    with open(out_file, "w") as f:
        f.write(f"   {n}\n")
        for i, name in enumerate(taxa_names):
            # PHYLIP requires taxon name padded to 10 chars, or followed by whitespace
            clean_name = name.strip()
            row_str = " ".join(f"{d_mat[i, j]:.6f}" for j in range(n))
            f.write(f"{clean_name:<10} {row_str}\n")


def run_nj(phylip_file: str, out_tree: str, tool: str = "rapidnj", threads: int = 1) -> None:
    """Executes RapidNJ or FastME on PHYLIP distance matrix."""
    bin_path = find_nj_binary(tool)

    if tool == "rapidnj":
        cmd = [bin_path, os.path.abspath(phylip_file), "-i", "pd", "-o", "n", "-c", str(threads)]
        with open(out_tree, "w") as f:
            res = subprocess.run(cmd, stdout=f, stderr=subprocess.PIPE, text=True)
        if res.returncode != 0:
            raise RuntimeError(f"RapidNJ execution failed (exit code {res.returncode}):\n{res.stderr}")
    elif tool == "fastme":
        # FastME output tree flag is -o
        cmd = [bin_path, "-i", os.path.abspath(phylip_file), "-o", os.path.abspath(out_tree), "-s"]
        res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if res.returncode != 0 or not os.path.exists(out_tree):
            raise RuntimeError(f"FastME execution failed (exit code {res.returncode}):\n{res.stderr}")


def create_star_tree(taxa_names: List[str], out_tree_path: str) -> None:
    """Generates fallback star tree if tree inference fails."""
    taxa_str = ",".join([f"{t}:0.001" for t in taxa_names])
    with open(out_tree_path, "w") as f:
        f.write(f"({taxa_str});\n")


def main():
    parser = argparse.ArgumentParser(
        description="MMseqs2 Pairwise Alignment + Neighbor-Joining (MMSEQS+NJ) Pipeline"
    )
    parser.add_argument("--fasta", required=True, help="Input unaligned FASTA file")
    parser.add_argument("--outtree", required=True, help="Output Newick tree file")
    parser.add_argument("--outmatrix", help="Output PHYLIP distance matrix file")
    parser.add_argument("--outjson", help="Output JSON metadata file")
    parser.add_argument(
        "--dist_model", choices=["log", "linear"], default="log",
        help="Distance mapping from similarity score: log (-ln(w)) or linear (1-w) (default: log)"
    )
    parser.add_argument(
        "--epsilon", type=float, default=1e-4,
        help="Epsilon floor for log transformation on 0-score pairs (default: 1e-4)"
    )
    parser.add_argument(
        "--tool", choices=["rapidnj", "fastme"], default="rapidnj",
        help="Neighbor-Joining implementation to use (default: rapidnj)"
    )
    parser.add_argument(
        "--sensitivity", type=float, default=7.5,
        help="MMseqs2 sensitivity parameter -s (default: 7.5)"
    )
    parser.add_argument(
        "--threads", type=int, default=1,
        help="CPU threads for MMseqs2 and RapidNJ (default: 1)"
    )
    args = parser.parse_args()

    t0 = time.time()

    # 1. Compute MMseqs2 pairwise similarity score matrix w_ij
    taxa_names, w_mat, summary_stats = compute_mmseqs_similarity_matrix(
        args.fasta,
        sensitivity=args.sensitivity,
        threads=args.threads,
    )

    # 2. Convert similarity matrix to distance matrix
    d_mat = similarity_to_distance_matrix(
        w_mat,
        dist_model=args.dist_model,
        epsilon=args.epsilon,
    )

    # 3. Export PHYLIP distance matrix
    phylip_out = args.outmatrix or (args.outtree + ".dist.phylip")
    write_phylip_matrix(taxa_names, d_mat, phylip_out)

    # 4. Infer Neighbor-Joining tree
    try:
        run_nj(phylip_out, args.outtree, tool=args.tool, threads=args.threads)
    except Exception as e:
        print(f"[WARNING] NJ reconstruction failed ({e}), creating fallback star tree.", file=sys.stderr)
        create_star_tree(taxa_names, args.outtree)

    elapsed = time.time() - t0

    # 5. Metadata JSON export
    metadata = {
        **summary_stats,
        "dist_model": args.dist_model,
        "epsilon": args.epsilon,
        "nj_tool": args.tool,
        "sensitivity": args.sensitivity,
        "runtime_sec": elapsed,
        "dist_mean": float(np.mean(d_mat[np.triu_indices(len(taxa_names), k=1)])),
        "dist_max": float(np.max(d_mat)),
        "dist_min": float(np.min(d_mat[np.triu_indices(len(taxa_names), k=1)])),
    }

    json_out = args.outjson or (args.outtree + ".json")
    with open(json_out, "w") as f:
        json.dump(metadata, f, indent=2)


if __name__ == "__main__":
    main()
