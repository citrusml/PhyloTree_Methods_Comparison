#!/usr/bin/env python3
"""
Sequence Similarity Score (SSS) Calculation Script: calculate_sss.py

Calculates the symmetrical relative Sequence Similarity Score (SSS) w_ij and its
average across all sequence pairs for simulated replicates, exactly matching the
benchmark definition of Matsui & Iwasaki (2020, Systematic Biology, syz049; Dufour et al. 2010)
and the Graph Splitting tool (gs2):

    w_ij = (S_bits(i, j) + S_bits(j, i)) / (S_bits(i, i) + S_bits(j, j))

    average SSS (w_bar) = 2 / (N * (N - 1)) * sum_{i < j} w_ij

Pairwise sequence alignment bit scores S_bits are calculated by all-to-all PSA using MMseqs2,
strictly adhering to the original paper and gs2 implementation (src/format.cpp bl2mat).

Outputs summary statistics (mean, median, min, max, fraction of zero-similarity pairs)
per replicate to CSV for benchmark aggregation and correlation with tree reconstruction error.
"""

from typing import List, Dict, Any, Tuple, Optional
import os
import sys
import glob
import shutil
import tempfile
import subprocess
import argparse
import numpy as np
import pandas as pd
from Bio import SeqIO


def find_mmseqs_binary() -> Optional[str]:
    """
    Finds mmseqs binary in PATH or conda/micromamba environments.

    Returns
    -------
    Optional[str]
        Path to mmseqs executable if found, else None.
    """
    bin_path = shutil.which("mmseqs")
    if bin_path:
        return bin_path
    py_dir = os.path.dirname(sys.executable)
    cand = os.path.join(py_dir, "mmseqs")
    if os.path.isfile(cand) and os.access(cand, os.X_OK):
        return cand
    env_dir = os.environ.get("CONDA_PREFIX")
    if env_dir:
        cand = os.path.join(env_dir, "bin", "mmseqs")
        if os.path.isfile(cand) and os.access(cand, os.X_OK):
            return cand
    user_home = os.path.expanduser("~")
    cand = os.path.join(user_home, ".micromamba", "envs", "phylomethod_env", "bin", "mmseqs")
    if os.path.isfile(cand) and os.access(cand, os.X_OK):
        return cand
    return None


def compute_sss_for_fasta(
    fasta_path: str,
    sensitivity: float = 7.5,
    threads: int = 1,
    mmseqs_bin: Optional[str] = None
) -> Dict[str, float]:
    """
    Computes pairwise Sequence Similarity Scores (SSS) w_ij using MMseqs2,
    strictly matching the benchmark implementation in Matsui & Iwasaki (2020)
    and gs2 (src/format.cpp bl2mat):

        w_ij = (S_bits(i, j) + S_bits(j, i)) / (S_bits(i, i) + S_bits(j, j))

    Parameters
    ----------
    fasta_path : str
        Path to input FASTA file.
    sensitivity : float, default=7.5
        MMseqs2 search sensitivity parameter (-s).
    threads : int, default=1
        Number of threads for MMseqs2.
    mmseqs_bin : Optional[str]
        Path to mmseqs binary. If None, searched automatically.

    Returns
    -------
    Dict[str, float]
        Dictionary with sss_mean, sss_median, sss_std, sss_min, sss_max, sss_zero_frac, num_taxa.
    """
    if mmseqs_bin is None:
        mmseqs_bin = find_mmseqs_binary()
    if not mmseqs_bin:
        raise RuntimeError(
            "ERROR!!: mmseqs binary not found. Please ensure MMseqs2 is installed and available in PATH."
        )

    records = list(SeqIO.parse(fasta_path, "fasta"))
    seqs = [str(r.seq).replace("-", "").upper() for r in records]
    if any(len(s) == 0 for s in seqs):
        raise ValueError("ERROR!!: this replication contains 0 length")

    n = len(seqs)
    if n < 2:
        return {
            "sss_mean": 0.0,
            "sss_median": 0.0,
            "sss_std": 0.0,
            "sss_min": 0.0,
            "sss_max": 0.0,
            "sss_zero_frac": 1.0,
            "num_taxa": n
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

    # Compute symmetrical relative sequence similarity score w_ij
    # exactly matching Matsui & Iwasaki (2020) and gs2 bl2mat:
    # w_ij = (S_bits(i, j) + S_bits(j, i)) / (S_bits(i, i) + S_bits(j, j))
    scores: List[float] = []
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


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Calculate SSS (Sequence Similarity Score) per Replicate strictly using MMseqs2 (Matsui & Iwasaki 2020)"
    )
    parser.add_argument("--fasta", help="Single FASTA file to evaluate")
    parser.add_argument("--fastas", nargs="*", help="List of FASTA files to evaluate")
    parser.add_argument("--rep_start", type=int, help="Start replicate ID for chunk")
    parser.add_argument("--rep_end", type=int, help="End replicate ID for chunk")
    parser.add_argument("--distance", type=float, required=True, help="Evolutionary distance D")
    parser.add_argument("--length", type=int, required=True, help="Initial sequence length L")
    parser.add_argument("--sensitivity", type=float, default=7.5, help="MMseqs2 sensitivity parameter (default: 7.5)")
    parser.add_argument("--threads", type=int, default=1, help="Threads for MMseqs2 (default: 1)")
    parser.add_argument("--outcsv", required=True, help="Output CSV path")
    args = parser.parse_args()

    records: List[Dict[str, Any]] = []

    # Case 1: Chunk execution with --rep_start and --rep_end
    if args.rep_start is not None and args.rep_end is not None:
        for rep in range(args.rep_start, args.rep_end + 1):
            fasta_candidate = f"seqs_{rep}.fasta"
            if not os.path.exists(fasta_candidate):
                alt = glob.glob(f"*_{rep}.fasta")
                if alt:
                    fasta_candidate = alt[0]
                else:
                    print(f"[Warning] FASTA not found for replicate {rep} (tried {fasta_candidate})")
                    continue

            stats = compute_sss_for_fasta(
                fasta_candidate,
                sensitivity=args.sensitivity,
                threads=args.threads
            )
            row = {
                "distance": args.distance,
                "length": args.length,
                "replicate": rep,
                **stats
            }
            records.append(row)

    # Case 2: List of fastas
    elif args.fastas:
        for fpath in args.fastas:
            fname = os.path.basename(fpath)
            rep = 1
            parts = fname.replace(".fasta", "").replace(".fa", "").split("_")
            for p in reversed(parts):
                if p.isdigit():
                    rep = int(p)
                    break
            stats = compute_sss_for_fasta(
                fpath,
                sensitivity=args.sensitivity,
                threads=args.threads
            )
            records.append({
                "distance": args.distance,
                "length": args.length,
                "replicate": rep,
                **stats
            })

    # Case 3: Single fasta
    elif args.fasta:
        stats = compute_sss_for_fasta(
            args.fasta,
            sensitivity=args.sensitivity,
            threads=args.threads
        )
        records.append({
            "distance": args.distance,
            "length": args.length,
            "replicate": 1,
            **stats
        })

    if not records:
        print("[Warning] No SSS records computed.")
        df = pd.DataFrame(columns=[
            "distance", "length", "replicate", "num_taxa",
            "sss_mean", "sss_median", "sss_std", "sss_min", "sss_max", "sss_zero_frac"
        ])
    else:
        df = pd.DataFrame(records)

    header = not os.path.exists(args.outcsv)
    os.makedirs(os.path.dirname(os.path.abspath(args.outcsv)), exist_ok=True)
    df.to_csv(args.outcsv, mode="a", index=False, header=header)
    print(f"Computed SSS (MMseqs2) for {len(records)} replicates -> {args.outcsv}")


if __name__ == "__main__":
    main()
