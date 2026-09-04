#!/usr/bin/env python3
"""
Sequence Similarity Score (SSS) Calculation Script: calculate_sss.py

Calculates the symmetrical relative Sequence Similarity Score (SSS) w_ij and its
average across all sequence pairs for simulated replicates, exactly matching the
benchmark definition of Matsui & Iwasaki (2020, Systematic Biology, syz049; Dufour et al. 2010):

    w_ij = max(0.0, S_bits(i, j)) / mean(S_bits(i, i), S_bits(j, j))

    average SSS (w_bar) = 2 / (N * (N - 1)) * sum_{i < j} w_ij

Outputs summary statistics (mean, median, min, max, fraction of zero-similarity pairs)
per replicate to CSV for benchmark aggregation and correlation with tree reconstruction error.
"""

from typing import List, Dict, Any, Tuple
import os
import sys
import glob
import argparse
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Align import PairwiseAligner, substitution_matrices


def get_blosum62_aligner(gap_open: float = 10.0, gap_extend: float = 0.5) -> PairwiseAligner:
    """
    Constructs global PairwiseAligner with BLOSUM62 matrix.

    Parameters
    ----------
    gap_open : float, default=10.0
        Gap opening penalty.
    gap_extend : float, default=0.5
        Gap extension penalty.

    Returns
    -------
    PairwiseAligner
        Configured aligner.
    """
    matrix = substitution_matrices.load("BLOSUM62")
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.substitution_matrix = matrix
    aligner.open_gap_score = -abs(gap_open)
    aligner.extend_gap_score = -abs(gap_extend)
    return aligner


def compute_sss_for_sequences(
    seqs: List[str],
    aligner: PairwiseAligner
) -> Dict[str, float]:
    """
    Computes pairwise Sequence Similarity Scores (SSS) w_ij and summary metrics.

    Parameters
    ----------
    seqs : List[str]
        List of unaligned amino acid sequences.
    aligner : PairwiseAligner
        Configured BLOSUM62 global aligner.

    Returns
    -------
    Dict[str, float]
        Dictionary with sss_mean, sss_median, sss_std, sss_min, sss_max, sss_zero_frac.
    """
    n = len(seqs)
    if n < 2:
        return {
            "sss_mean": 1.0,
            "sss_median": 1.0,
            "sss_std": 0.0,
            "sss_min": 1.0,
            "sss_max": 1.0,
            "sss_zero_frac": 0.0,
        }

    # Pre-calculate self-scores S(i, i)
    self_scores = [float(aligner.score(s, s)) for s in seqs]

    scores: List[float] = []
    for i in range(n):
        s1 = seqs[i]
        self_i = self_scores[i]
        for j in range(i + 1, n):
            s2 = seqs[j]
            self_j = self_scores[j]
            mean_self = (self_i + self_j) / 2.0
            if mean_self > 0.0:
                s_ij = float(aligner.score(s1, s2))
                w_ij = max(0.0, s_ij) / mean_self
                scores.append(min(1.0, w_ij))
            else:
                scores.append(0.0)

    arr = np.array(scores, dtype=np.float64)
    return {
        "sss_mean": round(float(np.mean(arr)), 6),
        "sss_median": round(float(np.median(arr)), 6),
        "sss_std": round(float(np.std(arr)), 6),
        "sss_min": round(float(np.min(arr)), 6),
        "sss_max": round(float(np.max(arr)), 6),
        "sss_zero_frac": round(float(np.count_nonzero(arr <= 1e-6) / len(arr)), 6),
    }


def compute_sss_for_fasta(
    fasta_path: str,
    aligner: PairwiseAligner
) -> Dict[str, float]:
    """
    Parses a FASTA file and computes SSS metrics.

    Parameters
    ----------
    fasta_path : str
        Path to FASTA file.
    aligner : PairwiseAligner
        Configured BLOSUM62 aligner.

    Returns
    -------
    Dict[str, float]
        Dictionary with SSS metrics and taxon count.
    """
    records = list(SeqIO.parse(fasta_path, "fasta"))
    seqs = [str(r.seq).replace("-", "").upper() for r in records]
    stats = compute_sss_for_sequences(seqs, aligner)
    stats["num_taxa"] = len(seqs)
    return stats


def main() -> None:
    parser = argparse.ArgumentParser(description="Calculate SSS (Sequence Similarity Score) per Replicate")
    parser.add_argument("--fasta", help="Single FASTA file to evaluate")
    parser.add_argument("--fastas", nargs="*", help="List of FASTA files to evaluate")
    parser.add_argument("--rep_start", type=int, help="Start replicate ID for chunk")
    parser.add_argument("--rep_end", type=int, help="End replicate ID for chunk")
    parser.add_argument("--distance", type=float, required=True, help="Evolutionary distance D")
    parser.add_argument("--length", type=int, required=True, help="Initial sequence length L")
    parser.add_argument("--gap_open", type=float, default=10.0, help="Gap open penalty (default: 10.0)")
    parser.add_argument("--gap_extend", type=float, default=0.5, help="Gap extend penalty (default: 0.5)")
    parser.add_argument("--outcsv", required=True, help="Output CSV path")
    args = parser.parse_args()

    aligner = get_blosum62_aligner(gap_open=args.gap_open, gap_extend=args.gap_extend)
    records: List[Dict[str, Any]] = []

    # Case 1: Chunk execution with --rep_start and --rep_end
    if args.rep_start is not None and args.rep_end is not None:
        for rep in range(args.rep_start, args.rep_end + 1):
            fasta_candidate = f"seqs_{rep}.fasta"
            if not os.path.exists(fasta_candidate):
                # Try unaligned fallback or general match
                alt = glob.glob(f"*_{rep}.fasta")
                if alt:
                    fasta_candidate = alt[0]
                else:
                    print(f"[Warning] FASTA not found for replicate {rep} (tried {fasta_candidate})")
                    continue

            stats = compute_sss_for_fasta(fasta_candidate, aligner)
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
            # Try to infer replicate number from filename
            fname = os.path.basename(fpath)
            rep = 1
            parts = fname.replace(".fasta", "").replace(".fa", "").split("_")
            for p in reversed(parts):
                if p.isdigit():
                    rep = int(p)
                    break
            stats = compute_sss_for_fasta(fpath, aligner)
            records.append({
                "distance": args.distance,
                "length": args.length,
                "replicate": rep,
                **stats
            })

    # Case 3: Single fasta
    elif args.fasta:
        stats = compute_sss_for_fasta(args.fasta, aligner)
        records.append({
            "distance": args.distance,
            "length": args.length,
            "replicate": 1,
            **stats
        })

    if not records:
        print("[Warning] No SSS records computed.")
        # Create empty DataFrame with expected headers
        df = pd.DataFrame(columns=[
            "distance", "length", "replicate", "num_taxa",
            "sss_mean", "sss_median", "sss_std", "sss_min", "sss_max", "sss_zero_frac"
        ])
    else:
        df = pd.DataFrame(records)

    header = not os.path.exists(args.outcsv)
    os.makedirs(os.path.dirname(os.path.abspath(args.outcsv)), exist_ok=True)
    df.to_csv(args.outcsv, mode="a", index=False, header=header)
    print(f"Computed SSS for {len(records)} replicates -> {args.outcsv}")


if __name__ == "__main__":
    main()
