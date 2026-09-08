#!/usr/bin/env python3
"""
PWA + NJ Pipeline Script
Performs pairwise Needleman-Wunsch alignment between amino acid sequences using BLOSUM62,
calculates gap-eliminated Poisson evolutionary distances (with Gamma-Poisson distance option),
and infers Neighbor-Joining (NJ) trees using RapidNJ or FastME.
"""

import sys
import os
import math
import argparse
import subprocess
import numpy as np
from Bio import SeqIO
from Bio.Align import PairwiseAligner, substitution_matrices

# Load standard BLOSUM62 substitution matrix
try:
    _BLOSUM62 = substitution_matrices.load("BLOSUM62")
except Exception as e:
    raise RuntimeError(f"Failed to load BLOSUM62 matrix from Biopython: {e}")

def create_aligner(gap_open=10.0, gap_extend=0.5):
    """
    Creates and configures a C-accelerated Bio.Align.PairwiseAligner.
    Uses affine gap scoring matching standard Needleman-Wunsch:
    Gap penalty = gap_open + k * gap_extend
    """
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    aligner.substitution_matrix = _BLOSUM62
    aligner.open_gap_score = -abs(gap_open)
    aligner.extend_gap_score = -abs(gap_extend)
    return aligner

_THREAD_LOCAL_ALIGNER = None

def get_aligner(gap_open=10.0, gap_extend=0.5):
    """Retrieves or creates a cached PairwiseAligner instance."""
    global _THREAD_LOCAL_ALIGNER
    if _THREAD_LOCAL_ALIGNER is None:
        _THREAD_LOCAL_ALIGNER = create_aligner(gap_open, gap_extend)
    return _THREAD_LOCAL_ALIGNER

def needleman_wunsch(seq1, seq2, aligner=None, gap_open=10.0, gap_extend=0.5):
    """
    Needleman-Wunsch global alignment using Biopython's C-accelerated PairwiseAligner.
    Returns (aligned_seq1, aligned_seq2).
    """
    if aligner is None:
        aligner = get_aligner(gap_open=gap_open, gap_extend=gap_extend)

    alignments = aligner.align(seq1, seq2)
    best_alignment = alignments[0]
    return str(best_alignment[0]), str(best_alignment[1])

def compute_pairwise_p_distance(
    aligned1: str | np.ndarray,
    aligned2: str | np.ndarray,
    min_overlap: int = 10
) -> float:
    """
    Calculates proportion of amino acid mismatches p excluding gap sites.

    Parameters
    ----------
    aligned1 : str or numpy.ndarray
        First aligned protein sequence.
    aligned2 : str or numpy.ndarray
        Second aligned protein sequence.
    min_overlap : int, default=10
        Minimum number of overlapping non-gap residue sites required.
        If valid overlapping sites < min_overlap, returns fallback max proportion (0.94)
        to prevent extreme distance fluctuations due to tiny overlaps.

    Returns
    -------
    float
        Mismatch proportion p (unitless, in range [0.0, 0.94]).
    """
    if isinstance(aligned1, str):
        arr1 = np.frombuffer(aligned1.upper().encode('ascii'), dtype=np.uint8)
        arr2 = np.frombuffer(aligned2.upper().encode('ascii'), dtype=np.uint8)
    else:
        arr1 = np.asarray(aligned1, dtype=np.uint8)
        arr2 = np.asarray(aligned2, dtype=np.uint8)

    gap_byte = 45  # ASCII code for '-'
    valid_mask = (arr1 != gap_byte) & (arr2 != gap_byte)
    valid_sites = int(np.count_nonzero(valid_mask))

    # Guard against zero or insufficient overlap (less than min_overlap residues)
    if valid_sites < min_overlap:
        return 0.94  # Max proportion fallback for insufficient overlap

    mismatches = int(np.count_nonzero(arr1[valid_mask] != arr2[valid_mask]))
    p = mismatches / valid_sites
    return min(p, 0.94)

def calculate_poisson_distance(
    aligned1: str | np.ndarray,
    aligned2: str | np.ndarray,
    min_overlap: int = 10
) -> float:
    """
    Computes standard Poisson evolutionary distance d from aligned sequence pair.

    Parameters
    ----------
    aligned1 : str or numpy.ndarray
        First aligned sequence.
    aligned2 : str or numpy.ndarray
        Second aligned sequence.
    min_overlap : int, default=10
        Minimum required overlapping residue sites.

    Returns
    -------
    float
        Estimated evolutionary distance d (expected substitutions per site, unitless).
    """
    p = compute_pairwise_p_distance(aligned1, aligned2, min_overlap=min_overlap)
    term = 1.0 - (20.0 / 19.0) * p
    if term <= 0.001:
        term = 0.001
    d = - (19.0 / 20.0) * math.log(term)
    return max(0.0, d)

def calculate_gamma_poisson_distance(
    aligned1: str | np.ndarray,
    aligned2: str | np.ndarray,
    alpha: float = 1.0,
    min_overlap: int = 10
) -> float:
    """
    Computes Gamma-corrected Poisson evolutionary distance d from aligned sequence pair.

    Parameters
    ----------
    aligned1 : str or numpy.ndarray
        First aligned sequence.
    aligned2 : str or numpy.ndarray
        Second aligned sequence.
    alpha : float, default=1.0
        Gamma shape parameter alpha (unitless).
    min_overlap : int, default=10
        Minimum required overlapping residue sites.

    Returns
    -------
    float
        Estimated evolutionary distance d (expected substitutions per site, unitless).
    """
    if alpha is None or alpha <= 0:
        return calculate_poisson_distance(aligned1, aligned2, min_overlap=min_overlap)

    p = compute_pairwise_p_distance(aligned1, aligned2, min_overlap=min_overlap)
    term = 1.0 - (20.0 / 19.0) * p
    if term <= 0.001:
        term = 0.001
    d = (19.0 / 20.0) * alpha * (math.pow(term, -1.0 / alpha) - 1.0)
    return max(0.0, d)

def calculate_distance(
    aligned1: str | np.ndarray,
    aligned2: str | np.ndarray,
    dist_model: str = "poisson",
    alpha: float = 1.0,
    min_overlap: int = 10
) -> float:
    """
    Dispatches to either Poisson or Gamma-Poisson distance calculation.

    Parameters
    ----------
    aligned1 : str or numpy.ndarray
        First aligned sequence.
    aligned2 : str or numpy.ndarray
        Second aligned sequence.
    dist_model : str, default="poisson"
        Distance model ('poisson' or 'gamma_poisson').
    alpha : float, default=1.0
        Gamma shape parameter alpha (used when dist_model='gamma_poisson').
    min_overlap : int, default=10
        Minimum required overlapping residue sites.

    Returns
    -------
    float
        Estimated evolutionary distance d (unitless).
    """
    if dist_model == "gamma_poisson":
        return calculate_gamma_poisson_distance(aligned1, aligned2, alpha=alpha, min_overlap=min_overlap)
    return calculate_poisson_distance(aligned1, aligned2, min_overlap=min_overlap)

def run_nj_tool(matrix_file, outtree_file, tool="rapidnj"):
    """
    Runs RapidNJ or FastME to infer NJ tree.
    Strictly raises RuntimeError if tool is unavailable or fails.
    """
    if tool == "rapidnj":
        try:
            res = subprocess.run(
                ["rapidnj", matrix_file, "-i", "pd", "-x", outtree_file],
                capture_output=True, text=True
            )
            if res.returncode != 0 or not os.path.exists(outtree_file) or os.path.getsize(outtree_file) == 0:
                raise RuntimeError(f"RapidNJ execution failed (exit code {res.returncode}):\nSTDOUT: {res.stdout}\nSTDERR: {res.stderr}")
            return
        except FileNotFoundError:
            raise RuntimeError("RapidNJ executable 'rapidnj' was not found in PATH. Please ensure rapidnj is installed.")

    elif tool == "fastme":
        try:
            res = subprocess.run(
                ["fastme", "-i", matrix_file, "-o", outtree_file],
                capture_output=True, text=True
            )
            if res.returncode != 0 or not os.path.exists(outtree_file) or os.path.getsize(outtree_file) == 0:
                raise RuntimeError(f"FastME execution failed (exit code {res.returncode}):\nSTDOUT: {res.stdout}\nSTDERR: {res.stderr}")
            return
        except FileNotFoundError:
            raise RuntimeError("FastME executable 'fastme' was not found in PATH. Please ensure fastme is installed.")

    else:
        raise ValueError(f"Unsupported NJ tool: {tool}")

def _compute_single_pair(args_tuple):
    """Worker function for parallel pairwise distance computation."""
    i, j, s1, s2, gap_open, gap_extend, dist_model, alpha = args_tuple
    al1, al2 = needleman_wunsch(s1, s2, gap_open=gap_open, gap_extend=gap_extend)
    d = calculate_distance(al1, al2, dist_model=dist_model, alpha=alpha)
    return (i, j, d)

def main():
    parser = argparse.ArgumentParser(description="PWA+NJ Pipeline Execution")
    parser.add_argument("--fasta", help="Input FASTA sequence file")
    parser.add_argument("--matrix", help="Input precomputed PHYLIP distance matrix (skips alignment if provided)")
    parser.add_argument("--outtree", required=True, help="Output Newick tree file")
    parser.add_argument("--outmatrix", help="Output PHYLIP distance matrix file")
    parser.add_argument("--gap_open", type=float, default=10.0, help="Gap open penalty (default: 10.0)")
    parser.add_argument("--gap_extend", type=float, default=0.5, help="Gap extension penalty (default: 0.5)")
    parser.add_argument("--dist_model", choices=["poisson", "gamma_poisson"], default="poisson", help="Distance model (default: poisson)")
    parser.add_argument("--alpha", type=float, default=1.0, help="Gamma shape parameter alpha for Gamma-Poisson distance (default: 1.0)")
    parser.add_argument("--tool", choices=["rapidnj", "fastme"], default="rapidnj", help="NJ software to use (default: rapidnj)")
    parser.add_argument("--threads", type=int, default=1, help="Number of parallel CPU workers for pairwise alignments (default: 1)")
    args = parser.parse_args()

    if args.matrix:
        # Precomputed matrix supplied directly
        run_nj_tool(args.matrix, args.outtree, tool=args.tool)
        print(f"Tree successfully inferred from {args.matrix} using {args.tool} -> {args.outtree}")
        return

    if not args.fasta:
        parser.error("Either --fasta or --matrix must be provided.")

    records = list(SeqIO.parse(args.fasta, "fasta"))
    names = [rec.id for rec in records]
    seqs = [str(rec.seq) for rec in records]
    N = len(names)

    dist_matrix = [[0.0] * N for _ in range(N)]

    if args.threads > 1 and N >= 16:
        # Multi-process parallel computation for large taxa counts
        pair_tasks = [
            (i, j, seqs[i], seqs[j], args.gap_open, args.gap_extend, args.dist_model, args.alpha)
            for i in range(N)
            for j in range(i + 1, N)
        ]
        from concurrent.futures import ProcessPoolExecutor
        with ProcessPoolExecutor(max_workers=args.threads) as executor:
            results = executor.map(_compute_single_pair, pair_tasks, chunksize=max(1, len(pair_tasks) // (args.threads * 4)))
            for i, j, d in results:
                dist_matrix[i][j] = d
                dist_matrix[j][i] = d
    else:
        # Direct C-accelerated single-process computation (avoids IPC/serialization overhead)
        aligner = get_aligner(gap_open=args.gap_open, gap_extend=args.gap_extend)
        for i in range(N):
            s1 = seqs[i]
            for j in range(i + 1, N):
                s2 = seqs[j]
                al1, al2 = needleman_wunsch(s1, s2, aligner=aligner, gap_open=args.gap_open, gap_extend=args.gap_extend)
                d = calculate_distance(al1, al2, dist_model=args.dist_model, alpha=args.alpha)
                dist_matrix[i][j] = d
                dist_matrix[j][i] = d

    # Write PHYLIP distance matrix
    matrix_str = f"   {N}\n"
    for i in range(N):
        row = f"{names[i]:<10}" + "".join(f"  {dist_matrix[i][j]:.6f}" for j in range(N)) + "\n"
        matrix_str += row

    if args.outmatrix:
        with open(args.outmatrix, "w") as f:
            f.write(matrix_str)

    tmp_matrix_file = args.outmatrix or (args.outtree + ".phylip")
    if not args.outmatrix:
        with open(tmp_matrix_file, "w") as f:
            f.write(matrix_str)

    # Run NJ tool (RapidNJ / FastME)
    run_nj_tool(tmp_matrix_file, args.outtree, tool=args.tool)

    if not args.outmatrix and os.path.exists(tmp_matrix_file):
        try:
            os.remove(tmp_matrix_file)
        except OSError:
            pass

    print(f"PWA+NJ tree successfully written to {args.outtree} (dist_model={args.dist_model}, gap_open={args.gap_open}, gap_extend={args.gap_extend}, tool={args.tool})")

if __name__ == "__main__":
    main()
