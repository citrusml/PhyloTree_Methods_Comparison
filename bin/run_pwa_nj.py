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
from Bio import SeqIO
from Bio.Align import substitution_matrices

# Load standard BLOSUM62 substitution matrix
try:
    _BLOSUM62 = substitution_matrices.load("BLOSUM62")
except Exception as e:
    raise RuntimeError(f"Failed to load BLOSUM62 matrix from Biopython: {e}")

def get_sub_score(a, b):
    """Retrieves substitution score for amino acid pair (a, b) from BLOSUM62."""
    a, b = a.upper(), b.upper()
    try:
        return _BLOSUM62[a, b]
    except KeyError:
        try:
            return _BLOSUM62[b, a]
        except KeyError:
            return 4.0 if a == b else -4.0

def needleman_wunsch(seq1, seq2, gap_open=-10.0, gap_extend=-0.5):
    """
    Needleman-Wunsch global alignment with affine gap penalty.
    Uses three DP tables (M, Ix, Iy) to distinguish gap-open from gap-extend.
    A gap of length k costs: gap_open + k * gap_extend
    Note: gap_open and gap_extend are assumed negative.
    """
    # Ensure penalties are negative
    gap_open = -abs(gap_open)
    gap_extend = -abs(gap_extend)

    m, n = len(seq1), len(seq2)
    NEG_INF = float('-inf')

    # M[i][j]: best score ending with seq1[i] and seq2[j] aligned (match/mismatch)
    # Ix[i][j]: best score ending with seq1[i] aligned to a gap (gap in seq2 = deletion)
    # Iy[i][j]: best score ending with seq2[j] aligned to a gap (gap in seq1 = insertion)
    M  = [[NEG_INF] * (n + 1) for _ in range(m + 1)]
    Ix = [[NEG_INF] * (n + 1) for _ in range(m + 1)]
    Iy = [[NEG_INF] * (n + 1) for _ in range(m + 1)]

    M[0][0] = 0.0
    for i in range(1, m + 1):
        Ix[i][0] = gap_open + i * gap_extend
    for j in range(1, n + 1):
        Iy[0][j] = gap_open + j * gap_extend

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            s = get_sub_score(seq1[i-1], seq2[j-1])
            M[i][j] = max(
                M[i-1][j-1]  + s,
                Ix[i-1][j-1] + s,
                Iy[i-1][j-1] + s
            )
            Ix[i][j] = max(
                M[i-1][j]  + gap_open + gap_extend,
                Ix[i-1][j] + gap_extend
            )
            Iy[i][j] = max(
                M[i][j-1]  + gap_open + gap_extend,
                Iy[i][j-1] + gap_extend
            )

    # Traceback: determine terminal state
    align1, align2 = [], []
    i, j = m, n

    # Find best terminal state
    terminal = max(
        (M[m][n],  'M'),
        (Ix[m][n], 'Ix'),
        (Iy[m][n], 'Iy'),
        key=lambda x: x[0]
    )
    state = terminal[1]

    while i > 0 or j > 0:
        if state == 'M':
            s = get_sub_score(seq1[i-1], seq2[j-1])
            if i > 0 and j > 0 and M[i][j] == M[i-1][j-1] + s:
                align1.append(seq1[i-1]); align2.append(seq2[j-1])
                i -= 1; j -= 1; state = 'M'
            elif i > 0 and j > 0 and M[i][j] == Ix[i-1][j-1] + s:
                align1.append(seq1[i-1]); align2.append(seq2[j-1])
                i -= 1; j -= 1; state = 'Ix'
            elif i > 0 and j > 0 and M[i][j] == Iy[i-1][j-1] + s:
                align1.append(seq1[i-1]); align2.append(seq2[j-1])
                i -= 1; j -= 1; state = 'Iy'
            else:
                break
        elif state == 'Ix':
            if i > 0 and Ix[i][j] == M[i-1][j] + gap_open + gap_extend:
                align1.append(seq1[i-1]); align2.append('-')
                i -= 1; state = 'M'
            elif i > 0 and Ix[i][j] == Ix[i-1][j] + gap_extend:
                align1.append(seq1[i-1]); align2.append('-')
                i -= 1; state = 'Ix'
            else:
                break
        elif state == 'Iy':
            if j > 0 and Iy[i][j] == M[i][j-1] + gap_open + gap_extend:
                align1.append('-'); align2.append(seq2[j-1])
                j -= 1; state = 'M'
            elif j > 0 and Iy[i][j] == Iy[i][j-1] + gap_extend:
                align1.append('-'); align2.append(seq2[j-1])
                j -= 1; state = 'Iy'
            else:
                break

    # Handle remaining unaligned residues (leading gaps)
    while i > 0:
        align1.append(seq1[i-1]); align2.append('-')
        i -= 1
    while j > 0:
        align1.append('-'); align2.append(seq2[j-1])
        j -= 1

    return "".join(reversed(align1)), "".join(reversed(align2))

def compute_pairwise_p_distance(aligned1, aligned2):
    """Calculates proportion of amino acid mismatches p excluding gap sites."""
    valid_sites = 0
    mismatches = 0
    for a, b in zip(aligned1, aligned2):
        if a == '-' or b == '-':
            continue
        valid_sites += 1
        if a.upper() != b.upper():
            mismatches += 1

    if valid_sites == 0:
        return 0.94  # Max proportion fallback for 0 overlap
    p = mismatches / valid_sites
    return min(p, 0.94)

def calculate_poisson_distance(aligned1, aligned2):
    """
    Computes standard Poisson evolutionary distance d from aligned sequence pair,
    excluding gap sites (pairwise deletion).
    Formula: d = -19/20 * ln(1 - 20/19 * p)
    """
    p = compute_pairwise_p_distance(aligned1, aligned2)
    term = 1.0 - (20.0 / 19.0) * p
    if term <= 0.001:
        term = 0.001
    d = - (19.0 / 20.0) * math.log(term)
    return max(0.0, d)

def calculate_gamma_poisson_distance(aligned1, aligned2, alpha=1.0):
    """
    Computes Gamma-corrected Poisson evolutionary distance d from aligned sequence pair.
    Formula: d = 19/20 * alpha * [(1 - 20/19 * p)^(-1/alpha) - 1]
    """
    if alpha is None or alpha <= 0:
        return calculate_poisson_distance(aligned1, aligned2)

    p = compute_pairwise_p_distance(aligned1, aligned2)
    term = 1.0 - (20.0 / 19.0) * p
    if term <= 0.001:
        term = 0.001
    d = (19.0 / 20.0) * alpha * (math.pow(term, -1.0 / alpha) - 1.0)
    return max(0.0, d)

def calculate_distance(aligned1, aligned2, dist_model="poisson", alpha=1.0):
    """Dispatches to either Poisson or Gamma-Poisson distance."""
    if dist_model == "gamma_poisson":
        return calculate_gamma_poisson_distance(aligned1, aligned2, alpha=alpha)
    return calculate_poisson_distance(aligned1, aligned2)

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

def main():
    parser = argparse.ArgumentParser(description="PWA+NJ Pipeline Execution")
    parser.add_argument("--fasta", required=True, help="Input FASTA sequence file")
    parser.add_argument("--outtree", required=True, help="Output Newick tree file")
    parser.add_argument("--outmatrix", help="Output PHYLIP distance matrix file")
    parser.add_argument("--gap_open", type=float, default=10.0, help="Gap open penalty (default: 10.0)")
    parser.add_argument("--gap_extend", type=float, default=0.5, help="Gap extension penalty (default: 0.5)")
    parser.add_argument("--dist_model", choices=["poisson", "gamma_poisson"], default="poisson", help="Distance model (default: poisson)")
    parser.add_argument("--alpha", type=float, default=1.0, help="Gamma shape parameter alpha for Gamma-Poisson distance (default: 1.0)")
    parser.add_argument("--tool", choices=["rapidnj", "fastme"], default="rapidnj", help="NJ software to use (default: rapidnj)")
    args = parser.parse_args()

    records = list(SeqIO.parse(args.fasta, "fasta"))
    names = [rec.id for rec in records]
    seqs = [str(rec.seq) for rec in records]
    N = len(names)

    # Pairwise alignment and distance matrix construction
    dist_matrix = [[0.0] * N for _ in range(N)]
    for i in range(N):
        for j in range(i + 1, N):
            al1, al2 = needleman_wunsch(seqs[i], seqs[j], gap_open=args.gap_open, gap_extend=args.gap_extend)
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
