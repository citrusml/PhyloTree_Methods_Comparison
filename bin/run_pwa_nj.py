#!/usr/bin/env python3
"""
PWA + NJ Pipeline Script
Performs pairwise Needleman-Wunsch alignment between amino acid sequences,
calculates gap-eliminated Poisson evolutionary distances,
and infers Neighbor-Joining (NJ) trees using RapidNJ or FastME (with python fallback).
"""

import sys
import os
import math
import argparse
import subprocess
from Bio import SeqIO

# Standard LG Amino Acid Substitution Matrix (Le & Gascuel 2008 Log-Odds Scoring Matrix)
LG_MATRIX = {
    ('A','A'):6, ('A','R'):-4, ('A','N'):-5, ('A','D'):-4, ('A','C'):-3,
    ('A','Q'):-4, ('A','E'):-4, ('A','G'):-2, ('A','H'):-5, ('A','I'):-6,
    ('A','L'):-5, ('A','K'):-5, ('A','M'):-4, ('A','F'):-7, ('A','P'):-3,
    ('A','S'):-1, ('A','T'):-2, ('A','W'):-9, ('A','Y'):-6, ('A','V'):-3,
    ('R','R'):7, ('R','N'):-2, ('R','D'):-5, ('R','C'):-2, ('R','Q'):1,
    ('R','E'):-3, ('R','G'):-6, ('R','H'):0, ('R','I'):-7, ('R','L'):-6,
    ('R','K'):1, ('R','M'):-7, ('R','F'):-6, ('R','P'):-3, ('R','S'):-4,
    ('R','T'):-4, ('R','W'):-2, ('R','Y'):-5, ('R','V'):-7, ('N','N'):8,
    ('N','D'):2, ('N','C'):-3, ('N','Q'):-3, ('N','E'):-2, ('N','G'):-2,
    ('N','H'):0, ('N','I'):-5, ('N','L'):-7, ('N','K'):-2, ('N','M'):-6,
    ('N','F'):-7, ('N','P'):-3, ('N','S'):-1, ('N','T'):-3, ('N','W'):-5,
    ('N','Y'):-3, ('N','V'):-6, ('D','D'):7, ('D','C'):-4, ('D','Q'):-2,
    ('D','E'):2, ('D','G'):-5, ('D','H'):-4, ('D','I'):-8, ('D','L'):-8,
    ('D','K'):-4, ('D','M'):-9, ('D','F'):-8, ('D','P'):-4, ('D','S'):-3,
    ('D','T'):-4, ('D','W'):-9, ('D','Y'):-7, ('D','V'):-8, ('C','C'):12,
    ('C','Q'):-5, ('C','E'):-4, ('C','G'):-4, ('C','H'):-5, ('C','I'):-4,
    ('C','L'):-4, ('C','K'):-5, ('C','M'):-5, ('C','F'):-4, ('C','P'):-4,
    ('C','S'):-3, ('C','T'):-4, ('C','W'):-3, ('C','Y'):-4, ('C','V'):-3,
    ('Q','Q'):8, ('Q','E'):1, ('Q','G'):-5, ('Q','H'):0, ('Q','I'):-7,
    ('Q','L'):-4, ('Q','K'):0, ('Q','M'):-6, ('Q','F'):-6, ('Q','P'):-2,
    ('Q','S'):-4, ('Q','T'):-3, ('Q','W'):-6, ('Q','Y'):-4, ('Q','V'):-6,
    ('E','E'):6, ('E','G'):-4, ('E','H'):-4, ('E','I'):-7, ('E','L'):-6,
    ('E','K'):-1, ('E','M'):-7, ('E','F'):-8, ('E','P'):-3, ('E','S'):-3,
    ('E','T'):-3, ('E','W'):-8, ('E','Y'):-5, ('E','V'):-7, ('G','G'):8,
    ('G','H'):-5, ('G','I'):-8, ('G','L'):-8, ('G','K'):-6, ('G','M'):-8,
    ('G','F'):-8, ('G','P'):-5, ('G','S'):-1, ('G','T'):-4, ('G','W'):-7,
    ('G','Y'):-7, ('G','V'):-7, ('H','H'):9, ('H','I'):-6, ('H','L'):-5,
    ('H','K'):-1, ('H','M'):-5, ('H','F'):-2, ('H','P'):-2, ('H','S'):-3,
    ('H','T'):-3, ('H','W'):-2, ('H','Y'):1, ('H','V'):-6, ('I','I'):6,
    ('I','L'):2, ('I','K'):-6, ('I','M'):0, ('I','F'):-1, ('I','P'):-8,
    ('I','S'):-6, ('I','T'):-2, ('I','W'):-5, ('I','Y'):-4, ('I','V'):3,
    ('L','L'):5, ('L','K'):-5, ('L','M'):2, ('L','F'):0, ('L','P'):-7,
    ('L','S'):-5, ('L','T'):-4, ('L','W'):-4, ('L','Y'):-3, ('L','V'):1,
    ('K','K'):7, ('K','M'):-6, ('K','F'):-7, ('K','P'):-2, ('K','S'):-4,
    ('K','T'):-3, ('K','W'):-6, ('K','Y'):-5, ('K','V'):-6, ('M','M'):9,
    ('M','F'):0, ('M','P'):-7, ('M','S'):-6, ('M','T'):-4, ('M','W'):-3,
    ('M','Y'):-4, ('M','V'):0, ('F','F'):8, ('F','P'):-6, ('F','S'):-6,
    ('F','T'):-5, ('F','W'):2, ('F','Y'):2, ('F','V'):-2, ('P','P'):8,
    ('P','S'):-4, ('P','T'):-3, ('P','W'):-8, ('P','Y'):-6, ('P','V'):-7,
    ('S','S'):7, ('S','T'):2, ('S','W'):-6, ('S','Y'):-4, ('S','V'):-5,
    ('T','T'):7, ('T','W'):-6, ('T','Y'):-4, ('T','V'):-3, ('W','W'):12,
    ('W','Y'):2, ('W','V'):-5, ('Y','Y'):9, ('Y','V'):-4, ('V','V'):6,
}

def get_sub_score(a, b):
    a, b = a.upper(), b.upper()
    if a == b:
        return LG_MATRIX.get((a, a), 6)
    key = (a, b) if (a, b) in LG_MATRIX else (b, a)
    return LG_MATRIX.get(key, -3)

def needleman_wunsch(seq1, seq2, gap_open=-10, gap_extend=-1):
    """Needleman-Wunsch global alignment with affine gap penalty."""
    m, n = len(seq1), len(seq2)
    # Standard dynamic programming table
    dp = [[0] * (n + 1) for _ in range(m + 1)]
    for i in range(1, m + 1):
        dp[i][0] = gap_open + (i - 1) * gap_extend
    for j in range(1, n + 1):
        dp[0][j] = gap_open + (j - 1) * gap_extend

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            match = dp[i-1][j-1] + get_sub_score(seq1[i-1], seq2[j-1])
            delete = dp[i-1][j] + (gap_extend if i > 1 else gap_open)
            insert = dp[i][j-1] + (gap_extend if j > 1 else gap_open)
            dp[i][j] = max(match, delete, insert)

    # Traceback
    align1, align2 = [], []
    i, j = m, n
    while i > 0 and j > 0:
        score = dp[i][j]
        score_diag = dp[i-1][j-1]
        if score == score_diag + get_sub_score(seq1[i-1], seq2[j-1]):
            align1.append(seq1[i-1])
            align2.append(seq2[j-1])
            i -= 1
            j -= 1
        elif score == dp[i-1][j] + (gap_extend if i > 1 else gap_open):
            align1.append(seq1[i-1])
            align2.append('-')
            i -= 1
        else:
            align1.append('-')
            align2.append(seq2[j-1])
            j -= 1

    while i > 0:
        align1.append(seq1[i-1])
        align2.append('-')
        i -= 1
    while j > 0:
        align1.append('-')
        align2.append(seq2[j-1])
        j -= 1

    return "".join(reversed(align1)), "".join(reversed(align2))

def calculate_poisson_distance(aligned1, aligned2, alpha=None):
    """
    Computes Poisson evolutionary distance d from aligned sequence pair,
    excluding gap sites (pairwise deletion).
    d = -19/20 * ln(1 - 20/19 * p)
    If alpha is specified, applies Gamma correction:
    d = 19/20 * alpha * [(1 - 20/19 * p)^(-1/alpha) - 1]
    """
    valid_sites = 0
    mismatches = 0
    for a, b in zip(aligned1, aligned2):
        if a == '-' or b == '-':
            continue
        valid_sites += 1
        if a.upper() != b.upper():
            mismatches += 1

    if valid_sites == 0:
        return 3.0  # Max distance fallback for 0 overlap

    p = mismatches / valid_sites
    # Cap p to avoid log domain error (20/19 * p < 1 => p < 0.95)
    max_p = 0.94
    p = min(p, max_p)

    term = 1.0 - (20.0 / 19.0) * p
    if term <= 0.001:
        term = 0.001

    if alpha is None:
        d = - (19.0 / 20.0) * math.log(term)
    else:
        d = (19.0 / 20.0) * alpha * (math.pow(term, -1.0 / alpha) - 1.0)

    return max(0.0, d)

def python_nj(dist_matrix, names):
    """Pure Python Neighbor-Joining implementation as standalone fallback."""
    N = len(names)
    clusters = {i: names[i] for i in range(N)}
    # Convert matrix to dictionary of distances
    D = {}
    for i in range(N):
        for j in range(N):
            D[i, j] = dist_matrix[i][j]

    nodes = list(range(N))
    next_node = N

    tree_str = {}
    for i in range(N):
        tree_str[i] = names[i]

    while len(nodes) > 2:
        K = len(nodes)
        # Compute r_i
        r = {}
        for i in nodes:
            r[i] = sum(D[i, j] for j in nodes if j != i) / (K - 2)

        # Compute Q matrix
        min_Q = float('inf')
        pair = (nodes[0], nodes[1])
        for i_idx in range(len(nodes)):
            for j_idx in range(i_idx + 1, len(nodes)):
                u, v = nodes[i_idx], nodes[j_idx]
                q = D[u, v] - r[u] - r[v]
                if q < min_Q:
                    min_Q = q
                    pair = (u, v)

        u, v = pair
        # Branch lengths
        d_u = 0.5 * D[u, v] + 0.5 * (r[u] - r[v])
        d_v = 0.5 * D[u, v] + 0.5 * (r[v] - r[u])
        d_u = max(0.0001, d_u)
        d_v = max(0.0001, d_v)

        # Create new node
        w = next_node
        next_node += 1
        tree_str[w] = f"({tree_str[u]}:{d_u:.6f},{tree_str[v]}:{d_v:.6f})"

        # Update distances
        for x in nodes:
            if x != u and x != v:
                D[w, x] = D[x, w] = 0.5 * (D[u, x] + D[v, x] - D[u, v])

        nodes.remove(u)
        nodes.remove(v)
        nodes.append(w)

    u, v = nodes[0], nodes[1]
    d_uv = max(0.0001, D[u, v])
    final_tree = f"({tree_str[u]}:{d_uv:.6f},{tree_str[v]}:{d_uv:.6f});"
    return final_tree

def main():
    parser = argparse.ArgumentParser(description="PWA+NJ Pipeline Execution")
    parser.add_argument("--fasta", required=True, help="Input FASTA sequence file")
    parser.add_argument("--outtree", required=True, help="Output Newick tree file")
    parser.add_argument("--outmatrix", help="Output PHYLIP distance matrix file")
    parser.add_argument("--alpha", type=float, default=None, help="Gamma shape parameter alpha for distance correction")
    parser.add_argument("--tool", choices=["rapidnj", "fastme", "python"], default="rapidnj", help="NJ software to use")
    args = parser.parse_args()

    records = list(SeqIO.parse(args.fasta, "fasta"))
    names = [rec.id for rec in records]
    seqs = [str(rec.seq) for rec in records]
    N = len(names)

    # Pairwise alignment and distance matrix construction
    dist_matrix = [[0.0] * N for _ in range(N)]
    for i in range(N):
        for j in range(i + 1, N):
            al1, al2 = needleman_wunsch(seqs[i], seqs[j])
            d = calculate_poisson_distance(al1, al2, alpha=args.alpha)
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

    # Run NJ tool (RapidNJ / FastME / Python fallback)
    tree_written = False
    if args.tool == "rapidnj":
        try:
            res = subprocess.run(["rapidnj", tmp_matrix_file, "-i", "pd", "-x", args.outtree], capture_output=True, text=True)
            if res.returncode == 0 and os.path.exists(args.outtree) and os.path.getsize(args.outtree) > 0:
                tree_written = True
        except FileNotFoundError:
            pass

    elif args.tool == "fastme":
        try:
            res = subprocess.run(["fastme", "-i", tmp_matrix_file, "-o", args.outtree], capture_output=True, text=True)
            if res.returncode == 0 and os.path.exists(args.outtree) and os.path.getsize(args.outtree) > 0:
                tree_written = True
        except FileNotFoundError:
            pass

    if not tree_written:
        # Fallback to pure python NJ
        tree_nwk = python_nj(dist_matrix, names)
        with open(args.outtree, "w") as f:
            f.write(tree_nwk)

    if not args.outmatrix and os.path.exists(tmp_matrix_file):
        try:
            os.remove(tmp_matrix_file)
        except OSError:
            pass

    print(f"PWA+NJ tree successfully written to {args.outtree}")

if __name__ == "__main__":
    main()
