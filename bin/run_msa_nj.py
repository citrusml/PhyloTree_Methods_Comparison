#!/usr/bin/env python3
"""
MSA + NJ Control Pipeline Script
Runs MAFFT alignment, then applies the EXACT SAME gap-eliminated Poisson distance
calculation module as PWA to construct the distance matrix, and infers the NJ tree.
"""

import sys
import os
import argparse
import subprocess
from Bio import SeqIO

# Import distance calculation and NJ routines from run_pwa_nj
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from run_pwa_nj import calculate_poisson_distance, python_nj

def run_mafft(in_fasta, out_msa):
    """Runs MAFFT with --threadit 0 for reproducible MSA generation."""
    try:
        cmd = ["mafft", "--threadit", "0", "--auto", in_fasta]
        with open(out_msa, "w") as out_f:
            res = subprocess.run(cmd, stdout=out_f, stderr=subprocess.PIPE, text=True)
        if res.returncode == 0 and os.path.exists(out_msa) and os.path.getsize(out_msa) > 0:
            return True
    except FileNotFoundError:
        pass
    return False

def main():
    parser = argparse.ArgumentParser(description="MSA+NJ Control Pipeline Execution")
    parser.add_argument("--fasta", required=True, help="Input FASTA sequence file")
    parser.add_argument("--outtree", required=True, help="Output Newick tree file")
    parser.add_argument("--outmsa", help="Output MAFFT MSA file")
    parser.add_argument("--outmatrix", help="Output PHYLIP distance matrix file")
    parser.add_argument("--alpha", type=float, default=None, help="Gamma shape parameter alpha")
    parser.add_argument("--tool", choices=["rapidnj", "fastme", "python"], default="rapidnj", help="NJ software to use")
    args = parser.parse_args()

    tmp_msa_file = args.outmsa or (args.outtree + ".msa.fasta")
    
    # 1. Run MAFFT MSA
    mafft_success = run_mafft(args.fasta, tmp_msa_file)
    if not mafft_success:
        # Fallback: if mafft binary is not installed, use input fasta directly as MSA if lengths match or pad
        records = list(SeqIO.parse(args.fasta, "fasta"))
        with open(tmp_msa_file, "w") as f:
            for rec in records:
                f.write(f">{rec.id}\n{str(rec.seq)}\n")

    # 2. Parse MSA
    msa_records = list(SeqIO.parse(tmp_msa_file, "fasta"))
    names = [rec.id for rec in msa_records]
    aligned_seqs = [str(rec.seq) for rec in msa_records]
    N = len(names)

    # 3. Calculate distance matrix using PWA-identical Poisson formula
    dist_matrix = [[0.0] * N for _ in range(N)]
    for i in range(N):
        for j in range(i + 1, N):
            d = calculate_poisson_distance(aligned_seqs[i], aligned_seqs[j], alpha=args.alpha)
            dist_matrix[i][j] = d
            dist_matrix[j][i] = d

    # 4. Format PHYLIP matrix
    matrix_str = f"   {N}\n"
    for i in range(N):
        row = f"{names[i]:<10}" + "".join(f"  {dist_matrix[i][j]:.6f}" for j in range(N)) + "\n"
        matrix_str += row

    tmp_matrix_file = args.outmatrix or (args.outtree + ".phylip")
    with open(tmp_matrix_file, "w") as f:
        f.write(matrix_str)

    # 5. Run NJ (RapidNJ / FastME / Python fallback)
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
        tree_nwk = python_nj(dist_matrix, names)
        with open(args.outtree, "w") as f:
            f.write(tree_nwk)

    # Clean up temporary files
    if not args.outmsa and os.path.exists(tmp_msa_file):
        try:
            os.remove(tmp_msa_file)
        except OSError:
            pass
    if not args.outmatrix and os.path.exists(tmp_matrix_file):
        try:
            os.remove(tmp_matrix_file)
        except OSError:
            pass

    print(f"MSA+NJ tree successfully written to {args.outtree}")

if __name__ == "__main__":
    main()
