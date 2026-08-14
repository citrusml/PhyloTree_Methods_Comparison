#!/usr/bin/env python3
"""
MSA + NJ Control Pipeline Script
Runs MAFFT alignment, then applies the EXACT SAME gap-eliminated Poisson distance
calculation module as PWA to construct the distance matrix, and infers the NJ tree using RapidNJ or FastME.
"""

import sys
import os
import argparse
import subprocess
from Bio import SeqIO

# Import distance calculation and NJ execution from run_pwa_nj
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from run_pwa_nj import calculate_distance, run_nj_tool

def run_mafft(in_fasta, out_msa):
    """Runs MAFFT with --threadit 0 for reproducible MSA generation."""
    try:
        cmd = ["mafft", "--threadit", "0", "--auto", in_fasta]
        with open(out_msa, "w") as out_f:
            res = subprocess.run(cmd, stdout=out_f, stderr=subprocess.PIPE, text=True)
        if res.returncode != 0 or not os.path.exists(out_msa) or os.path.getsize(out_msa) == 0:
            raise RuntimeError(f"MAFFT alignment failed (exit code {res.returncode}):\n{res.stderr}")
    except FileNotFoundError:
        raise RuntimeError("MAFFT executable 'mafft' was not found in PATH. Please ensure mafft is installed.")

def main():
    parser = argparse.ArgumentParser(description="MSA+NJ Control Pipeline Execution")
    parser.add_argument("--fasta", required=True, help="Input FASTA sequence file")
    parser.add_argument("--outtree", required=True, help="Output Newick tree file")
    parser.add_argument("--outmsa", help="Output MAFFT MSA file")
    parser.add_argument("--outmatrix", help="Output PHYLIP distance matrix file")
    parser.add_argument("--dist_model", choices=["poisson", "gamma_poisson"], default="poisson", help="Distance model (default: poisson)")
    parser.add_argument("--alpha", type=float, default=1.0, help="Gamma shape parameter alpha (default: 1.0)")
    parser.add_argument("--tool", choices=["rapidnj", "fastme"], default="rapidnj", help="NJ software to use (default: rapidnj)")
    args = parser.parse_args()

    tmp_msa_file = args.outmsa or (args.outtree + ".msa.fasta")
    
    # 1. Run MAFFT MSA
    run_mafft(args.fasta, tmp_msa_file)

    # 2. Parse MSA
    msa_records = list(SeqIO.parse(tmp_msa_file, "fasta"))
    names = [rec.id for rec in msa_records]
    aligned_seqs = [str(rec.seq) for rec in msa_records]
    N = len(names)

    # 3. Calculate distance matrix using PWA-identical Poisson formula
    dist_matrix = [[0.0] * N for _ in range(N)]
    for i in range(N):
        for j in range(i + 1, N):
            d = calculate_distance(aligned_seqs[i], aligned_seqs[j], dist_model=args.dist_model, alpha=args.alpha)
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

    # 5. Run NJ (RapidNJ / FastME)
    run_nj_tool(tmp_matrix_file, args.outtree, tool=args.tool)

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

    print(f"MSA+NJ tree successfully written to {args.outtree} (dist_model={args.dist_model}, tool={args.tool})")

if __name__ == "__main__":
    main()
