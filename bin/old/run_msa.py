#!/usr/bin/env python3
"""
Standalone Multiple Sequence Alignment (MSA) Script
Runs MAFFT to align unaligned amino acid sequences.
Allows MSA results to be shared and reused across downstream pipelines (MSA+NJ, MSA+ML).
"""

import sys
import os
import argparse
import subprocess

def run_mafft(fasta_file, outmsa_file, threads=1, auto=True):
    """Executes MAFFT to generate multiple sequence alignment with --threadit 0."""
    cmd = ["mafft", "--threadit", "0"]
    if auto:
        cmd.append("--auto")
    if threads > 1:
        cmd.extend(["--thread", str(threads)])
    cmd.append("--quiet")
    cmd.append(fasta_file)

    try:
        with open(outmsa_file, "w") as out_f:
            res = subprocess.run(cmd, stdout=out_f, stderr=subprocess.PIPE, text=True)
            if res.returncode != 0 or not os.path.exists(outmsa_file) or os.path.getsize(outmsa_file) == 0:
                raise RuntimeError(f"MAFFT alignment failed (exit code {res.returncode}):\n{res.stderr}")
    except FileNotFoundError:
        raise RuntimeError("MAFFT executable 'mafft' was not found in PATH. Please ensure MAFFT is installed.")

def main():
    parser = argparse.ArgumentParser(description="Multiple Sequence Alignment Execution (MAFFT)")
    parser.add_argument("--fasta", required=True, help="Input unaligned FASTA file")
    parser.add_argument("--outmsa", required=True, help="Output aligned FASTA file")
    parser.add_argument("--threads", type=int, default=1, help="Number of CPU threads for MAFFT (default: 1)")
    args = parser.parse_args()

    run_mafft(args.fasta, args.outmsa, threads=args.threads)
    print(f"MSA successfully written to {args.outmsa} (tool=mafft, threads={args.threads})")

if __name__ == "__main__":
    main()
