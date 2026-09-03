#!/usr/bin/env python3
"""
Run GS (Graph Splitting Method v2.5 by Motomu Matsui & Wataru Iwasaki)
Performs alignment-free / pairwise graph-based phylogenetic tree reconstruction
following the benchmark protocol of Matsui & Iwasaki (2020, Systematic Biology).

GS settings:
  gs2 -s -l -t <cpus> -m 7.5 <fasta> > <outtree>
"""

import os
import sys
import shutil
import argparse
import subprocess
import time
import dendropy

def find_gs_binary() -> str:
    """Locate the gs2 executable."""
    script_dir = os.path.dirname(os.path.abspath(__file__))
    candidates = [
        os.path.join(script_dir, "gs2"),
        shutil.which("gs2"),
        os.path.expanduser("~/bin/gs2"),
    ]
    for c in candidates:
        if c and os.path.isfile(c) and os.access(c, os.X_OK):
            return c
    
    # Try running install_gs.sh if present
    install_script = os.path.join(script_dir, "install_gs.sh")
    if os.path.exists(install_script):
        print("[GS] Attempting to auto-build gs2 via install_gs.sh...", file=sys.stderr)
        subprocess.run(["bash", install_script], check=True)
        gs_bin = os.path.join(script_dir, "gs2")
        if os.path.isfile(gs_bin) and os.access(gs_bin, os.X_OK):
            return gs_bin

    raise FileNotFoundError("gs2 binary not found. Please run bin/install_gs.sh to build it.")

def main():
    parser = argparse.ArgumentParser(description="Run Graph Splitting (GS) Method (syz049 benchmark)")
    parser.add_argument("--fasta", required=True, help="Input unaligned FASTA file")
    parser.add_argument("--outtree", required=True, help="Output Newick tree file (.nwk)")
    parser.add_argument("--threads", type=int, default=2, help="Number of CPU threads for MMseqs2 (default: 2)")
    parser.add_argument("--sensitivity", type=float, default=7.5, help="MMseqs2 sensitivity (default: 7.5)")
    args = parser.parse_args()

    gs_bin = find_gs_binary()
    
    # Verify mmseqs is available
    if not shutil.which("mmseqs"):
        raise EnvironmentError("mmseqs executable not found in PATH. Please install mmseqs2 (micromamba install -c bioconda mmseqs2).")

    outdir = os.path.dirname(os.path.abspath(args.outtree)) or "."
    fasta_path = os.path.abspath(args.fasta)
    
    t0 = time.time()
    cmd = [
        gs_bin,
        "-s",
        "-l",
        "-t", str(args.threads),
        "-m", str(args.sensitivity),
        fasta_path
    ]
    
    proc = subprocess.run(cmd, cwd=outdir, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"gs2 execution failed (exit code {proc.returncode}):\n{proc.stderr}")

    raw_newick = proc.stdout.strip()
    if not raw_newick or "(" not in raw_newick:
        raise ValueError(f"Invalid Newick output from gs2: {raw_newick}\nStderr: {proc.stderr}")

    # Remove double quotes around taxon names to conform to standard Newick (DendroPy requirement)
    cleaned_newick = raw_newick.replace('"', '')
    if not cleaned_newick.endswith(";"):
        cleaned_newick += ";"

    # Verify validity via DendroPy
    try:
        dendropy.Tree.get(string=cleaned_newick, schema="newick")
    except Exception as e:
        raise ValueError(f"Failed to parse cleaned GS Newick tree with DendroPy: {e}\nTree string: {cleaned_newick}")

    with open(args.outtree, "w") as f:
        f.write(cleaned_newick + "\n")

    elapsed = time.time() - t0
    print(f"[GS] Successfully reconstructed tree in {elapsed:.2f}s -> {args.outtree}")

if __name__ == "__main__":
    main()
