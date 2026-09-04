#!/usr/bin/env python3
"""
MSA + RAxML Pipeline Script (-f d)

Performs Maximum Likelihood tree reconstruction using RAxML (v8.2.x) with the
rapid hill-climbing algorithm (-f d), matching the benchmark protocol of
Matsui & Iwasaki (2020, Systematic Biology).
"""

from typing import Dict, Any, List, Optional
import os
import sys
import shutil
import argparse
import subprocess
import time
import json
import tempfile
from Bio import SeqIO


def run_mafft(in_fasta: str, out_msa: str) -> bool:
    """
    Runs MAFFT alignment with --auto and --threadit 0.

    Parameters
    ----------
    in_fasta : str
        Input unaligned FASTA file path.
    out_msa : str
        Output aligned FASTA file path.

    Returns
    -------
    bool
        True if alignment succeeded, False otherwise.
    """
    try:
        cmd = ["mafft", "--threadit", "0", "--auto", in_fasta]
        with open(out_msa, "w") as out_f:
            res = subprocess.run(cmd, stdout=out_f, stderr=subprocess.PIPE, text=True)
        if res.returncode == 0 and os.path.exists(out_msa) and os.path.getsize(out_msa) > 0:
            return True
    except FileNotFoundError:
        pass
    return False


def find_raxml_binary() -> str:
    """
    Locates the RAxML executable in conda environment or system PATH.

    Returns
    -------
    str
        Path to raxml executable.

    Raises
    ------
    RuntimeError
        If no RAxML executable is found.
    """
    py_dir = os.path.dirname(sys.executable)
    candidates = [
        "raxmlHPC-PTHREADS",
        "raxmlHPC-AVX2",
        "raxmlHPC",
        "raxml",
    ]
    # Check current python environment first
    for name in candidates:
        full_path = os.path.join(py_dir, name)
        if os.path.isfile(full_path) and os.access(full_path, os.X_OK):
            return full_path

    # Check system PATH
    for name in candidates:
        path = shutil.which(name)
        if path:
            return path

    raise RuntimeError(
        "RAxML executable ('raxmlHPC', 'raxmlHPC-PTHREADS', or 'raxml') not found in PATH or environment."
    )


def create_star_tree(fasta_path: str, out_tree_path: str) -> bool:
    """
    Generates a trivial star tree fallback if all sequences are completely identical.

    Parameters
    ----------
    fasta_path : str
        Input FASTA file path.
    out_tree_path : str
        Output Newick tree path.

    Returns
    -------
    bool
        True if star tree was created, False otherwise.
    """
    records = list(SeqIO.parse(fasta_path, "fasta"))
    if not records:
        return False
    taxa_str = ",".join([f"{r.id}:0.001" for r in records])
    star_newick = f"({taxa_str});\n"
    with open(out_tree_path, "w") as f:
        f.write(star_newick)
    return True


def run_raxml(
    msa_file: str,
    out_tree: str,
    model: str = "PROTGAMMALGX",
    seed: int = 12345,
    threads: int = 1,
) -> Dict[str, Any]:
    """
    Executes RAxML -f d tree inference on the specified MSA.

    Parameters
    ----------
    msa_file : str
        Input aligned FASTA file path.
    out_tree : str
        Destination Newick tree file path.
    model : str, default="PROTGAMMALGX"
        Substitution model (e.g. PROTGAMMALGX or PROTGAMMAIWAGX).
    seed : int, default=12345
        Random seed for parsimony starting tree.
    threads : int, default=1
        Number of CPU threads.

    Returns
    -------
    Dict[str, Any]
        Metadata dictionary with runtime and model info.
    """
    raxml_bin = find_raxml_binary()
    run_id = f"raxml_run_{os.getpid()}_{int(time.time() * 1000) % 1000000}"
    abs_msa = os.path.abspath(msa_file)
    abs_out_tree = os.path.abspath(out_tree)

    with tempfile.TemporaryDirectory(prefix="raxml_work_") as tmpdir:
        abs_tmpdir = os.path.abspath(tmpdir)
        cmd: List[str] = [
            raxml_bin,
            "-m", model,
            "-f", "d",
            "-p", str(seed),
            "-s", abs_msa,
            "-n", run_id,
            "-w", abs_tmpdir,
        ]
        if threads > 1 and "PTHREADS" in os.path.basename(raxml_bin):
            cmd.extend(["-T", str(threads)])

        t0 = time.time()
        res = subprocess.run(cmd, capture_output=True, text=True)
        elapsed = time.time() - t0

        best_tree_file = os.path.join(abs_tmpdir, f"RAxML_bestTree.{run_id}")
        if not os.path.exists(best_tree_file) or os.path.getsize(best_tree_file) == 0:
            # Check if all sequences are identical
            if create_star_tree(abs_msa, abs_out_tree):
                print(f"[RAxML] Identical sequences detected, fallback to star tree -> {out_tree}")
                return {
                    "tool": "RAxML_star_fallback",
                    "model": model,
                    "algorithm": "-f d",
                    "seed": seed,
                    "elapsed_sec": elapsed,
                    "status": "star_fallback",
                }
            raise RuntimeError(
                f"RAxML execution failed (code {res.returncode}):\nSTDOUT: {res.stdout}\nSTDERR: {res.stderr}"
            )

        # Copy best tree to target destination
        os.makedirs(os.path.dirname(abs_out_tree) or ".", exist_ok=True)
        shutil.copyfile(best_tree_file, abs_out_tree)

        metadata = {
            "tool": "RAxML",
            "version": os.path.basename(raxml_bin),
            "model": model,
            "algorithm": "-f d",
            "seed": seed,
            "elapsed_sec": elapsed,
            "status": "success",
        }
        return metadata


def main() -> None:
    parser = argparse.ArgumentParser(description="MSA + RAxML Pipeline Execution (-f d)")
    parser.add_argument("--fasta", help="Input unaligned FASTA sequence file")
    parser.add_argument("--msa", help="Input pre-aligned MSA FASTA file (skips MAFFT if provided)")
    parser.add_argument("--outtree", required=True, help="Output Newick tree file")
    parser.add_argument("--outmsa", help="Output MAFFT MSA file (if MAFFT is run)")
    parser.add_argument("--outjson", help="Output JSON metadata file")
    parser.add_argument("--model", default="PROTGAMMALGX", help="RAxML amino acid model (default: PROTGAMMALGX)")
    parser.add_argument("--seed", type=int, default=12345, help="Random seed for RAxML parsimony tree (default: 12345)")
    parser.add_argument("--threads", type=int, default=1, help="Number of CPU threads for RAxML")
    args = parser.parse_args()

    if not args.msa and not args.fasta:
        parser.error("Either --msa (pre-aligned) or --fasta (unaligned) must be provided.")

    tmp_msa = args.msa
    created_tmp_msa = False
    if not tmp_msa:
        tmp_msa = args.outmsa or (args.outtree + ".msa.fasta")
        created_tmp_msa = True
        mafft_ok = run_mafft(args.fasta, tmp_msa)
        if not mafft_ok:
            # Fallback copy
            records = list(SeqIO.parse(args.fasta, "fasta"))
            with open(tmp_msa, "w") as f:
                for rec in records:
                    f.write(f">{rec.id}\n{str(rec.seq)}\n")

    metadata = run_raxml(
        msa_file=tmp_msa,
        out_tree=args.outtree,
        model=args.model,
        seed=args.seed,
        threads=args.threads,
    )

    json_out = args.outjson or (args.outtree + ".json")
    with open(json_out, "w") as f:
        json.dump(metadata, f, indent=2)

    # Clean up intermediate MSA if it was temporary
    if created_tmp_msa and not args.outmsa and os.path.exists(tmp_msa):
        try:
            os.remove(tmp_msa)
        except OSError:
            pass

    print(f"RAxML tree successfully written to {args.outtree} (model={args.model}, time={metadata.get('elapsed_sec', 0):.2f}s)")


if __name__ == "__main__":
    main()
