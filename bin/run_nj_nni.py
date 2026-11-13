#!/usr/bin/env python3
"""
Refine a Neighbor-Joining (NJ) phylogenetic tree using FastME NNI
(Nearest Neighbor Interchange) from a PHYLIP distance matrix.
"""

import argparse
import os
import shutil
import subprocess
import sys
import tempfile
from typing import Optional


def find_fastme() -> str:
    """
    Locate the FastME executable in PATH or the active Python environment.

    Returns
    -------
    str
        Absolute path or command name for the FastME binary.

    Raises
    ------
    FileNotFoundError
        If the `fastme` executable cannot be found.
    """
    exe = shutil.which("fastme")
    if exe:
        return exe
    env_candidate = os.path.join(os.path.dirname(sys.executable), "fastme")
    if os.path.isfile(env_candidate) and os.access(env_candidate, os.X_OK):
        return env_candidate
    raise FileNotFoundError("Executable 'fastme' not found in PATH or active environment.")


def _sanitize_newick_for_fastme(nj_tree_path: str, tmp_dir: str) -> str:
    """
    Strip single quotes from taxon labels in a RapidNJ/Newick tree so FastME `-u`
    can parse it cleanly.

    Parameters
    ----------
    nj_tree_path : str
        Path to the starting NJ Newick file.
    tmp_dir : str
        Temporary directory to store the sanitized Newick file.

    Returns
    -------
    str
        Path to the sanitized Newick file.
    """
    with open(nj_tree_path, "r", encoding="utf-8") as fh:
        nwk = fh.read().strip().replace("'", "")
    sanitized_path = os.path.join(tmp_dir, "start_nj.nwk")
    with open(sanitized_path, "w", encoding="utf-8") as fh:
        fh.write(nwk + "\n")
    return sanitized_path


def run_fastme_nni(
    matrix_file: str,
    outtree_file: str,
    nj_tree_file: Optional[str] = None,
    nni_type: str = "B",
) -> None:
    """
    Run FastME NNI topology refinement starting from an existing NJ tree
    (or building an NJ tree from the distance matrix if omitted/incompatible).

    Parameters
    ----------
    matrix_file : str
        Path to the input PHYLIP distance matrix (pairwise distances, unitless).
    outtree_file : str
        Path to write the NNI-refined Newick tree (branch lengths in substitutions/site).
    nj_tree_file : Optional[str], optional
        Path to an existing starting NJ Newick tree (default: None).
    nni_type : str, optional
        FastME NNI mode: 'B' (BalME, default) or 'O' (OLS).

    Returns
    -------
    None

    Raises
    ------
    RuntimeError
        If FastME fails or produces an empty output tree.
    """
    fastme_exe = find_fastme()
    out_dir = os.path.dirname(os.path.abspath(outtree_file))
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    with tempfile.TemporaryDirectory() as tmp_dir:
        ran_ok = False
        last_err = ""

        if nj_tree_file and os.path.isfile(nj_tree_file) and os.path.getsize(nj_tree_file) > 0:
            clean_nj = _sanitize_newick_for_fastme(nj_tree_file, tmp_dir)
            cmd_user = [
                fastme_exe,
                "-i", matrix_file,
                "-u", clean_nj,
                "-o", outtree_file,
                "-n", nni_type,
                "-q",
            ]
            res = subprocess.run(cmd_user, capture_output=True, text=True)
            if res.returncode == 0 and os.path.isfile(outtree_file) and os.path.getsize(outtree_file) > 0:
                ran_ok = True
            else:
                last_err = res.stderr.strip()

        if not ran_ok:
            cmd_nj_nni = [
                fastme_exe,
                "-i", matrix_file,
                "-o", outtree_file,
                "-m", "N",
                "-n", nni_type,
                "-q",
            ]
            res = subprocess.run(cmd_nj_nni, capture_output=True, text=True)
            if res.returncode != 0 or not os.path.isfile(outtree_file) or os.path.getsize(outtree_file) == 0:
                raise RuntimeError(
                    f"FastME NNI execution failed (exit code {res.returncode}). "
                    f"Prior -u stderr: {last_err}; -m N stderr: {res.stderr.strip()}"
                )


def parse_args() -> argparse.Namespace:
    """
    Parse command-line arguments for NJ+NNI tree refinement.

    Returns
    -------
    argparse.Namespace
        Parsed CLI arguments.
    """
    parser = argparse.ArgumentParser(
        description="Refine an NJ tree with FastME NNI using a PHYLIP distance matrix."
    )
    parser.add_argument("--matrix", required=True, help="Input PHYLIP distance matrix file")
    parser.add_argument("--nj_tree", default=None, help="Optional starting NJ Newick tree file")
    parser.add_argument("--outtree", required=True, help="Output NNI-refined Newick tree file")
    parser.add_argument(
        "--nni_type",
        choices=["B", "O"],
        default="B",
        help="FastME NNI type: B (BalME, default) or O (OLS)",
    )
    return parser.parse_args()


def main() -> None:
    """
    CLI entry point for FastME NJ+NNI refinement.

    Returns
    -------
    None
    """
    args = parse_args()
    run_fastme_nni(
        matrix_file=args.matrix,
        outtree_file=args.outtree,
        nj_tree_file=args.nj_tree,
        nni_type=args.nni_type,
    )


if __name__ == "__main__":
    main()
