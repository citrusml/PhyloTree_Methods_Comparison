"""
Module for Invariant Category Sites (ICS Full) sequence simulation and alignment.

Provides functions to generate synthetic phylogenetic trees under the backward Yule model,
simulate amino acid sequences under the 100% ICS+G4 model with Zipfian power-law indels
using AliSim, and generate estimated MSAs using MAFFT.
"""

from typing import Dict, Any, Tuple, Optional
import os
import sys
import subprocess
from pathlib import Path


def generate_paper_tree(
    taxa: int,
    scale: float,
    seed: int,
    outtree_path: str,
    rate_sd: float = 0.0,
    lba_ratio: float = 1.0,
    py_executable: str = sys.executable
) -> str:
    """
    Generate a synthetic phylogenetic tree under the backward Yule process (Matsui & Iwasaki 2020).

    Parameters
    ----------
    taxa : int
        Number of taxa (leaves) in the generated tree [unitless].
    scale : float
        Evolutionary distance scaling factor D [substitutions/site].
    seed : int
        Pseudorandom number generator seed [unitless].
    outtree_path : str
        Target file path for saving the Newick tree.
    rate_sd : float, optional
        Standard deviation for log-normal rate heterogeneity [unitless], default is 0.0.
    lba_ratio : float, optional
        Long-branch attraction ratio b/a [unitless], default is 1.0.
    py_executable : str, optional
        Python interpreter path to execute generate_tree.py.

    Returns
    -------
    str
        Path to the generated Newick tree file.

    Raises
    ------
    RuntimeError
        If tree generation fails.
    """
    cmd = [
        py_executable, "bin/generate_tree.py",
        "--taxa", str(taxa),
        "--scale", str(scale),
        "--seed", str(seed),
        "--model", "paper_yule",
        "--rate_sd", str(rate_sd),
        "--lba_ratio", str(lba_ratio),
        "--outtree", outtree_path
    ]
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        raise RuntimeError(f"Tree generation failed (code {res.returncode}): {res.stderr}")
    if not os.path.exists(outtree_path) or os.path.getsize(outtree_path) == 0:
        raise RuntimeError(f"Generated tree file is empty or missing: {outtree_path}")
    return outtree_path


def simulate_ics_full_sequences(
    tree_path: str,
    length: int,
    out_prefix: str,
    seed: int,
    model_nex_path: str = "models/ics_model.nex",
    insert_rate: float = 0.10,
    delete_rate: float = 0.10,
    indel_size: str = "POW{1.7/50},POW{1.7/50}",
    iqtree_bin: str = "bin/iqtree2"
) -> Tuple[str, str]:
    """
    Simulate protein sequences under 100% Invariant Category Sites (ICS+G4{1.0}) with Zipfian indels.

    Parameters
    ----------
    tree_path : str
        Path to the reference Newick tree file.
    length : int
        Root ancestral sequence length L [amino acid residues / sites].
    out_prefix : str
        Prefix for AliSim output files.
    seed : int
        Pseudorandom number generator seed [unitless].
    model_nex_path : str, optional
        Path to custom NEXUS model definition file containing ICS model definition.
    insert_rate : float, optional
        Insertion rate relative to substitution rate [unitless], default is 0.10.
    delete_rate : float, optional
        Deletion rate relative to substitution rate [unitless], default is 0.10.
    indel_size : str, optional
        Zipfian power-law indel length distribution spec, default is 'POW{1.7/50},POW{1.7/50}'.
    iqtree_bin : str, optional
        Path to iqtree / alisim binary.

    Returns
    -------
    Tuple[str, str]
        (unaligned_fasta_path, true_msa_fasta_path)

    Raises
    ------
    RuntimeError
        If sequence simulation fails or outputs are invalid.
    """
    cmd = [
        iqtree_bin, "--alisim", out_prefix,
        "--mdef", model_nex_path, "--seqtype", "AA",
        "-m", "ICS+G4{1.0}",
        "--length", str(length),
        "-t", tree_path,
        "--indel", f"{insert_rate},{delete_rate}",
        "--indel-size", indel_size,
        "-af", "fasta",
        "-seed", str(seed),
        "--redo"
    ]
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        raise RuntimeError(f"AliSim simulation failed (code {res.returncode}): {res.stderr}")

    unaligned_path = f"{out_prefix}.unaligned.fa"
    true_msa_path = f"{out_prefix}.fa"

    if not os.path.exists(unaligned_path) or not os.path.exists(true_msa_path):
        raise RuntimeError(f"Expected AliSim outputs not found: {unaligned_path}, {true_msa_path}")

    return unaligned_path, true_msa_path


def align_with_mafft(
    unaligned_fasta: str,
    out_msa_path: str,
    mafft_bin: str = "mafft",
    threads: int = 2
) -> str:
    """
    Align unaligned sequences using MAFFT (--auto --quiet).

    Parameters
    ----------
    unaligned_fasta : str
        Path to input unaligned FASTA file.
    out_msa_path : str
        Path to output estimated MSA FASTA file.
    mafft_bin : str, optional
        Path or command name for mafft binary.
    threads : int, optional
        Number of parallel CPU threads, default is 2.

    Returns
    -------
    str
        Path to output estimated MSA FASTA file.

    Raises
    ------
    RuntimeError
        If MAFFT execution fails.
    """
    cmd = f"{mafft_bin} --threadit 0 --auto --thread {threads} --quiet {unaligned_fasta} > {out_msa_path}"
    res = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if res.returncode != 0:
        raise RuntimeError(f"MAFFT failed (code {res.returncode}): {res.stderr}")
    if not os.path.exists(out_msa_path) or os.path.getsize(out_msa_path) == 0:
        raise RuntimeError(f"Estimated MSA file is empty or missing: {out_msa_path}")
    return out_msa_path
