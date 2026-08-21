#!/usr/bin/env python3
"""
Simulate Data Script
Generates true phylogenetic trees using IQ-TREE AliSim's RANDOM Birth-Death model
and simulates protein sequence alignments with substitutions (LG+G4 with customizable alpha)
and indels (customizable insert_rate and delete_rate).
"""

import sys
import os
import math
import random
import argparse
import subprocess
import shutil
from Bio import SeqIO

def format_model_string(model, alpha=None):
    """
    Formats the substitution model string for AliSim.
    If the model string already specifies parameters (contains '{' and '}'), returns as is.
    If the model contains '+G' or '+G4' (Gamma rate heterogeneity) and alpha is provided, formats as 'MODEL+G4{alpha}'.
    If the model does not have '+G' (e.g., 'LG', 'WAG'), returns model string without alpha.
    """
    if not model:
        model = "LG+G4"

    # Already has parameter specification (e.g. "LG+G4{0.5}")
    if "{" in model and "}" in model:
        return model

    # If Gamma rate variation is in model, attach alpha if provided
    if "+G" in model:
        if alpha is not None:
            import re
            return re.sub(r'(\+G\d*)(\{[^}]*\})?', rf'\1{{{alpha}}}', model)
        return model
    else:
        # Non-gamma model (e.g., "LG", "WAG", "JC")
        return model

def find_iqtree_cmd():
    """Locates iqtree executable (iqtree3, iqtree2, or iqtree)."""
    for bin_name in ["iqtree3", "iqtree2", "iqtree"]:
        try:
            res = subprocess.run([bin_name, "--version"], capture_output=True, text=True)
            if res.returncode == 0:
                return bin_name
        except FileNotFoundError:
            pass
    return None

def simulate_tree_and_sequences_with_alisim(
    out_tree,
    out_fasta,
    out_true_msa=None,
    num_taxa=32,
    distance=1.0,
    length=300,
    birth_rate=0.1,
    death_rate=0.05,
    insert_rate=0.05,
    delete_rate=0.10,
    model="LG+G4",
    alpha=1.0,
    seed=None
):
    """
    Simulates a random Birth-Death phylogenetic tree, sequences, and true MSA using IQ-TREE AliSim.
    Tree is generated via -t "RANDOM{bd{birth_rate/death_rate}/num_taxa}" with exponential branch lengths
    and scaled by branch scale (distance).
    """
    iqtree_cmd = find_iqtree_cmd()
    if not iqtree_cmd:
        raise RuntimeError("Error: AliSim (iqtree3 / iqtree2 / iqtree) is not found in PATH. Please ensure IQ-TREE is installed and accessible.")

    prefix = out_fasta + ".alisim"
    model_str = format_model_string(model, alpha)
    tree_random_spec = f"RANDOM{{bd{{{birth_rate}/{death_rate}}}/{num_taxa}}}"

    cmd = [
        iqtree_cmd,
        "--alisim", prefix,
        "-m", model_str,
        "--length", str(length),
        "-t", tree_random_spec,
        "--indel", f"{insert_rate},{delete_rate}",
        "--branch-scale", str(distance),
        "--redo"
    ]
    if seed is not None:
        cmd.extend(["-seed", str(seed)])

    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        raise RuntimeError(f"Error: AliSim execution failed with returncode {res.returncode}.\nCommand: {' '.join(cmd)}\nSTDOUT:\n{res.stdout}\nSTDERR:\n{res.stderr}")

    # 1. Process True Tree file
    alisim_treefile = prefix + ".treefile"
    if not (os.path.exists(alisim_treefile) and os.path.getsize(alisim_treefile) > 0):
        alt_treefile = prefix + ".tree"
        if os.path.exists(alt_treefile) and os.path.getsize(alt_treefile) > 0:
            alisim_treefile = alt_treefile
        else:
            raise RuntimeError(f"Error: AliSim did not produce expected output treefile: {alisim_treefile}\nCommand: {' '.join(cmd)}\nSTDOUT:\n{res.stdout}\nSTDERR:\n{res.stderr}")

    with open(alisim_treefile, "r") as f_in, open(out_tree, "w") as f_out:
        tree_content = f_in.read()
        f_out.write(tree_content)

    # 2. Process Sequences and True Alignment files
    # AliSim outputs .fa or .phy for True MSA, and .unaligned.fa for unaligned sequences
    alisim_aligned_fa = prefix + ".fa"
    alisim_unaligned_fa = prefix + ".unaligned.fa"
    alisim_phy = prefix + ".phy"

    aligned_records = []
    if os.path.exists(alisim_aligned_fa) and os.path.getsize(alisim_aligned_fa) > 0:
        aligned_records = list(SeqIO.parse(alisim_aligned_fa, "fasta"))
    elif os.path.exists(alisim_phy) and os.path.getsize(alisim_phy) > 0:
        aligned_records = list(SeqIO.parse(alisim_phy, "phylip"))

    if os.path.exists(alisim_unaligned_fa) and os.path.getsize(alisim_unaligned_fa) > 0:
        unaligned_records = list(SeqIO.parse(alisim_unaligned_fa, "fasta"))
    elif aligned_records:
        unaligned_records = aligned_records
    else:
        raise RuntimeError(f"Error: AliSim did not produce expected sequence files.\nCommand: {' '.join(cmd)}\nSTDOUT:\n{res.stdout}\nSTDERR:\n{res.stderr}")

    # Check taxa count
    min_required_len = 1
    if len(unaligned_records) < num_taxa:
        raise RuntimeError(f"Error: AliSim generated only {len(unaligned_records)} taxa, but expected {num_taxa}.")

    # Write unaligned FASTA (gaps stripped)
    with open(out_fasta, "w") as out_f:
        for rec in unaligned_records:
            clean_id = rec.id.strip("_")
            clean_seq = str(rec.seq).replace("-", "")
            if len(clean_seq) < min_required_len:
                raise RuntimeError(f"Error: Taxon {rec.id} has empty sequence in AliSim output.")
            out_f.write(f">{clean_id}\n{clean_seq}\n")

    # Write True MSA FASTA if requested
    if out_true_msa:
        with open(out_true_msa, "w") as out_m:
            if aligned_records:
                for rec in aligned_records:
                    clean_id = rec.id.strip("_")
                    out_m.write(f">{clean_id}\n{str(rec.seq)}\n")
            else:
                # If aligned FASTA wasn't present, write unaligned sequences as fallback
                for rec in unaligned_records:
                    clean_id = rec.id.strip("_")
                    out_m.write(f">{clean_id}\n{str(rec.seq)}\n")

    # Clean up temporary AliSim generated files
    for ext in [".phy", ".unaligned.fa", ".fa", ".treefile", ".tree", ".log", ".iqtree"]:
        fpath = prefix + ext
        if os.path.exists(fpath):
            try:
                os.remove(fpath)
            except OSError:
                pass

def main():
    parser = argparse.ArgumentParser(description="Simulate True Tree and Protein Sequences with Indels using AliSim")
    parser.add_argument("--distance", type=float, default=1.0, help="Branch length scale / Evolutionary distance D (default: 1.0)")
    parser.add_argument("--length", type=int, default=300, help="Sequence length L (aa) (default: 300)")
    parser.add_argument("--num_taxa", type=int, default=32, help="Taxon count N (default: 32)")
    parser.add_argument("--birth_rate", type=float, default=0.1, help="Birth rate for Birth-Death tree model (default: 0.1)")
    parser.add_argument("--death_rate", type=float, default=0.05, help="Death rate for Birth-Death tree model (default: 0.05)")
    parser.add_argument("--insert_rate", type=float, default=0.05, help="Insertion rate relative to substitutions (default: 0.05)")
    parser.add_argument("--delete_rate", type=float, default=0.10, help="Deletion rate relative to substitutions (default: 0.10)")
    parser.add_argument("--indel_rate", type=float, default=None, help="Symmetric indel rate for both insertion and deletion (overrides insert/delete if specified)")
    parser.add_argument("--alpha", type=float, default=1.0, help="Gamma shape parameter alpha for LG+G4 (default: 1.0)")
    parser.add_argument("--model", type=str, default="LG+G4", help="Substitution model prefix (default: LG+G4)")
    parser.add_argument("--outtree", required=True, help="Output true Newick tree file")
    parser.add_argument("--outfasta", required=True, help="Output unaligned FASTA file")
    parser.add_argument("--outtrue_msa", help="Output True MSA (aligned) FASTA file")
    parser.add_argument("--seed", type=int, default=None, help="Random seed for reproducibility")
    args = parser.parse_args()

    insert_rate = args.insert_rate if args.indel_rate is None else args.indel_rate
    delete_rate = args.delete_rate if args.indel_rate is None else args.indel_rate

    base_seed = args.seed if args.seed is not None else random.randint(1, 1000000)
    max_trials = 10
    last_error = None

    for trial in range(max_trials):
        current_seed = base_seed + trial * 10007
        random.seed(current_seed)

        try:
            simulate_tree_and_sequences_with_alisim(
                out_tree=args.outtree,
                out_fasta=args.outfasta,
                out_true_msa=args.outtrue_msa,
                num_taxa=args.num_taxa,
                distance=args.distance,
                length=args.length,
                birth_rate=args.birth_rate,
                death_rate=args.death_rate,
                insert_rate=insert_rate,
                delete_rate=delete_rate,
                model=args.model,
                alpha=args.alpha,
                seed=current_seed
            )
            used_model_str = format_model_string(args.model, args.alpha)
            print(f"Simulated data generated (trial {trial+1}/{max_trials}): Tree -> {args.outtree}, FASTA -> {args.outfasta} (N={args.num_taxa}, D={args.distance}, L={args.length}, bd={{{args.birth_rate}/{args.death_rate}}}, model={used_model_str}, indel={{{insert_rate},{delete_rate}}})")
            return
        except Exception as e:
            last_error = e
            # Clean up partial output if any
            for p in [args.outfasta, args.outtree, args.outtrue_msa]:
                if p and os.path.exists(p):
                    try:
                        os.remove(p)
                    except OSError:
                        pass

    # If all trials failed, raise the last exception
    raise RuntimeError(f"Failed to generate valid sequence simulation after {max_trials} attempts. Last error: {last_error}")

if __name__ == "__main__":
    main()
