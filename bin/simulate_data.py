#!/usr/bin/env python3
"""
Simulate Data Script
Generates true phylogenetic trees (relaxed clock log-normal branch length variation, customizable sigma)
and protein sequence alignments with substitutions (LG+G with customizable alpha) and indels (customizable indel_rate).
"""

import sys
import os
import math
import random
import argparse
import subprocess
import dendropy
from Bio import SeqIO



def generate_random_tree(num_taxa, total_depth_substitutions, sigma=0.5):
    """
    Generates a random birth-death tree with N taxa, scaled to total depth D,
    and applies log-normal relaxed clock rate variation across branches (sigma).
    """
    tns = dendropy.TaxonNamespace([f"Taxon_{i+1}" for i in range(num_taxa)])
    tree = dendropy.simulate.treesim.birth_death_tree(birth_rate=1.0, death_rate=0.2, num_extant_tips=num_taxa, taxon_namespace=tns)

    # Scale max depth to D (substitutions/site)
    max_depth = max(nd.distance_from_root() for nd in tree.leaf_nodes())
    if max_depth > 0:
        scale_factor = total_depth_substitutions / max_depth
        for edge in tree.postorder_edge_iter():
            if edge.length is not None:
                edge.length *= scale_factor

    # Apply relaxed clock lognormal noise (sigma)
    for edge in tree.postorder_edge_iter():
        if edge.length is not None:
            multiplier = random.lognormvariate(0.0, sigma)
            edge.length = max(0.001, edge.length * multiplier)

    return tree

def format_model_string(model, alpha=None):
    """
    Formats the substitution model string for AliSim.
    If the model string already specifies parameters (contains '{' and '}'), returns as is.
    If the model contains '+G' (Gamma rate heterogeneity) and alpha is provided, formats as 'MODEL+G{alpha}'.
    If the model does not have '+G' (e.g., 'LG', 'WAG'), returns model string without alpha.
    """
    if not model:
        model = "LG+G"

    # Already has parameter specification (e.g. "LG+G{0.5}")
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

def simulate_with_alisim(tree_file, length, out_fasta, seed=None, indel_rate=0.05, model="LG+G", alpha=1.0, expected_num_taxa=16):
    """
    Simulates sequences using IQ-TREE AliSim with specified substitution model, Gamma alpha, and indel rate.
    Raises RuntimeError if AliSim is unavailable, fails, or produces incomplete/empty sequences.
    """
    iqtree_cmd = None
    for bin_name in ["iqtree2", "iqtree"]:
        try:
            res = subprocess.run([bin_name, "--version"], capture_output=True, text=True)
            if res.returncode == 0:
                iqtree_cmd = bin_name
                break
        except FileNotFoundError:
            pass

    if not iqtree_cmd:
        raise RuntimeError("Error: AliSim (iqtree2 / iqtree) is not found in PATH. Please ensure IQ-TREE 2 is installed and accessible.")

    prefix = out_fasta + ".alisim"
    model_str = format_model_string(model, alpha)
    cmd = [
        iqtree_cmd,
        "--alisim", prefix,
        "-m", model_str,
        "--length", str(length),
        "-t", tree_file,
        "--indel", f"{indel_rate},{indel_rate}",
        "--redo"
    ]
    if seed is not None:
        cmd.extend(["-seed", str(seed)])

    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        raise RuntimeError(f"Error: AliSim execution failed with returncode {res.returncode}.\nCommand: {' '.join(cmd)}\nSTDOUT:\n{res.stdout}\nSTDERR:\n{res.stderr}")

    unaligned_fa = prefix + ".unaligned.fa"
    if not (os.path.exists(unaligned_fa) and os.path.getsize(unaligned_fa) > 0):
        raise RuntimeError(f"Error: AliSim did not produce expected output FASTA: {unaligned_fa}\nCommand: {' '.join(cmd)}\nSTDOUT:\n{res.stdout}\nSTDERR:\n{res.stderr}")

    # Standardize taxon names (strip leading/trailing underscores added by AliSim)
    records = list(SeqIO.parse(unaligned_fa, "fasta"))
    
    # Check that ALL taxa exist and have non-empty sequences (at least 5 amino acids)
    min_required_len = min(5, max(1, length // 10))
    if len(records) < expected_num_taxa:
        raise RuntimeError(f"Error: AliSim generated only {len(records)} taxa, but expected {expected_num_taxa}.")
    
    for rec in records:
        clean_seq = str(rec.seq).replace("-", "")
        if len(clean_seq) < min_required_len:
            raise RuntimeError(f"Error: Taxon {rec.id} has too short/empty sequence ({len(clean_seq)} aa < {min_required_len}) in AliSim output.")

    with open(out_fasta, "w") as out_f:
        for rec in records:
            clean_id = rec.id.strip("_")
            out_f.write(f">{clean_id}\n{str(rec.seq)}\n")

    # Clean up temporary AliSim files
    for ext in [".phy", ".unaligned.fa"]:
        fpath = prefix + ext
        if os.path.exists(fpath):
            try:
                os.remove(fpath)
            except OSError:
                pass

def main():
    parser = argparse.ArgumentParser(description="Simulate True Tree and Protein Sequences with Indels")
    parser.add_argument("--distance", type=float, default=1.0, help="Evolutionary distance D (substitutions/site)")
    parser.add_argument("--length", type=int, default=300, help="Sequence length L (aa)")
    parser.add_argument("--num_taxa", type=int, default=16, help="Taxon count N")
    parser.add_argument("--sigma", type=float, default=0.5, help="Lognormal relaxed clock rate variation sigma (default: 0.5)")
    parser.add_argument("--alpha", type=float, default=1.0, help="Gamma shape parameter alpha for LG+G (default: 1.0)")
    parser.add_argument("--model", type=str, default="LG+G", help="Substitution model prefix (default: LG+G)")
    parser.add_argument("--indel_rate", type=float, default=0.05, help="Indel rate for insertions and deletions (default: 0.05)")
    parser.add_argument("--outtree", required=True, help="Output true Newick tree file")
    parser.add_argument("--outfasta", required=True, help="Output unaligned FASTA file")
    parser.add_argument("--seed", type=int, default=None, help="Random seed for reproducibility")
    args = parser.parse_args()

    base_seed = args.seed if args.seed is not None else random.randint(1, 1000000)
    max_trials = 10
    last_error = None

    for trial in range(max_trials):
        current_seed = base_seed + trial * 10007
        random.seed(current_seed)

        # 1. Generate True Tree
        tree = generate_random_tree(args.num_taxa, args.distance, sigma=args.sigma)
        tree.write(path=args.outtree, schema="newick")

        # 2. Simulate Sequences with AliSim
        try:
            simulate_with_alisim(
                tree_file=args.outtree,
                length=args.length,
                out_fasta=args.outfasta,
                seed=current_seed,
                indel_rate=args.indel_rate,
                model=args.model,
                alpha=args.alpha,
                expected_num_taxa=args.num_taxa
            )
            used_model_str = format_model_string(args.model, args.alpha)
            print(f"Simulated data generated (trial {trial+1}/{max_trials}): Tree -> {args.outtree}, FASTA -> {args.outfasta} (N={args.num_taxa}, D={args.distance}, L={args.length}, sigma={args.sigma}, model={used_model_str}, indel_rate={args.indel_rate})")
            return
        except Exception as e:
            last_error = e
            # Clean up partial output if any
            if os.path.exists(args.outfasta):
                try:
                    os.remove(args.outfasta)
                except OSError:
                    pass

    # If all trials failed, raise the last exception
    raise RuntimeError(f"Failed to generate valid sequence simulation after {max_trials} attempts. Last error: {last_error}")

if __name__ == "__main__":
    main()
