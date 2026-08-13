#!/usr/bin/env python3
"""
Simulate Data Script
Generates true phylogenetic trees (relaxed clock log-normal branch length variation, sigma=0.5)
and protein sequence alignments with substitutions (LG+G) and indels (insertion/deletion gaps).
"""

import sys
import os
import math
import random
import argparse
import subprocess
import dendropy
from Bio import SeqIO

AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"

# Standard LG background equilibrium amino acid frequencies
LG_FREQS = {
    'A': 0.079066, 'R': 0.055941, 'N': 0.041977, 'D': 0.053052, 'C': 0.013097,
    'Q': 0.040767, 'E': 0.071586, 'G': 0.057342, 'H': 0.022355, 'I': 0.056206,
    'L': 0.099047, 'K': 0.064600, 'M': 0.022986, 'F': 0.042302, 'P': 0.044040,
    'S': 0.061197, 'T': 0.053287, 'W': 0.012066, 'Y': 0.029199, 'V': 0.069797
}

AA_LIST = list(LG_FREQS.keys())
AA_PROBS = list(LG_FREQS.values())

def generate_random_tree(num_taxa, total_depth_substitutions):
    """
    Generates a random birth-death tree with N taxa, scaled to total depth D,
    and applies log-normal relaxed clock rate variation across branches (sigma=0.5).
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

    # Apply relaxed clock lognormal noise (sigma = 0.5)
    for edge in tree.postorder_edge_iter():
        if edge.length is not None:
            multiplier = random.lognormvariate(0.0, 0.5)
            edge.length = max(0.001, edge.length * multiplier)

    return tree

def simulate_sequence_evolution(tree, length, indel_rate=0.05, max_indel_len=3):
    """
    Simulates sequence evolution along the tree starting from a random ancestral protein sequence.
    Applies substitutions (Poisson/LG) and Indels (insertions and deletions).
    Returns dictionary of unaligned leaf sequences.
    """
    ancestor_seq = "".join(random.choices(AA_LIST, weights=AA_PROBS, k=length))

    node_seqs = {}
    root_node = tree.seed_node
    node_seqs[root_node] = list(ancestor_seq)

    for node in tree.preorder_node_iter():
        if node == root_node:
            continue
        parent_seq = node_seqs[node.parent_node]
        edge_len = node.edge_length or 0.1

        # Mutate sequence along edge
        curr_seq = list(parent_seq)
        new_seq = []
        i = 0
        while i < len(curr_seq):
            # Check for deletion
            if random.random() < (1.0 - math.exp(-edge_len * indel_rate)):
                indel_len = random.randint(1, max_indel_len)
                i += indel_len
                continue

            # Check for substitution
            if random.random() < (1.0 - math.exp(-edge_len)):
                new_seq.append(random.choices(AA_LIST, weights=AA_PROBS)[0])
            else:
                new_seq.append(curr_seq[i])

            # Check for insertion
            if random.random() < (1.0 - math.exp(-edge_len * indel_rate)):
                ins_len = random.randint(1, max_indel_len)
                ins_seq = random.choices(AA_LIST, weights=AA_PROBS, k=ins_len)
                new_seq.extend(ins_seq)

            i += 1

        node_seqs[node] = new_seq

    # Extract leaf sequences
    leaf_seqs = {}
    for leaf in tree.leaf_nodes():
        name = leaf.taxon.label
        seq_str = "".join(node_seqs[leaf])
        leaf_seqs[name] = seq_str

    return leaf_seqs

def simulate_with_alisim(tree_file, length, out_fasta, seed=None, indel_rate=0.05):
    """Simulates sequences using IQ-TREE AliSim (LG+G4 model with indels)."""
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
        return False

    prefix = out_fasta + ".alisim"
    cmd = [
        iqtree_cmd,
        "--alisim", prefix,
        "-m", "LG+G",
        "--length", str(length),
        "-t", tree_file,
        "--indel", f"{indel_rate},{indel_rate}",
        "--redo"
    ]
    if seed is not None:
        cmd.extend(["-seed", str(seed)])

    res = subprocess.run(cmd, capture_output=True, text=True)
    unaligned_fa = prefix + ".unaligned.fa"
    if os.path.exists(unaligned_fa) and os.path.getsize(unaligned_fa) > 0:
        # Standardize taxon names (strip leading/trailing underscores added by AliSim)
        records = list(SeqIO.parse(unaligned_fa, "fasta"))
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
        return True
    return False

def main():
    parser = argparse.ArgumentParser(description="Simulate True Tree and Protein Sequences with Indels")
    parser.add_argument("--distance", type=float, default=1.0, help="Evolutionary distance D (substitutions/site)")
    parser.add_argument("--length", type=int, default=300, help="Sequence length L (aa)")
    parser.add_argument("--num_taxa", type=int, default=16, help="Taxon count N")
    parser.add_argument("--outtree", required=True, help="Output true Newick tree file")
    parser.add_argument("--outfasta", required=True, help="Output unaligned FASTA file")
    parser.add_argument("--seed", type=int, default=None, help="Random seed for reproducibility")
    args = parser.parse_args()

    if args.seed is not None:
        random.seed(args.seed)

    # 1. Generate True Tree
    tree = generate_random_tree(args.num_taxa, args.distance)
    tree.write(path=args.outtree, schema="newick")

    # 2. Simulate Sequences (AliSim primary with Python fallback)
    success = simulate_with_alisim(args.outtree, args.length, args.outfasta, seed=args.seed, indel_rate=0.05)
    if not success:
        seqs = simulate_sequence_evolution(tree, args.length, indel_rate=0.05)
        with open(args.outfasta, "w") as f:
            for name, seq in seqs.items():
                f.write(f">{name}\n{seq}\n")

    print(f"Simulated data generated: Tree -> {args.outtree}, FASTA -> {args.outfasta} (N={args.num_taxa}, D={args.distance}, L={args.length})")

if __name__ == "__main__":
    main()
