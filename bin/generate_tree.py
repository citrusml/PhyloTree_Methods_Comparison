#!/usr/bin/env python3
"""
Generate Phylogenetic Tree CLI Script: generate_tree.py

Self-contained script to generate synthetic Newick tree files according to
Matsui & Iwasaki (2020, Systematic Biology) simulation specifications:
1. Topology: Backward Yule process (Kingman's coalescent pairing).
2. Branch Lengths: Sampled from logarithmic distribution l = 1 - ln(u * (e - 1) + 1) and scaled by D.
3. Taxa: Terminal leaves named T1, T2, ..., TN for AliSim compatibility.

Completely standalone with zero external dependencies (uses standard library only).
"""

from typing import List, Optional
import os
import sys
import math
import random
import argparse


class TreeNode:
    """
    Binary tree node representation for phylogenetic trees.

    Parameters
    ----------
    name : str, default=""
        Taxon name for leaf nodes (e.g., 'T1', 'T2'), or empty string for internal nodes.
    length : float, default=0.0
        Branch length leading to the parent node (substitutions per site).
    left : TreeNode or None, default=None
        Left child node.
    right : TreeNode or None, default=None
        Right child node.
    """

    def __init__(
        self,
        name: str = "",
        length: float = 0.0,
        left: Optional["TreeNode"] = None,
        right: Optional["TreeNode"] = None,
    ) -> None:
        self.name: str = name
        self.length: float = length
        self.left: Optional["TreeNode"] = left
        self.right: Optional["TreeNode"] = right

    def is_leaf(self) -> bool:
        """
        Checks whether this node is a leaf (terminal node).

        Returns
        -------
        bool
            True if node has no children, False otherwise.
        """
        return self.left is None and self.right is None

    def to_newick(self) -> str:
        """
        Recursively converts the subtree rooted at this node to Newick format string.

        Returns
        -------
        str
            Newick formatted representation of the subtree.
        """
        if self.is_leaf():
            return f"{self.name}:{self.length:.6f}"

        assert self.left is not None and self.right is not None, "Binary tree node must have both children"
        children_str = f"({self.left.to_newick()},{self.right.to_newick()})"
        if self.length > 0.0:
            return f"{children_str}:{self.length:.6f}"
        return children_str


def sample_paper_branch_length(
    rng: random.Random,
    scale: float = 1.0,
    rate_sd: float = 0.0,
    min_length: float = 1e-6,
) -> float:
    """
    Samples a single branch length according to the Matsui & Iwasaki (2020) model.

    The base branch length follows the distribution:
        l = 1 - ln(u * (e - 1) + 1)
    where u ~ Uniform(0, 1), and e is Euler's number (math.e).

    Parameters
    ----------
    rng : random.Random
        Random number generator instance for reproducible sampling.
    scale : float, default=1.0
        Global evolutionary distance scaling factor D (substitutions per site scaling).
    rate_sd : float, default=0.0
        Standard deviation of Gaussian rate variation across branches (dimensionless).
    min_length : float, default=1e-6
        Minimum allowable branch length (substitutions per site).

    Returns
    -------
    float
        Sampled branch length (substitutions per site).
    """
    u = rng.random()
    # Matsui & Iwasaki (2020) formula: l = 1 - ln(u * (e - 1) + 1)
    base_l = 1.0 - math.log(u * (math.e - 1.0) + 1.0)

    # Optional rate heterogeneity (Gaussian perturbation)
    if rate_sd > 0.0:
        perturbation = rng.gauss(0.0, rate_sd)
        base_l = max(0.0, base_l + perturbation)

    branch_length = max(min_length, base_l * scale)
    return branch_length


def generate_paper_yule_tree(
    taxa: int = 32,
    scale: float = 1.0,
    seed: Optional[int] = None,
    rate_sd: float = 0.0,
    lba_ratio: float = 1.0,
) -> str:
    """
    Generates a phylogenetic tree based on the backward Yule process and paper branch lengths.

    Parameters
    ----------
    taxa : int, default=32
        Number of terminal taxa (leaves) in the generated tree. Must be >= 2.
    scale : float, default=1.0
        Evolutionary distance scaling factor D (substitutions per site scaling).
    seed : int or None, default=None
        Random seed for deterministic, reproducible tree generation.
    rate_sd : float, default=0.0
        Rate heterogeneity standard deviation across branches (dimensionless).
    lba_ratio : float, default=1.0
        Long-branch attraction ratio (dimensionless, b/a >= 1.0).

    Returns
    -------
    str
        Valid Newick tree string terminated with a semicolon ';'.

    Raises
    ------
    ValueError
        If taxa < 2 or scale <= 0.
    """
    if taxa < 2:
        raise ValueError(f"Taxa count must be at least 2, got {taxa}")
    if scale <= 0.0:
        raise ValueError(f"Scale must be strictly positive, got {scale}")

    rng = random.Random(seed)

    # Initialize leaves: T1, T2, ..., TN
    nodes: List[TreeNode] = [TreeNode(name=f"T{i + 1}") for i in range(taxa)]

    # Backward Yule process: randomly pair nodes until only root remains
    while len(nodes) > 1:
        # Pick 2 distinct nodes uniformly at random without replacement
        idx1, idx2 = sorted(rng.sample(range(len(nodes)), 2), reverse=True)
        # Pop in descending order so indices remain valid
        v_node = nodes.pop(idx1)
        u_node = nodes.pop(idx2)

        # Sample independent branch lengths
        lu = sample_paper_branch_length(rng, scale=scale, rate_sd=rate_sd)
        lv = sample_paper_branch_length(rng, scale=scale, rate_sd=rate_sd)

        # Apply LBA scaling if requested for specific taxa (T1 and TN)
        if lba_ratio > 1.0:
            if u_node.name in ("T1", f"T{taxa}"):
                lu *= lba_ratio
            if v_node.name in ("T1", f"T{taxa}"):
                lv *= lba_ratio

        u_node.length = lu
        v_node.length = lv

        # Create parent internal node
        parent = TreeNode(left=u_node, right=v_node)
        nodes.append(parent)

    root = nodes[0]
    return root.to_newick() + ";"


def main() -> None:
    """
    Parses CLI arguments, generates the specified tree, and writes to output file.
    """
    parser = argparse.ArgumentParser(
        description="Generate synthetic phylogenetic trees (Matsui & Iwasaki 2020 paper model)"
    )
    parser.add_argument(
        "--taxa",
        type=int,
        default=32,
        help="Number of terminal taxa (default: 32)",
    )
    parser.add_argument(
        "--scale",
        type=float,
        default=1.0,
        help="Evolutionary distance scaling factor D (default: 1.0)",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Random seed for reproducibility (default: None)",
    )
    parser.add_argument(
        "--model",
        choices=["paper_yule", "yule"],
        default="paper_yule",
        help="Tree generation model (default: paper_yule)",
    )
    parser.add_argument(
        "--rate_sd",
        type=float,
        default=0.0,
        help="Rate heterogeneity standard deviation (default: 0.0)",
    )
    parser.add_argument(
        "--lba_ratio",
        type=float,
        default=1.0,
        help="Long-branch attraction ratio b/a (default: 1.0)",
    )
    parser.add_argument(
        "--outtree",
        type=str,
        required=True,
        help="Output Newick tree filepath",
    )

    args = parser.parse_args()

    # Generate tree
    if args.model in ("paper_yule", "yule"):
        tree_newick = generate_paper_yule_tree(
            taxa=args.taxa,
            scale=args.scale,
            seed=args.seed,
            rate_sd=args.rate_sd,
            lba_ratio=args.lba_ratio,
        )
    else:
        raise ValueError(f"Unsupported model: {args.model}")

    # Ensure parent directory exists
    out_dir = os.path.dirname(os.path.abspath(args.outtree))
    if out_dir and not os.path.exists(out_dir):
        os.makedirs(out_dir, exist_ok=True)

    with open(args.outtree, "w") as f:
        f.write(tree_newick + "\n")

    print(f"Generated {args.model} tree ({args.taxa} taxa, scale={args.scale}, seed={args.seed}) -> {args.outtree}")


if __name__ == "__main__":
    main()
