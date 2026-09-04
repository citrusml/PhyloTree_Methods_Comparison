#!/usr/bin/env python3
"""
Generate Phylogenetic Tree CLI Script.

Generates synthetic Newick tree files according to selected tree simulation models,
primarily the Matsui & Iwasaki (2020) backward Yule process with logarithmic branch lengths.
"""

from typing import List, Optional
import os
import sys
import math
import random
import argparse

# Enable importing from src/
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

try:
    from src.tree_generator import generate_paper_yule_tree
except ImportError:
    # Embedded fallback to ensure execution even if src/ is missing on remote nodes
    class TreeNode:
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
            return self.left is None and self.right is None

        def to_newick(self) -> str:
            if self.is_leaf():
                return f"{self.name}:{self.length:.6f}"
            assert self.left is not None and self.right is not None
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
        u = rng.random()
        base_l = 1.0 - math.log(u * (math.e - 1.0) + 1.0)
        if rate_sd > 0.0:
            base_l = max(0.0, base_l + rng.gauss(0.0, rate_sd))
        return max(min_length, base_l * scale)

    def generate_paper_yule_tree(
        taxa: int = 32,
        scale: float = 1.0,
        seed: Optional[int] = None,
        rate_sd: float = 0.0,
        lba_ratio: float = 1.0,
    ) -> str:
        if taxa < 2:
            raise ValueError(f"Taxa count must be at least 2, got {taxa}")
        if scale <= 0.0:
            raise ValueError(f"Scale must be strictly positive, got {scale}")

        rng = random.Random(seed)
        nodes: List[TreeNode] = [TreeNode(name=f"T{i + 1}") for i in range(taxa)]

        while len(nodes) > 1:
            idx1, idx2 = sorted(rng.sample(range(len(nodes)), 2), reverse=True)
            v_node = nodes.pop(idx1)
            u_node = nodes.pop(idx2)

            lu = sample_paper_branch_length(rng, scale=scale, rate_sd=rate_sd)
            lv = sample_paper_branch_length(rng, scale=scale, rate_sd=rate_sd)

            if lba_ratio > 1.0:
                if u_node.name in ("T1", f"T{taxa}"):
                    lu *= lba_ratio
                if v_node.name in ("T1", f"T{taxa}"):
                    lv *= lba_ratio

            u_node.length = lu
            v_node.length = lv
            parent = TreeNode(left=u_node, right=v_node)
            nodes.append(parent)

        return nodes[0].to_newick() + ";"


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
