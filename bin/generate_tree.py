#!/usr/bin/env python3
"""
Generate Phylogenetic Tree CLI Script.

Generates synthetic Newick tree files according to selected tree simulation models,
primarily the Matsui & Iwasaki (2020) backward Yule process with logarithmic branch lengths.
"""

from typing import Optional
import os
import sys
import argparse

# Enable importing from src/
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

from src.tree_generator import generate_paper_yule_tree


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
