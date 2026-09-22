#!/usr/bin/env python3
"""
Generate Phylogenetic Tree CLI Script: generate_tree.py

Generates synthetic Newick tree files according to Matsui & Iwasaki (2020, Systematic Biology)
simulation specifications:
1. Primary engine: BioPerl 1.6.924 (Bio::Tree::RandomFactory, rand_yule_c_tree) via bin/generate_tree.pl
2. Topology: Backward Yule process with constant birth rate lambda = 1.0.
3. Branch Lengths: Sampled from logarithmic distribution l = 1 - ln(u * (e - 1) + 1) and scaled by D.
4. Taxa: Terminal leaves named T1, T2, ..., TN for AliSim and INDELible compatibility.
5. Robust fallback: Pure Python implementation matching BioPerl specifications if Perl/BioPerl is absent.
"""

from typing import List, Optional
import os
import sys
import math
import random
import shutil
import argparse
import tempfile
import subprocess


class TreeNode:
    """
    Binary tree node representation for phylogenetic trees (Python fallback engine).

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
        Recursively converts the subtree rooted at this node to Newick format string without internal node labels.

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
    base_l = 1.0 - math.log(u * (math.e - 1.0) + 1.0)

    if rate_sd > 0.0:
        perturbation = rng.gauss(0.0, rate_sd)
        base_l = max(0.0, base_l + perturbation)

    branch_length = max(min_length, base_l * scale)
    return branch_length


def _generate_paper_yule_tree_python(
    taxa: int = 32,
    scale: float = 1.0,
    seed: Optional[int] = None,
    rate_sd: float = 0.0,
    lba_ratio: float = 1.0,
    min_length: float = 1e-6,
) -> str:
    """
    Generates a phylogenetic tree using pure Python, strictly implementing
    BioPerl's Bio::Tree::RandomFactory rand_yule_c_tree coalescent pairing algorithm.
    """
    rng = random.Random(seed)

    # 1. Sample N - 1 coalescence times
    times: List[float] = []
    for _ in range(taxa - 1):
        u = rng.random()
        base_t = 1.0 - math.log(u * (math.e - 1.0) + 1.0)
        times.append(base_t * scale)
    times.sort()

    # 2. Initialize leaves: T1, T2, ..., TN
    nodes: List[TreeNode] = [TreeNode(name=f"T{i + 1}") for i in range(taxa)]

    # 3. Backward Yule coalescence
    while len(nodes) > 1:
        t = times.pop(0)

        idx1, idx2 = sorted(rng.sample(range(len(nodes)), 2), reverse=True)
        v_node = nodes.pop(idx1)
        u_node = nodes.pop(idx2)

        lu = t
        lv = t

        if rate_sd > 0.0:
            lu = max(min_length, lu + rng.gauss(0.0, rate_sd * scale))
            lv = max(min_length, lv + rng.gauss(0.0, rate_sd * scale))

        if lba_ratio > 1.0:
            if u_node.name in ("T1", f"T{taxa}"):
                lu *= lba_ratio
            if v_node.name in ("T1", f"T{taxa}"):
                lv *= lba_ratio

        u_node.length = max(min_length, lu)
        v_node.length = max(min_length, lv)

        parent = TreeNode(left=u_node, right=v_node)
        nodes.append(parent)

    root = nodes[0]
    return root.to_newick() + ";"


def get_perl_executable() -> Optional[str]:
    """
    Finds a perl executable capable of running BioPerl (Bio::Tree::RandomFactory).
    Checks conda/micromamba environment first, then PATH and standard locations.
    """
    candidates = []

    # 1. Same directory as current python interpreter (conda/micromamba env)
    env_perl = os.path.join(os.path.dirname(sys.executable), "perl")
    if os.path.isfile(env_perl) and os.access(env_perl, os.X_OK):
        candidates.append(env_perl)

    # 2. CONDA_PREFIX if set
    conda_prefix = os.environ.get("CONDA_PREFIX")
    if conda_prefix:
        cp_perl = os.path.join(conda_prefix, "bin", "perl")
        if os.path.isfile(cp_perl) and os.access(cp_perl, os.X_OK):
            candidates.append(cp_perl)

    # 3. perl from PATH
    which_perl = shutil.which("perl")
    if which_perl and which_perl not in candidates:
        candidates.append(which_perl)

    # Test each candidate for Bio::Tree::RandomFactory
    for p in candidates:
        try:
            res = subprocess.run(
                [p, "-MBio::Tree::RandomFactory", "-e", "1"],
                capture_output=True,
                text=True,
                timeout=5,
            )
            if res.returncode == 0:
                return p
        except Exception:
            continue

    # Fallback to first available perl candidate if any
    return candidates[0] if candidates else (which_perl or "perl")


def is_bioperl_available() -> bool:
    """
    Checks if BioPerl (Bio::Tree::RandomFactory) is available in any perl executable.
    """
    perl_bin = get_perl_executable()
    if not perl_bin:
        return False
    try:
        res = subprocess.run(
            [perl_bin, "-MBio::Tree::RandomFactory", "-e", "1"],
            capture_output=True,
            text=True,
            timeout=5,
        )
        return res.returncode == 0
    except Exception:
        return False


def generate_tree_with_bioperl(
    taxa: int = 32,
    scale: float = 1.0,
    seed: Optional[int] = None,
    rate_sd: float = 0.0,
    lba_ratio: float = 1.0,
    outtree: Optional[str] = None,
) -> str:
    """
    Generates a tree strictly using BioPerl Bio::Tree::RandomFactory via bin/generate_tree.pl.
    Raises RuntimeError on failure without fallback.
    """
    script_dir = os.path.dirname(os.path.abspath(__file__))
    pl_script = os.path.join(script_dir, "generate_tree.pl")
    if not os.path.exists(pl_script):
        raise FileNotFoundError(f"BioPerl tree generator script not found at {pl_script}")

    cleanup_temp = False
    if outtree is None:
        fd, target_path = tempfile.mkstemp(suffix=".nwk")
        os.close(fd)
        cleanup_temp = True
    else:
        target_path = outtree

    perl_bin = get_perl_executable() or "perl"
    cmd = [
        perl_bin,
        pl_script,
        "--taxa", str(taxa),
        "--scale", str(scale),
        "--rate_sd", str(rate_sd),
        "--lba_ratio", str(lba_ratio),
        "--outtree", target_path,
    ]
    if seed is not None:
        cmd.extend(["--seed", str(seed)])

    try:
        proc = subprocess.run(cmd, capture_output=True, text=True)
        if proc.returncode != 0:
            err_msg = proc.stderr.strip() or proc.stdout.strip()
            hint = ""
            if "Bio/Tree/RandomFactory.pm" in err_msg or "IO/String.pm" in err_msg:
                hint = "\n[HINT] BioPerl is missing or not accessible. Run: micromamba install -y -c bioconda -c conda-forge perl-bioperl"
            raise RuntimeError(
                f"BioPerl execution failed (code {proc.returncode}) using {perl_bin}:\n{err_msg}{hint}"
            )

        with open(target_path, "r") as f:
            newick = f.read().strip()
        return newick
    finally:
        if cleanup_temp and os.path.exists(target_path):
            try:
                os.remove(target_path)
            except OSError:
                pass


def generate_paper_yule_tree(
    taxa: int = 32,
    scale: float = 1.0,
    seed: Optional[int] = None,
    rate_sd: float = 0.0,
    lba_ratio: float = 1.0,
    min_length: float = 1e-6,
    engine: str = "bioperl",
) -> str:
    """
    Generates a phylogenetic tree based on the backward Yule process matching
    Matsui & Iwasaki (2020, Systematic Biology) and BioPerl Bio::Tree::RandomFactory.

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
    min_length : float, default=1e-6
        Minimum allowable branch length (substitutions per site).
    engine : str, default='bioperl'
        Backend engine: 'bioperl' (strictly runs BioPerl without fallback)
        or 'python' (pure python implementation).

    Returns
    -------
    str
        Valid Newick tree string terminated with a semicolon ';'.
    """
    if taxa < 2:
        raise ValueError(f"Taxa count must be at least 2, got {taxa}")
    if scale <= 0.0:
        raise ValueError(f"Scale must be strictly positive, got {scale}")

    if engine == "bioperl":
        return generate_tree_with_bioperl(
            taxa=taxa,
            scale=scale,
            seed=seed,
            rate_sd=rate_sd,
            lba_ratio=lba_ratio,
        )
    elif engine == "python":
        return _generate_paper_yule_tree_python(
            taxa=taxa,
            scale=scale,
            seed=seed,
            rate_sd=rate_sd,
            lba_ratio=lba_ratio,
            min_length=min_length,
        )
    else:
        raise ValueError(f"Unsupported engine: '{engine}'. Must be 'bioperl' or 'python'.")


def main() -> None:
    """
    Parses CLI arguments, generates the specified tree, and writes to output file.
    """
    parser = argparse.ArgumentParser(
        description="Generate synthetic phylogenetic trees (Matsui & Iwasaki 2020 paper BioPerl model)"
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
        "--engine",
        choices=["auto", "bioperl", "python"],
        default="auto",
        help="Tree generation engine (default: auto; uses BioPerl if available, falls back to Python with warning).",
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

    # Ensure parent directory exists
    out_dir = os.path.dirname(os.path.abspath(args.outtree))
    if out_dir and not os.path.exists(out_dir):
        os.makedirs(out_dir, exist_ok=True)

    # Determine execution engine
    chosen_engine = args.engine
    if chosen_engine == "auto":
        if is_bioperl_available():
            chosen_engine = "bioperl"
        else:
            sys.stderr.write(
                "[WARNING] BioPerl (Bio::Tree::RandomFactory) not found. "
                "Using mathematically equivalent Python backward Yule generator.\n"
            )
            chosen_engine = "python"

    # Generate tree strictly with selected engine
    if chosen_engine == "bioperl":
        generate_tree_with_bioperl(
            taxa=args.taxa,
            scale=args.scale,
            seed=args.seed,
            rate_sd=args.rate_sd,
            lba_ratio=args.lba_ratio,
            outtree=args.outtree,
        )
        print(f"Generated {args.model} tree [BioPerl] ({args.taxa} taxa, scale={args.scale}, seed={args.seed}) -> {args.outtree}")
    elif chosen_engine == "python":
        tree_newick = _generate_paper_yule_tree_python(
            taxa=args.taxa,
            scale=args.scale,
            seed=args.seed,
            rate_sd=args.rate_sd,
            lba_ratio=args.lba_ratio,
        )
        with open(args.outtree, "w") as f:
            f.write(tree_newick + "\n")
        print(f"Generated {args.model} tree [Python] ({args.taxa} taxa, scale={args.scale}, seed={args.seed}) -> {args.outtree}")
    else:
        raise ValueError(f"Unsupported engine: {args.engine}")


if __name__ == "__main__":
    main()

