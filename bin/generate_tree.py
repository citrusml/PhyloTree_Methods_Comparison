#!/usr/bin/env python3
"""
Generate Phylogenetic Tree CLI Script: generate_tree.py

Generates synthetic Newick tree files according to Matsui & Iwasaki (2020, Systematic Biology)
simulation specifications:
1. Primary engine: BioPerl 1.6.924 (Bio::Tree::RandomFactory, rand_yule_c_tree) via bin/generate_tree.pl
2. Topology: Backward Yule process with constant birth rate lambda = 1.0.
3. Branch Lengths: Sampled from logarithmic distribution l = 1 - ln(u * (e - 1) + 1) and scaled by D.
4. NELSI Autocorrelated Models: Yule chronogram (ape::rphylo) + NELSI (simulate.autocor.kishino / simulate.autocor.thorne) via bin/generate_tree_nelsi.R.
5. Taxa: Terminal leaves named T1, T2, ..., TN for AliSim and INDELible compatibility.
6. Robust fallback: Pure Python implementation matching BioPerl specifications if Perl/BioPerl is absent (for paper_yule).
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

    # 0. Explicit BIOPERL_PERL environment variable or dedicated bioperl_env
    bioperl_env_perl = os.environ.get("BIOPERL_PERL")
    if bioperl_env_perl and os.path.isfile(bioperl_env_perl) and os.access(bioperl_env_perl, os.X_OK):
        candidates.append(bioperl_env_perl)

    conda_prefix = os.environ.get("CONDA_PREFIX")
    mamba_root = os.environ.get("MAMBA_ROOT_PREFIX")
    home_dir = os.path.expanduser("~")
    possible_env_roots = []
    if conda_prefix:
        possible_env_roots.append(os.path.dirname(conda_prefix))
    if mamba_root:
        possible_env_roots.append(os.path.join(mamba_root, "envs"))
    possible_env_roots.extend([
        os.path.join(home_dir, ".micromamba", "envs"),
        os.path.join(home_dir, "micromamba", "envs"),
        os.path.join(home_dir, ".conda", "envs"),
        os.path.join(home_dir, "miniconda3", "envs"),
    ])
    for env_root in possible_env_roots:
        bp_perl = os.path.join(env_root, "bioperl_env", "bin", "perl")
        if os.path.isfile(bp_perl) and os.access(bp_perl, os.X_OK) and bp_perl not in candidates:
            candidates.append(bp_perl)

    # 1. Same directory as current python interpreter (conda/micromamba env)
    env_perl = os.path.join(os.path.dirname(sys.executable), "perl")
    if os.path.isfile(env_perl) and os.access(env_perl, os.X_OK) and env_perl not in candidates:
        candidates.append(env_perl)

    # 2. CONDA_PREFIX if set
    if conda_prefix:
        cp_perl = os.path.join(conda_prefix, "bin", "perl")
        if os.path.isfile(cp_perl) and os.access(cp_perl, os.X_OK) and cp_perl not in candidates:
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


def get_rscript_executable() -> Optional[str]:
    """
    Finds an Rscript executable capable of loading the NELSI R package.

    Checks ``NELSI_RSCRIPT``, the active Python/Conda environment, ``PATH``,
    and standard macOS/Linux R installation paths.

    Returns
    -------
    Optional[str]
        Path to an Rscript binary (preferring one with NELSI installed), or None if not found.
    """
    candidates: List[str] = []

    env_rscript = os.environ.get("NELSI_RSCRIPT")
    if env_rscript and os.path.isfile(env_rscript) and os.access(env_rscript, os.X_OK):
        candidates.append(env_rscript)

    local_rscript = os.path.join(os.path.dirname(sys.executable), "Rscript")
    if os.path.isfile(local_rscript) and os.access(local_rscript, os.X_OK) and local_rscript not in candidates:
        candidates.append(local_rscript)

    conda_prefix = os.environ.get("CONDA_PREFIX")
    if conda_prefix:
        cp_rscript = os.path.join(conda_prefix, "bin", "Rscript")
        if os.path.isfile(cp_rscript) and os.access(cp_rscript, os.X_OK) and cp_rscript not in candidates:
            candidates.append(cp_rscript)

    which_rscript = shutil.which("Rscript")
    if which_rscript and which_rscript not in candidates:
        candidates.append(which_rscript)

    for std_path in ("/usr/local/bin/Rscript", "/opt/homebrew/bin/Rscript", "/usr/bin/Rscript"):
        if os.path.isfile(std_path) and os.access(std_path, os.X_OK) and std_path not in candidates:
            candidates.append(std_path)

    for rbin in candidates:
        try:
            res = subprocess.run(
                [rbin, "-e", "if (!requireNamespace('NELSI', quietly=TRUE)) quit(status=1)"],
                capture_output=True,
                text=True,
                timeout=10,
            )
            if res.returncode == 0:
                return rbin
        except Exception:
            continue

    return candidates[0] if candidates else None


def is_nelsi_available() -> bool:
    """
    Checks whether Rscript and the NELSI R package are available.

    Returns
    -------
    bool
        True if NELSI can be loaded in Rscript, False otherwise.
    """
    rscript_bin = get_rscript_executable()
    if not rscript_bin:
        return False
    try:
        res = subprocess.run(
            [rscript_bin, "-e", "if (!requireNamespace('NELSI', quietly=TRUE)) quit(status=1)"],
            capture_output=True,
            text=True,
            timeout=10,
        )
        return res.returncode == 0
    except Exception:
        return False


def generate_tree_with_nelsi(
    taxa: int = 32,
    scale: float = 1.0,
    seed: Optional[int] = None,
    model: str = "nelsi_kishino",
    v: float = 0.3,
    initial_rate: float = 1.0,
    scale_mode: str = "root_to_tip",
    lba_ratio: float = 1.0,
    min_length: float = 1e-6,
    outtree: Optional[str] = None,
) -> str:
    """
    Generates a Yule-topology phylogenetic tree with autocorrelated branch lengths
    strictly using the R package NELSI via ``bin/generate_tree_nelsi.R`` (no fallback).

    Parameters
    ----------
    taxa : int, default=32
        Number of terminal taxa (leaves) in the tree [taxa]. Must be >= 2.
    scale : float, default=1.0
        Evolutionary distance scaling factor D [substitutions/site]. Must be > 0.
    seed : Optional[int], default=None
        Random seed for R RNG reproducibility [unitless].
    model : str, default='nelsi_kishino'
        NELSI autocorrelation model: ``'nelsi_kishino'`` (or ``'nelsi'``, ``'kishino'``)
        or ``'nelsi_thorne'`` (or ``'thorne'``).
    v : float, default=0.3
        NELSI autocorrelation variance parameter v [unitless]. Must be >= 0.
    initial_rate : float, default=1.0
        NELSI initial root substitution rate [substitutions/site/unit time]. Must be > 0.
    scale_mode : str, default='root_to_tip'
        How ``scale`` is applied to the NELSI phylogram: ``'root_to_tip'`` (mean
        root-to-tip distance equals ``scale``), ``'mean_branch'`` (mean branch length
        equals ``scale``), or ``'raw'`` (multiplies raw NELSI lengths by ``scale``).
    lba_ratio : float, default=1.0
        Long-branch attraction multiplier applied to T1 and TN [unitless].
    min_length : float, default=1e-6
        Minimum branch length floor [substitutions/site].
    outtree : Optional[str], default=None
        Optional output Newick file path. If None, uses a temporary file.

    Returns
    -------
    str
        Valid Newick tree string terminated with a semicolon ';'.
    """
    if taxa < 2:
        raise ValueError(f"Taxa count must be at least 2, got {taxa}")
    if scale <= 0.0:
        raise ValueError(f"Scale must be strictly positive, got {scale}")
    if v < 0.0:
        raise ValueError(f"NELSI autocorrelation parameter v must be non-negative, got {v}")
    if initial_rate <= 0.0:
        raise ValueError(f"NELSI initial_rate must be strictly positive, got {initial_rate}")
    if scale_mode not in ("root_to_tip", "mean_branch", "raw"):
        raise ValueError(f"Unsupported scale_mode '{scale_mode}'. Must be 'root_to_tip', 'mean_branch', or 'raw'.")

    script_dir = os.path.dirname(os.path.abspath(__file__))
    r_script = os.path.join(script_dir, "generate_tree_nelsi.R")
    if not os.path.exists(r_script):
        raise FileNotFoundError(f"NELSI R tree generator script not found at {r_script}")

    rscript_bin = get_rscript_executable()
    if not rscript_bin:
        raise RuntimeError(
            "Rscript executable not found. Install R (>= 3.6.0) along with packages "
            "'ape', 'castor', 'digest', and 'NELSI' (sebastianduchene/NELSI)."
        )

    cleanup_temp = False
    if outtree is None:
        fd, target_path = tempfile.mkstemp(suffix=".nwk")
        os.close(fd)
        cleanup_temp = True
    else:
        target_path = outtree

    cmd = [
        rscript_bin,
        r_script,
        "--taxa", str(taxa),
        "--scale", str(scale),
        "--model", str(model),
        "--v", str(v),
        "--initial_rate", str(initial_rate),
        "--scale_mode", str(scale_mode),
        "--lba_ratio", str(lba_ratio),
        "--min_length", str(min_length),
        "--outtree", target_path,
    ]
    if seed is not None:
        cmd.extend(["--seed", str(seed)])

    try:
        proc = subprocess.run(cmd, capture_output=True, text=True)
        if proc.returncode != 0:
            err_msg = proc.stderr.strip() or proc.stdout.strip()
            raise RuntimeError(
                f"NELSI R execution failed (code {proc.returncode}) using {rscript_bin}:\n{err_msg}"
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


def generate_nelsi_yule_tree(
    taxa: int = 32,
    scale: float = 1.0,
    seed: Optional[int] = None,
    model: str = "nelsi_kishino",
    v: float = 0.3,
    initial_rate: float = 1.0,
    scale_mode: str = "root_to_tip",
    lba_ratio: float = 1.0,
    min_length: float = 1e-6,
) -> str:
    """
    Public wrapper to generate a Yule-topology tree with NELSI autocorrelated branch lengths.

    Parameters
    ----------
    taxa : int, default=32
        Number of terminal taxa [taxa].
    scale : float, default=1.0
        Evolutionary distance scaling factor D [substitutions/site].
    seed : Optional[int], default=None
        Random seed for reproducibility [unitless].
    model : str, default='nelsi_kishino'
        NELSI model: ``'nelsi_kishino'`` or ``'nelsi_thorne'``.
    v : float, default=0.3
        NELSI autocorrelation variance parameter v [unitless].
    initial_rate : float, default=1.0
        NELSI initial root substitution rate [substitutions/site/unit time].
    scale_mode : str, default='root_to_tip'
        Scaling mode: ``'root_to_tip'``, ``'mean_branch'``, or ``'raw'``.
    lba_ratio : float, default=1.0
        Long-branch attraction ratio b/a [unitless].
    min_length : float, default=1e-6
        Minimum branch length [substitutions/site].

    Returns
    -------
    str
        Valid Newick tree string terminated with ';'.
    """
    return generate_tree_with_nelsi(
        taxa=taxa,
        scale=scale,
        seed=seed,
        model=model,
        v=v,
        initial_rate=initial_rate,
        scale_mode=scale_mode,
        lba_ratio=lba_ratio,
        min_length=min_length,
    )


def main() -> None:
    """
    Parses CLI arguments, generates the specified tree, and writes to output file.
    """
    parser = argparse.ArgumentParser(
        description="Generate synthetic phylogenetic trees (BioPerl paper_yule or R NELSI autocorrelated Yule)"
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
        choices=["paper_yule", "yule", "nelsi_kishino", "nelsi_thorne", "nelsi"],
        default="paper_yule",
        help="Tree generation model (default: paper_yule)",
    )
    parser.add_argument(
        "--engine",
        choices=["auto", "bioperl", "python", "nelsi"],
        default="auto",
        help="Tree generation engine (default: auto; uses BioPerl for paper_yule or R NELSI for nelsi_* models).",
    )
    parser.add_argument(
        "--rate_sd",
        type=float,
        default=0.0,
        help="Rate heterogeneity standard deviation for paper_yule (default: 0.0)",
    )
    parser.add_argument(
        "--nelsi_v",
        type=float,
        default=0.3,
        help="NELSI autocorrelation variance parameter v (default: 0.3)",
    )
    parser.add_argument(
        "--nelsi_initial_rate",
        type=float,
        default=1.0,
        help="NELSI initial root substitution rate (default: 1.0)",
    )
    parser.add_argument(
        "--nelsi_scale_mode",
        choices=["root_to_tip", "mean_branch", "raw"],
        default="root_to_tip",
        help="Scaling mode for NELSI phylogram relative to --scale (default: root_to_tip)",
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

    is_nelsi_model = args.model in ("nelsi_kishino", "nelsi_thorne", "nelsi")

    # Determine execution engine
    chosen_engine = args.engine
    if is_nelsi_model:
        if chosen_engine in ("bioperl", "python"):
            raise ValueError(
                f"Engine '{chosen_engine}' is not compatible with NELSI model '{args.model}'. "
                "Use '--engine nelsi' or '--engine auto'."
            )
        chosen_engine = "nelsi"
    elif chosen_engine == "nelsi":
        raise ValueError(
            f"Engine 'nelsi' requires a NELSI model (--model nelsi_kishino or --model nelsi_thorne), "
            f"got '--model {args.model}'."
        )
    elif chosen_engine == "auto":
        if is_bioperl_available():
            chosen_engine = "bioperl"
        else:
            sys.stderr.write(
                "[WARNING] BioPerl (Bio::Tree::RandomFactory) not found. "
                "Using mathematically equivalent Python backward Yule generator.\n"
            )
            chosen_engine = "python"

    # Generate tree strictly with selected engine
    if chosen_engine == "nelsi":
        generate_tree_with_nelsi(
            taxa=args.taxa,
            scale=args.scale,
            seed=args.seed,
            model=args.model,
            v=args.nelsi_v,
            initial_rate=args.nelsi_initial_rate,
            scale_mode=args.nelsi_scale_mode,
            lba_ratio=args.lba_ratio,
            outtree=args.outtree,
        )
        print(
            f"Generated {args.model} tree [NELSI] ({args.taxa} taxa, scale={args.scale}, "
            f"v={args.nelsi_v}, scale_mode={args.nelsi_scale_mode}, seed={args.seed}) -> {args.outtree}"
        )
    elif chosen_engine == "bioperl":
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

