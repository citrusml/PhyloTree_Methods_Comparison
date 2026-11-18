#!/usr/bin/env python3
"""
Evaluate Long-Branch Attraction (LBA)-specific metrics on Nextflow replication
directories (e.g. results_nelsi_yule, results_nelsi_lba3, results_nelsi_lba5).

Metrics implemented (following JoiningBasedOnAdditivity/src/evaluation.py):
  - rf, nrf: Standard unrooted Robinson-Foulds distance and normalized RF [unitless]
  - ds_true, ds_est, delta_ds: Shortest topological edge path between LBA pair (u, v) [edges]
  - is_sister: 1.0 if (u, v) are inferred as immediate sisters (ds_est == 2), else 0.0 [unitless]
  - nta: Normalized Topological Attraction = (ds_true - ds_est) / (ds_true - 2) [unitless]
  - pair_nrf, pair_split_recall, pair_split_fnr: RF restricted to splits separating (u, v) [unitless]
  - lba_quartet_attraction: Fraction of separated quartets (u,v,k,l) attracted to uv|kl [unitless]
  - lba_quartet_error: Fraction of all quartets (u,v,k,l) mismatched vs true tree [unitless]
  - matrix_lba_margin: Mean normalized quartet margin in [-1, 1] from PHYLIP distance matrix [unitless]
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
import itertools
from pathlib import Path
import re
from typing import Dict, List, Optional, Set, Tuple

import dendropy
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


PIPELINE_FILE_MAP: Dict[str, Tuple[str, Optional[str]]] = {
    "PSA+NJ": ("pwa_nj.nwk", "pwa_matrix.phylip"),
    "PSA+NJ+NNI": ("pwa_nj_nni.nwk", "pwa_matrix.phylip"),
    "MMSEQS+NJ": ("mmseqs_nj.nwk", "mmseqs_matrix.phylip"),
    "MMSEQS+NJ+NNI": ("mmseqs_nj_nni.nwk", "mmseqs_matrix.phylip"),
    "MSA+NJ": ("msa_nj.nwk", "msa_matrix.phylip"),
    "MSA+NJ+NNI": ("msa_nj_nni.nwk", "msa_matrix.phylip"),
    "MSA+ML": ("msa_ml.nwk", None),
    "MSA+RAXML": ("msa_raxml.nwk", None),
    "GS": ("gs.nwk", None),
    "TRUE_PSA+NJ": ("true_pwa_nj.nwk", "true_pwa_matrix.phylip"),
    "TRUE_PSA+NJ+NNI": ("true_pwa_nj_nni.nwk", "true_pwa_matrix.phylip"),
    "TRUE_MSA+NJ": ("true_msa_nj.nwk", "true_msa_matrix.phylip"),
    "TRUE_MSA+NJ+NNI": ("true_msa_nj_nni.nwk", "true_msa_matrix.phylip"),
    "TRUE_MSA+ML": ("true_msa_ml.nwk", None),
    "TRUE_MSA+RAXML": ("true_msa_raxml.nwk", None),
}

_DIR_PATTERN = re.compile(r"^D(?P<dist>[0-9.]+)_L(?P<len>\d+)_rep(?P<rep>\d+)$")


def _parse_unrooted_tree(tree_newick: str, tns: dendropy.TaxonNamespace) -> dendropy.Tree:
    """
    Parse a Newick string into an unrooted DendroPy tree with encoded bipartitions.

    Parameters
    ----------
    tree_newick : str
        Newick tree string (unitless).
    tns : dendropy.TaxonNamespace
        Shared TaxonNamespace across trees being compared.

    Returns
    -------
    dendropy.Tree
        Unrooted DendroPy tree with bipartitions encoded.
    """
    tree = dendropy.Tree.get(
        data=tree_newick,
        schema="newick",
        taxon_namespace=tns,
        preserve_underscores=True,
    )
    tree.is_rooted = False
    tree.deroot()
    tree.encode_bipartitions()
    return tree


def _compute_tree_patristic_matrix(
    tree: dendropy.Tree,
    taxa_names: List[str],
    unit_lengths: bool = False,
) -> np.ndarray:
    """
    Compute the pairwise patristic (or topological edge-count) distance matrix.

    Parameters
    ----------
    tree : dendropy.Tree
        DendroPy tree object.
    taxa_names : List[str]
        Ordered list of terminal taxon labels.
    unit_lengths : bool, optional
        If True, treat every edge as length 1.0 (topological distance, [edges]).
        If False, use branch lengths ([substitutions/site]).

    Returns
    -------
    np.ndarray
        Symmetric (N, N) distance matrix.
    """
    n_taxa = len(taxa_names)
    name_to_idx = {name: i for i, name in enumerate(taxa_names)}
    D = np.zeros((n_taxa, n_taxa), dtype=np.float64)

    work_tree = tree.clone(depth=1)
    if unit_lengths:
        for edge in work_tree.preorder_edge_iter():
            if edge.tail_node is not None:
                edge.length = 1.0
    else:
        for edge in work_tree.preorder_edge_iter():
            if edge.tail_node is not None and edge.length is None:
                edge.length = 0.0

    pdm = work_tree.phylogenetic_distance_matrix()
    taxa = list(work_tree.taxon_namespace)
    for i, t_i in enumerate(taxa):
        idx_i = name_to_idx.get(t_i.label)
        if idx_i is None:
            continue
        for j in range(i + 1, len(taxa)):
            t_j = taxa[j]
            idx_j = name_to_idx.get(t_j.label)
            if idx_j is None:
                continue
            val = float(pdm.distance(t_i, t_j))
            D[idx_i, idx_j] = val
            D[idx_j, idx_i] = val
    return D


def _extract_pair_separating_splits(
    tree: dendropy.Tree,
    mask_u: int,
    mask_v: int,
    all_ones: int,
) -> Set[int]:
    """
    Extract canonical split bitmasks of nontrivial internal edges separating (u, v).

    Parameters
    ----------
    tree : dendropy.Tree
        Unrooted DendroPy tree with bipartitions encoded.
    mask_u : int
        Taxon bitmask for first target taxon.
    mask_v : int
        Taxon bitmask for second target taxon.
    all_ones : int
        Full bitmask `(1 << n_taxa) - 1`.

    Returns
    -------
    Set[int]
        Canonical bitmasks for internal splits separating u and v.
    """
    splits: Set[int] = set()
    for edge in tree.postorder_edge_iter():
        if edge.bipartition is None:
            continue
        bm = edge.bipartition.split_bitmask
        if bm == 0 or bm == all_ones:
            continue
        pop = bin(bm).count("1")
        if pop <= 1 or pop >= bin(all_ones).count("1") - 1:
            continue
        if bool(bm & mask_u) != bool(bm & mask_v):
            splits.add(min(bm, all_ones ^ bm))
    return splits


def evaluate_tree_lba(
    est_tree_newick: str,
    true_tree_newick: str,
    taxa_names: List[str],
    lba_pair: Optional[Tuple[str, str]] = None,
    auto_max_pair: bool = False,
) -> Dict[str, float]:
    """
    Evaluate global RF/nRF and LBA-specific topological attraction metrics for
    a target pair (u, v) between an estimated tree and a true tree.

    Parameters
    ----------
    est_tree_newick : str
        Newick string of the inferred tree.
    true_tree_newick : str
        Newick string of the ground-truth tree.
    taxa_names : List[str]
        Complete ordered list of leaf taxon labels.
    lba_pair : Optional[Tuple[str, str]], optional
        Target LBA taxon pair (u, v). Defaults to (taxa_names[0], taxa_names[-1])
        unless `auto_max_pair=True`.
    auto_max_pair : bool, optional
        If True and `lba_pair` is None, select the pair maximizing true patristic
        distance [substitutions/site].

    Returns
    -------
    Dict[str, float]
        Dictionary of RF, topological edge path, sisterhood, pair-split RF,
        and quartet attraction metrics.
    """
    n_taxa = len(taxa_names)
    tns = dendropy.TaxonNamespace(taxa_names)
    t_true = _parse_unrooted_tree(true_tree_newick, tns)
    t_est = _parse_unrooted_tree(est_tree_newick, tns)

    rf = float(dendropy.calculate.treecompare.symmetric_difference(t_true, t_est))
    max_rf = float(2 * (n_taxa - 3)) if n_taxa > 3 else 1.0
    nrf = rf / max_rf if max_rf > 0.0 else 0.0

    if lba_pair is not None:
        t1, t2 = lba_pair
    elif auto_max_pair:
        D_pat = _compute_tree_patristic_matrix(t_true, taxa_names, unit_lengths=False)
        triu_i, triu_j = np.triu_indices(n_taxa, k=1)
        max_flat = int(np.argmax(D_pat[triu_i, triu_j]))
        t1 = taxa_names[int(triu_i[max_flat])]
        t2 = taxa_names[int(triu_j[max_flat])]
    else:
        t1, t2 = taxa_names[0], taxa_names[-1]

    name_to_idx = {name: i for i, name in enumerate(taxa_names)}
    idx_u = name_to_idx[t1]
    idx_v = name_to_idx[t2]

    D_topo_true = _compute_tree_patristic_matrix(t_true, taxa_names, unit_lengths=True)
    D_topo_est = _compute_tree_patristic_matrix(t_est, taxa_names, unit_lengths=True)

    ds_true = int(round(D_topo_true[idx_u, idx_v]))
    ds_est = int(round(D_topo_est[idx_u, idx_v]))
    delta_ds = ds_true - ds_est
    is_sister = 1.0 if ds_est == 2 else 0.0
    nta = float(delta_ds) / float(ds_true - 2) if ds_true > 2 else 0.0

    # Pair-separating split metrics
    mask_u = tns.taxon_bitmask(tns.get_taxon(t1))
    mask_v = tns.taxon_bitmask(tns.get_taxon(t2))
    all_ones = (1 << n_taxa) - 1
    s_true = _extract_pair_separating_splits(t_true, mask_u, mask_v, all_ones)
    s_est = _extract_pair_separating_splits(t_est, mask_u, mask_v, all_ones)

    n_true_s = len(s_true)
    n_est_s = len(s_est)
    n_shared_s = len(s_true & s_est)
    pair_recall = float(n_shared_s) / float(n_true_s) if n_true_s > 0 else 1.0
    pair_fnr = 1.0 - pair_recall
    denom_s = n_true_s + n_est_s
    pair_nrf = float(len(s_true ^ s_est)) / float(denom_s) if denom_s > 0 else 0.0

    # LBA Quartet Attraction metrics over all (n-2 choose 2) quartets (u, v, k, l)
    others = [x for x in range(n_taxa) if x != idx_u and x != idx_v]
    if len(others) >= 2:
        pairs_kl = np.array(list(itertools.combinations(others, 2)), dtype=np.int64)
        k = pairs_kl[:, 0]
        l = pairs_kl[:, 1]

        st1 = D_topo_true[idx_u, idx_v] + D_topo_true[k, l]
        st2 = D_topo_true[idx_u, k] + D_topo_true[idx_v, l]
        st3 = D_topo_true[idx_u, l] + D_topo_true[idx_v, k]
        true_choice = np.argmin(np.stack([st1, st2, st3], axis=1), axis=1)
        separated_mask = st1 > np.minimum(st2, st3)

        se1 = D_topo_est[idx_u, idx_v] + D_topo_est[k, l]
        se2 = D_topo_est[idx_u, k] + D_topo_est[idx_v, l]
        se3 = D_topo_est[idx_u, l] + D_topo_est[idx_v, k]
        est_choice = np.argmin(np.stack([se1, se2, se3], axis=1), axis=1)
        est_attracted = se1 < np.minimum(se2, se3)

        n_sep = int(np.sum(separated_mask))
        quartet_attraction = float(np.mean(est_attracted[separated_mask])) if n_sep > 0 else 0.0
        quartet_error = float(np.mean(est_choice != true_choice))
    else:
        quartet_attraction = 0.0
        quartet_error = 0.0

    return {
        "rf": rf,
        "nrf": nrf,
        "ds_est": float(ds_est),
        "ds_true": float(ds_true),
        "delta_ds": float(delta_ds),
        "is_sister": is_sister,
        "nta": nta,
        "pair_nrf": pair_nrf,
        "pair_split_recall": pair_recall,
        "pair_split_fnr": pair_fnr,
        "lba_quartet_attraction": quartet_attraction,
        "lba_quartet_error": quartet_error,
    }


def _read_phylip_matrix(matrix_path: Path) -> Tuple[List[str], np.ndarray]:
    """
    Read a square PHYLIP distance matrix file.

    Parameters
    ----------
    matrix_path : Path
        Path to PHYLIP distance matrix.

    Returns
    -------
    Tuple[List[str], np.ndarray]
        Ordered taxon labels and (N, N) distance matrix [substitutions/site].
    """
    lines = [ln.strip() for ln in matrix_path.read_text(encoding="utf-8").splitlines() if ln.strip()]
    n_taxa = int(lines[0].split()[0])
    names: List[str] = []
    mat = np.zeros((n_taxa, n_taxa), dtype=np.float64)
    for i, ln in enumerate(lines[1 : n_taxa + 1]):
        parts = ln.split()
        names.append(parts[0])
        mat[i, :n_taxa] = [max(0.0, float(x)) for x in parts[1 : n_taxa + 1]]
    return names, mat


def compute_matrix_lba_margin(
    matrix_path: Path,
    true_tree_newick: str,
    taxa_names: List[str],
    lba_pair: Tuple[str, str],
    eps: float = 1e-7,
) -> float:
    """
    Compute the mean normalized four-point quartet margin for separated quartets
    (u, v, k, l) from a PHYLIP distance matrix.

    Positive values indicate the distance matrix separates (u, v); negative
    values indicate entry into the Felsenstein attraction zone.

    Parameters
    ----------
    matrix_path : Path
        Path to PHYLIP distance matrix.
    true_tree_newick : str
        Newick string of ground-truth tree.
    taxa_names : List[str]
        Ordered taxon labels matching the true tree.
    lba_pair : Tuple[str, str]
        Target LBA taxon pair (u, v).
    eps : float, optional
        Numerical stabilizer (default: 1e-7).

    Returns
    -------
    float
        Mean quartet margin in [-1, 1] [unitless].
    """
    mat_names, D_raw = _read_phylip_matrix(matrix_path)
    mat_idx = {nm: i for i, nm in enumerate(mat_names)}
    order = [mat_idx[nm] for nm in taxa_names]
    D = D_raw[np.ix_(order, order)]

    n_taxa = len(taxa_names)
    name_to_idx = {nm: i for i, nm in enumerate(taxa_names)}
    idx_u = name_to_idx[lba_pair[0]]
    idx_v = name_to_idx[lba_pair[1]]

    others = [x for x in range(n_taxa) if x != idx_u and x != idx_v]
    if len(others) < 2:
        return 0.0

    tns = dendropy.TaxonNamespace(taxa_names)
    t_true = _parse_unrooted_tree(true_tree_newick, tns)
    D_topo_true = _compute_tree_patristic_matrix(t_true, taxa_names, unit_lengths=True)

    pairs_kl = np.array(list(itertools.combinations(others, 2)), dtype=np.int64)
    k = pairs_kl[:, 0]
    l = pairs_kl[:, 1]
    st1 = D_topo_true[idx_u, idx_v] + D_topo_true[k, l]
    st2 = D_topo_true[idx_u, k] + D_topo_true[idx_v, l]
    st3 = D_topo_true[idx_u, l] + D_topo_true[idx_v, k]
    separated_mask = st1 > np.minimum(st2, st3)
    if int(np.sum(separated_mask)) == 0:
        return 0.0

    ks, ls = k[separated_mask], l[separated_mask]
    s1 = D[idx_u, idx_v] + D[ks, ls]
    s2 = D[idx_u, ks] + D[idx_v, ls]
    s3 = D[idx_u, ls] + D[idx_v, ks]
    s_max = np.maximum(s1, np.maximum(s2, s3))
    margin = (s1 - np.minimum(s2, s3)) / (s_max + eps)
    return float(np.mean(margin))


def _evaluate_rep_dir(
    task: Tuple[str, Path, Tuple[str, str], bool],
) -> List[Dict[str, object]]:
    """
    Evaluate all available inferred trees in a single `D*_L*_rep*` directory.

    Parameters
    ----------
    task : Tuple[str, Path, Tuple[str, str], bool]
        Tuple of (experiment_label, rep_dir, lba_pair, auto_max_pair).

    Returns
    -------
    List[Dict[str, object]]
        Per-pipeline metric records for this replicate.
    """
    exp_label, rep_dir, lba_pair, auto_max_pair = task
    true_path = rep_dir / "true_tree.nwk"
    if not true_path.is_file():
        return []

    m = _DIR_PATTERN.match(rep_dir.name)
    if not m:
        return []
    dist = float(m.group("dist"))
    length = int(m.group("len"))
    rep = int(m.group("rep"))

    true_nwk = true_path.read_text(encoding="utf-8").strip()
    t_tmp = dendropy.Tree.get(data=true_nwk, schema="newick", preserve_underscores=True)
    taxa_names = sorted(
        [lf.taxon.label for lf in t_tmp.leaf_nodes()],
        key=lambda x: (int(x[1:]) if x.startswith("T") and x[1:].isdigit() else 999999, x),
    )

    pair_arg = None if auto_max_pair else lba_pair
    if not auto_max_pair and (lba_pair[0] not in taxa_names or lba_pair[1] not in taxa_names):
        pair_arg = (taxa_names[0], taxa_names[-1])

    rows: List[Dict[str, object]] = []
    margin_cache: Dict[str, float] = {}

    for pipeline_name, (tree_fname, mat_fname) in PIPELINE_FILE_MAP.items():
        tree_path = rep_dir / tree_fname
        if not tree_path.is_file() or tree_path.stat().st_size == 0:
            continue

        est_nwk = tree_path.read_text(encoding="utf-8").strip()
        metrics = evaluate_tree_lba(
            est_tree_newick=est_nwk,
            true_tree_newick=true_nwk,
            taxa_names=taxa_names,
            lba_pair=pair_arg,
            auto_max_pair=auto_max_pair,
        )

        mat_margin = np.nan
        if mat_fname is not None:
            mat_path = rep_dir / mat_fname
            if mat_path.is_file() and mat_path.stat().st_size > 0:
                if mat_fname not in margin_cache:
                    eff_pair = pair_arg if pair_arg is not None else (taxa_names[0], taxa_names[-1])
                    margin_cache[mat_fname] = compute_matrix_lba_margin(
                        mat_path, true_nwk, taxa_names, eff_pair
                    )
                mat_margin = margin_cache[mat_fname]

        rows.append({
            "experiment": exp_label,
            "num_taxa": len(taxa_names),
            "distance": dist,
            "length": length,
            "replicate": rep,
            "pipeline": pipeline_name,
            "matrix_lba_margin": mat_margin,
            **metrics,
        })
    return rows


def _discover_rep_dirs(result_dir: Path) -> List[Path]:
    """
    Locate all `D*_L*_rep*` directories under `<result_dir>/replications` or `<result_dir>`.

    Parameters
    ----------
    result_dir : Path
        Root experiment output directory.

    Returns
    -------
    List[Path]
        Sorted list of replicate directories.
    """
    rep_root = result_dir / "replications"
    base = rep_root if rep_root.is_dir() else result_dir
    if not base.is_dir():
        return []
    return sorted([p for p in base.iterdir() if p.is_dir() and _DIR_PATTERN.match(p.name)])


def plot_lba_summary(df: pd.DataFrame, plot_dir: Path) -> None:
    """
    Generate summary comparison plots for LBA-specific metrics across experiments
    and inference pipelines.

    Parameters
    ----------
    df : pd.DataFrame
        Per-replicate LBA metrics dataframe.
    plot_dir : Path
        Directory to save output PNGs.

    Returns
    -------
    None
    """
    plot_dir.mkdir(parents=True, exist_ok=True)
    key_pipes = [
        p for p in ["PSA+NJ", "MMSEQS+NJ", "MSA+NJ", "MSA+ML", "MSA+RAXML", "GS"]
        if p in df["pipeline"].unique()
    ]
    sub = df[df["pipeline"].isin(key_pipes)].copy()
    if sub.empty:
        return

    metrics_to_plot = [
        ("is_sister", "False Sister Rate (ds_est == 2)"),
        ("nta", "Normalized Topological Attraction (NTA)"),
        ("lba_quartet_attraction", "LBA Quartet Attraction Rate"),
        ("pair_nrf", "Pair-Separating Split nRF"),
    ]

    fig, axes = plt.subplots(2, 2, figsize=(16, 11), sharex=True)
    for ax, (col, title_txt) in zip(axes.ravel(), metrics_to_plot):
        sns.barplot(
            data=sub,
            x="pipeline",
            y=col,
            hue="experiment",
            order=key_pipes,
            errorbar=("ci", 95),
            ax=ax,
        )
        ax.set_title(title_txt, fontweight="bold", fontsize=11.5)
        ax.set_xlabel("")
        ax.set_ylabel(col, fontweight="bold")
        ax.tick_params(axis="x", rotation=18)
        ax.grid(True, axis="y", linestyle="--", alpha=0.4)

    fig.suptitle("LBA-Specific Topological Evaluation Across Experiments", fontweight="bold", fontsize=13.5, y=1.01)
    plt.tight_layout()
    fig.savefig(plot_dir / "lba_specific_metrics_comparison.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    """
    Parse command-line arguments.

    Returns
    -------
    argparse.Namespace
        Parsed CLI arguments.
    """
    parser = argparse.ArgumentParser(
        description="Evaluate LBA-specific topological metrics on Nextflow replications."
    )
    parser.add_argument(
        "--results_dirs",
        nargs="+",
        type=Path,
        required=True,
        help="One or more experiment directories (e.g. results/results_nelsi_yule results/results_nelsi_lba3).",
    )
    parser.add_argument(
        "--lba_pair",
        nargs=2,
        default=["T1", "T20"],
        help="Target LBA taxon pair (default: T1 T20).",
    )
    parser.add_argument(
        "--auto_max_pair",
        action="store_true",
        help="Automatically choose the pair with maximum true patristic distance instead of --lba_pair.",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=4,
        help="Number of parallel worker processes (default: 4).",
    )
    parser.add_argument(
        "--outcsv",
        type=Path,
        default=Path("results/nelsi_lba_specific_metrics.csv"),
        help="Output CSV path for per-replicate LBA metrics.",
    )
    parser.add_argument(
        "--summary_csv",
        type=Path,
        default=Path("results/nelsi_lba_specific_summary.csv"),
        help="Output CSV path for aggregated summary by experiment x pipeline x distance.",
    )
    parser.add_argument(
        "--plot_dir",
        type=Path,
        default=Path("results/nelsi_lba_specific_plots"),
        help="Output directory for LBA metric comparison PNGs.",
    )
    return parser.parse_args()


def main() -> None:
    """
    CLI entry point for LBA-specific evaluation across replication directories.

    Returns
    -------
    None
    """
    args = parse_args()
    lba_pair = (args.lba_pair[0], args.lba_pair[1])

    tasks: List[Tuple[str, Path, Tuple[str, str], bool]] = []
    for rdir in args.results_dirs:
        rep_dirs = _discover_rep_dirs(rdir)
        print(f"Discovered {len(rep_dirs)} replicate directories in {rdir}")
        for rp in rep_dirs:
            tasks.append((rdir.name, rp, lba_pair, args.auto_max_pair))

    if not tasks:
        raise FileNotFoundError("No replication directories (D*_L*_rep*) found in --results_dirs.")

    all_rows: List[Dict[str, object]] = []
    if args.workers > 1 and len(tasks) > 1:
        with ProcessPoolExecutor(max_workers=args.workers) as ex:
            for chunk_rows in ex.map(_evaluate_rep_dir, tasks, chunksize=20):
                all_rows.extend(chunk_rows)
    else:
        for t in tasks:
            all_rows.extend(_evaluate_rep_dir(t))

    df = pd.DataFrame(all_rows)
    df.sort_values(by=["experiment", "distance", "length", "replicate", "pipeline"], inplace=True)

    args.outcsv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.outcsv, index=False)
    print(f"Saved per-replicate LBA metrics ({len(df)} rows) -> {args.outcsv}")

    metric_cols = [
        "nrf", "ds_est", "ds_true", "delta_ds", "is_sister", "nta",
        "pair_nrf", "pair_split_recall", "pair_split_fnr",
        "lba_quartet_attraction", "lba_quartet_error", "matrix_lba_margin",
    ]
    summary = (
        df.groupby(["experiment", "pipeline", "distance"], as_index=False)[metric_cols]
        .mean()
    )
    summary.to_csv(args.summary_csv, index=False)
    print(f"Saved summary CSV -> {args.summary_csv}")

    plot_lba_summary(df, args.plot_dir)
    print(f"Saved LBA plots -> {args.plot_dir}")


if __name__ == "__main__":
    main()
