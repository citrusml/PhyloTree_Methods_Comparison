#!/usr/bin/env python3
"""
Evaluate additivity-based and FastME distance methods on existing Nextflow
replication directories without re-running sequence simulation, alignment, ML,
or Standard NJ.

For each replication directory (`replications/D*_L*_rep*/`), reads the three
pre-computed PHYLIP distance matrices:
  - `mmseqs_matrix.phylip`     (sss_linear -> MMSEQS+)
  - `mmseqs_log_matrix.phylip` (sss_log    -> MMSEQS_LOG+)
  - `pwa_matrix.phylip`        (poisson    -> PSA+)

And runs the 7 distance-based tree inference methods:
  - BioNJ
  - FastME_NNI
  - FastME_SPR
  - AddNJ2
  - AddNJ3
  - AddNJ4
  - BioNJ_Add_Hybrid2

Outputs:
  - `<results_dir>/additivity_methods_summary.csv`
  - `<results_dir>/combined_benchmark_summary.csv` (when `benchmark_summary.csv` exists)
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
import os
from pathlib import Path
import re
import sys
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

# Ensure repo root is importable and Python's bin dir (with fastme) is on PATH
_REPO_ROOT = Path(__file__).resolve().parent.parent
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

_PY_BIN = str(Path(sys.executable).resolve().parent)
if _PY_BIN not in os.environ.get("PATH", "").split(os.pathsep):
    os.environ["PATH"] = f"{_PY_BIN}{os.pathsep}{os.environ.get('PATH', '')}"

from analysis.evaluate_lba_metrics import _read_phylip_matrix, evaluate_tree_lba
from bin.additivity.add_nj import add_nj2, add_nj3, add_nj4
from bin.additivity.fastme_wrapper import run_fastme
from bin.additivity.topology_search import TreeDistanceScorer, hill_climb_search

_DIR_PATTERN = re.compile(r"^D(?P<dist>[0-9.]+)_L(?P<len>\d+)_rep(?P<rep>\d+)$")

DIST_MATRICES: Dict[str, Tuple[str, str]] = {
    "sss_linear": ("mmseqs_matrix.phylip", "MMSEQS+"),
    "sss_log": ("mmseqs_log_matrix.phylip", "MMSEQS_LOG+"),
    "poisson": ("pwa_matrix.phylip", "PSA+"),
}

ADDITIVITY_METHODS: List[str] = [
    "BioNJ",
    "FastME_NNI",
    "FastME_SPR",
    "AddNJ2",
    "AddNJ3",
    "AddNJ4",
    "BioNJ_Add_Hybrid2",
]


def infer_tree_from_matrix(
    method: str,
    dist_matrix: np.ndarray,
    taxa_names: List[str],
    max_iter: int = 20,
) -> str:
    """
    Infer an unrooted phylogenetic tree Newick string from a distance matrix.

    Parameters
    ----------
    method : str
        Inference method name (`BioNJ`, `FastME_NNI`, `FastME_SPR`, `AddNJ2`,
        `AddNJ3`, `AddNJ4`, or `BioNJ_Add_Hybrid2`) (unitless).
    dist_matrix : np.ndarray
        Symmetric pairwise distance matrix of shape `(N, N)` [substitutions/site].
    taxa_names : List[str]
        Ordered terminal taxon labels corresponding to rows/columns (unitless).
    max_iter : int, optional
        Maximum hill-climbing iterations per phase for `BioNJ_Add_Hybrid2`
        (unitless). Default is 20.

    Returns
    -------
    str
        Inferred Newick tree string (unitless).
    """
    fastme_dispatch: Dict[str, str] = {
        "BioNJ": "bionj",
        "FastME_NNI": "nni",
        "FastME_SPR": "spr",
    }
    if method in fastme_dispatch:
        return run_fastme(dist_matrix, taxa_names, method=fastme_dispatch[method])

    addnj_dispatch = {
        "AddNJ2": add_nj2,
        "AddNJ3": add_nj3,
        "AddNJ4": add_nj4,
    }
    if method in addnj_dispatch:
        return addnj_dispatch[method](dist_matrix, taxa_names)

    if method == "BioNJ_Add_Hybrid2":
        bionj_nwk = run_fastme(dist_matrix, taxa_names, method="bionj")
        scorer = TreeDistanceScorer(dist_matrix, taxa_names, mode="ols_rmse")
        best_nwk, _ = hill_climb_search(
            bionj_nwk, scorer, mode="hybrid", max_iter=max_iter
        )
        return best_nwk

    raise ValueError(f"Unsupported method: {method}")


def _process_one_replication(
    task: Tuple[Path, str, int, Optional[Tuple[str, str]], bool, bool, bool],
) -> List[Dict[str, object]]:
    """
    Process a single replication directory across all 3 distance matrices and 7 methods.

    Parameters
    ----------
    task : Tuple[Path, str, int, Optional[Tuple[str, str]], bool, bool, bool]
        Tuple of `(rep_dir, condition, max_iter, lba_pair, auto_max_pair,
        save_trees, overwrite)`.

    Returns
    -------
    List[Dict[str, object]]
        List of evaluation metric dictionaries (one per distance x method).
    """
    rep_dir, condition, max_iter, lba_pair, auto_max_pair, save_trees, overwrite = task
    match = _DIR_PATTERN.match(rep_dir.name)
    if match is None:
        return []

    dist_val = float(match.group("dist"))
    length_val = int(match.group("len"))
    rep_val = int(match.group("rep"))

    true_tree_path = rep_dir / "true_tree.nwk"
    if not true_tree_path.exists():
        return []
    true_nwk = true_tree_path.read_text(encoding="utf-8").strip()

    rows: List[Dict[str, object]] = []
    for dist_key, (matrix_fname, pipeline_prefix) in DIST_MATRICES.items():
        matrix_path = rep_dir / matrix_fname
        if not matrix_path.exists():
            continue

        taxa_names, dist_matrix = _read_phylip_matrix(matrix_path)
        dist_matrix = 0.5 * (dist_matrix + dist_matrix.T)
        np.fill_diagonal(dist_matrix, 0.0)

        for method in ADDITIVITY_METHODS:
            tree_fname = f"{dist_key}_{method.lower()}.nwk"
            tree_path = rep_dir / tree_fname

            if tree_path.exists() and not overwrite:
                est_nwk = tree_path.read_text(encoding="utf-8").strip()
            else:
                est_nwk = infer_tree_from_matrix(
                    method=method,
                    dist_matrix=dist_matrix,
                    taxa_names=taxa_names,
                    max_iter=max_iter,
                )
                if save_trees:
                    tree_path.write_text(est_nwk + "\n", encoding="utf-8")

            effective_lba_pair = lba_pair
            if effective_lba_pair is None and ("T1" in taxa_names and "T20" in taxa_names):
                effective_lba_pair = ("T1", "T20")

            metrics = evaluate_tree_lba(
                est_tree_newick=est_nwk,
                true_tree_newick=true_nwk,
                taxa_names=taxa_names,
                lba_pair=effective_lba_pair,
                auto_max_pair=auto_max_pair,
            )

            pipeline_name = f"{pipeline_prefix}{method}"
            n_taxa = len(taxa_names)
            total_internal = max(n_taxa - 3, 1)
            rf_val = float(metrics["rf"])
            nrf_val = float(metrics["nrf"])

            rows.append(
                {
                    "condition": condition,
                    "distance": dist_val,
                    "length": length_val,
                    "replicate": rep_val,
                    "distance_matrix": dist_key,
                    "method": method,
                    "pipeline": pipeline_name,
                    "num_taxa": n_taxa,
                    "rf_distance": int(round(rf_val)),
                    "nrf_distance": nrf_val,
                    "rcrb": 1.0 - nrf_val,
                    "true_branches": total_internal,
                    "est_branches": total_internal,
                    "ds_est": metrics["ds_est"],
                    "ds_true": metrics["ds_true"],
                    "delta_ds": metrics["delta_ds"],
                    "is_sister": metrics["is_sister"],
                    "nta": metrics["nta"],
                    "pair_nrf": metrics["pair_nrf"],
                    "pair_split_recall": metrics["pair_split_recall"],
                    "pair_split_fnr": metrics["pair_split_fnr"],
                    "lba_quartet_attraction": metrics["lba_quartet_attraction"],
                    "lba_quartet_error": metrics["lba_quartet_error"],
                }
            )

    return rows


def merge_with_benchmark_summary(
    results_dir: Path,
    additivity_df: pd.DataFrame,
    condition: str,
) -> Optional[Path]:
    """
    Merge `additivity_methods_summary.csv` with existing `benchmark_summary.csv`
    (which contains Standard NJ, MSA+ML, MSA+ML_WAG_G4, GS) into
    `combined_benchmark_summary.csv`.

    Parameters
    ----------
    results_dir : Path
        Nextflow results directory containing `benchmark_summary.csv`.
    additivity_df : pd.DataFrame
        DataFrame returned by additivity evaluation.
    condition : str
        Condition label attached to rows (unitless).

    Returns
    -------
    Optional[Path]
        Path to `combined_benchmark_summary.csv`, or None if `benchmark_summary.csv`
        does not exist.
    """
    bench_path = results_dir / "benchmark_summary.csv"
    if not bench_path.exists():
        return None

    bench_df = pd.read_csv(bench_path)
    if "condition" not in bench_df.columns:
        bench_df["condition"] = condition

    pipeline_meta: Dict[str, Tuple[str, str]] = {
        "MMSEQS+NJ": ("sss_linear", "Standard_NJ"),
        "MMSEQS_LOG+NJ": ("sss_log", "Standard_NJ"),
        "PSA+NJ": ("poisson", "Standard_NJ"),
        "MSA+ML": ("msa", "MSA_ML_MFP"),
        "MSA+ML_WAG_G4": ("msa", "MSA_ML_WAG_G4"),
        "GS": ("gs", "GS"),
    }
    bench_df["distance_matrix"] = bench_df["pipeline"].map(
        lambda p: pipeline_meta.get(str(p), ("other", str(p)))[0]
    )
    bench_df["method"] = bench_df["pipeline"].map(
        lambda p: pipeline_meta.get(str(p), ("other", str(p)))[1]
    )

    combined = pd.concat([bench_df, additivity_df], ignore_index=True, sort=False)
    out_path = results_dir / "combined_benchmark_summary.csv"
    combined.to_csv(out_path, index=False)
    return out_path


def evaluate_results_dir(
    results_dir: Path,
    condition: str,
    workers: int = 4,
    max_iter: int = 20,
    lba_pair: Optional[Tuple[str, str]] = None,
    auto_max_pair: bool = False,
    save_trees: bool = True,
    overwrite: bool = False,
    limit: Optional[int] = None,
    output_csv: Optional[Path] = None,
    merge_benchmark: bool = True,
) -> pd.DataFrame:
    """
    Run additivity-based distance methods across all replication directories in
    `results_dir` and save summary CSV files.

    Parameters
    ----------
    results_dir : Path
        Root output directory of a Nextflow run (containing `replications/`).
    condition : str
        Condition label for summary rows (e.g. `v0.3_lba1`) (unitless).
    workers : int, optional
        Number of parallel worker processes (unitless). Default is 4.
    max_iter : int, optional
        Max iterations for `BioNJ_Add_Hybrid2` hill-climb search (unitless).
    lba_pair : Optional[Tuple[str, str]], optional
        Explicit LBA taxon pair `(u, v)` (unitless).
    auto_max_pair : bool, optional
        If True and `lba_pair` is None, pick max patristic distance pair (unitless).
    save_trees : bool, optional
        Whether to save inferred Newick files into each replication dir (unitless).
    overwrite : bool, optional
        Whether to overwrite existing inferred Newick files (unitless).
    limit : Optional[int], optional
        Optional cap on number of replication directories to process (for testing).
    output_csv : Optional[Path], optional
        Destination CSV path. Defaults to `<results_dir>/additivity_methods_summary.csv`.
    merge_benchmark : bool, optional
        Whether to merge with `benchmark_summary.csv` if present (unitless).

    Returns
    -------
    pd.DataFrame
        DataFrame of all evaluated `(distance, length, replicate, distance_matrix, method)` rows.
    """
    rep_root = results_dir / "replications"
    if not rep_root.exists():
        raise FileNotFoundError(f"Replications directory not found: {rep_root}")

    rep_dirs = sorted(
        [d for d in rep_root.iterdir() if d.is_dir() and _DIR_PATTERN.match(d.name)]
    )
    if limit is not None:
        rep_dirs = rep_dirs[:limit]

    tasks = [
        (d, condition, max_iter, lba_pair, auto_max_pair, save_trees, overwrite)
        for d in rep_dirs
    ]

    all_rows: List[Dict[str, object]] = []
    if workers <= 1 or len(tasks) <= 1:
        for task in tasks:
            all_rows.extend(_process_one_replication(task))
    else:
        with ProcessPoolExecutor(max_workers=workers) as ex:
            for chunk_rows in ex.map(_process_one_replication, tasks):
                all_rows.extend(chunk_rows)

    df = pd.DataFrame(all_rows)
    if not df.empty:
        df = df.sort_values(
            ["distance", "length", "replicate", "distance_matrix", "method"]
        ).reset_index(drop=True)

    dest_csv = output_csv or (results_dir / "additivity_methods_summary.csv")
    df.to_csv(dest_csv, index=False)

    if merge_benchmark and not df.empty:
        merge_with_benchmark_summary(results_dir, df, condition)

    return df


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Evaluate BioNJ, FastME (NNI/SPR), AddNJ2/3/4, and BioNJ+Add-Hybrid2 on PHYLIP matrices."
    )
    parser.add_argument("--results_dir", type=Path, required=True, help="Nextflow results directory")
    parser.add_argument("--condition", type=str, required=True, help="Condition label (e.g. v0.3_lba1)")
    parser.add_argument("--workers", type=int, default=4, help="Parallel worker count")
    parser.add_argument("--max_iter", type=int, default=20, help="Max iterations for BioNJ_Add_Hybrid2")
    parser.add_argument("--lba_pair", nargs=2, default=None, help="Explicit LBA taxon pair (e.g. T1 T20)")
    parser.add_argument("--auto_max_pair", action="store_true", help="Auto-select max patristic pair for LBA")
    parser.add_argument("--no_save_trees", action="store_true", help="Do not write .nwk files to replications/")
    parser.add_argument("--overwrite", action="store_true", help="Re-infer even if .nwk already exists")
    parser.add_argument("--limit", type=int, default=None, help="Process at most N replications (testing)")
    parser.add_argument("--output_csv", type=Path, default=None, help="Override output CSV path")
    parser.add_argument("--no_merge", action="store_true", help="Skip merging with benchmark_summary.csv")
    return parser.parse_args()


def main() -> None:
    args = _parse_args()
    lba_pair: Optional[Tuple[str, str]] = tuple(args.lba_pair) if args.lba_pair else None
    df = evaluate_results_dir(
        results_dir=args.results_dir,
        condition=args.condition,
        workers=args.workers,
        max_iter=args.max_iter,
        lba_pair=lba_pair,
        auto_max_pair=args.auto_max_pair,
        save_trees=not args.no_save_trees,
        overwrite=args.overwrite,
        limit=args.limit,
        output_csv=args.output_csv,
        merge_benchmark=not args.no_merge,
    )
    print(f"[eval_additivity] Wrote {len(df)} rows for condition={args.condition}")


if __name__ == "__main__":
    main()
