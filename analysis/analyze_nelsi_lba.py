#!/usr/bin/env python3
"""
Analyze the effect of Long-Branch Attraction (LBA) under the NELSI Yule model
(taxa=20) by comparing lba_ratio=1.0 (baseline), 3.0, and 5.0 across methods.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


DEFAULT_DATASETS: List[Tuple[str, float, Path]] = [
    ("LBA 1x (Baseline)", 1.0, Path("results/results_nelsi_yule/benchmark_summary.csv")),
    ("LBA 3x (T1 & T20)", 3.0, Path("results/results_nelsi_lba3/benchmark_summary.csv")),
    ("LBA 5x (T1 & T20)", 5.0, Path("results/results_nelsi_lba5/benchmark_summary.csv")),
]


def load_lba_datasets(datasets: List[Tuple[str, float, Path]]) -> pd.DataFrame:
    """
    Load and combine benchmark_summary.csv files across LBA ratio conditions.

    Parameters
    ----------
    datasets : List[Tuple[str, float, Path]]
        List of (label, lba_ratio [unitless], csv_path) tuples.

    Returns
    -------
    pd.DataFrame
        Combined dataframe with normalized pipeline names and lba_ratio column.

    Raises
    ------
    FileNotFoundError
        If none of the specified benchmark CSV files exist.
    """
    frames: List[pd.DataFrame] = []
    for label, ratio, path in datasets:
        if not path.is_file():
            print(f"[SKIP] Missing dataset: {path}")
            continue
        df = pd.read_csv(path)
        df["pipeline_disp"] = df["pipeline"].astype(str).str.replace("PWA", "PSA", regex=False)
        df["lba_label"] = label
        df["lba_ratio"] = float(ratio)
        frames.append(df)

    if not frames:
        raise FileNotFoundError("No LBA benchmark_summary.csv files found.")
    return pd.concat(frames, ignore_index=True)


def summarize_lba_penalty(df: pd.DataFrame) -> pd.DataFrame:
    """
    Compute mean nRF (unitless) and Delta-nRF relative to baseline (lba_ratio=1.0)
    for each pipeline, sequence length (aa), and evolutionary distance (substitutions/site).

    Parameters
    ----------
    df : pd.DataFrame
        Combined dataframe returned by `load_lba_datasets`.

    Returns
    -------
    pd.DataFrame
        Summary table with mean_nrf, sd_nrf, and delta_nrf_vs_1x.
    """
    agg = (
        df.groupby(["lba_ratio", "lba_label", "length", "distance", "pipeline_disp"], as_index=False)
        .agg(
            mean_nrf=("nrf_distance", "mean"),
            sd_nrf=("nrf_distance", "std"),
            n_reps=("nrf_distance", "count"),
        )
    )
    base = agg[agg["lba_ratio"] == 1.0][["length", "distance", "pipeline_disp", "mean_nrf"]].rename(
        columns={"mean_nrf": "baseline_nrf"}
    )
    merged = agg.merge(base, on=["length", "distance", "pipeline_disp"], how="left")
    merged["delta_nrf_vs_1x"] = merged["mean_nrf"] - merged["baseline_nrf"]
    return merged


def plot_lba_comparison(df: pd.DataFrame, summary: pd.DataFrame, outdir: Path) -> None:
    """
    Generate publication-quality figures showing nRF vs Distance across LBA ratios
    and Delta-nRF LBA penalty per method.

    Parameters
    ----------
    df : pd.DataFrame
        Combined raw replicate dataframe.
    summary : pd.DataFrame
        Aggregated summary dataframe from `summarize_lba_penalty`.
    outdir : Path
        Output directory for saving PNG plots and summary CSV.

    Returns
    -------
    None
    """
    outdir.mkdir(parents=True, exist_ok=True)
    summary.to_csv(outdir / "nelsi_lba_summary.csv", index=False)

    plt.rcParams.update({"font.size": 10, "axes.titlesize": 11, "axes.labelsize": 10})
    methods = [
        m for m in [
            "PSA+NJ", "PSA+NJ+NNI", "MMSEQS+NJ", "MMSEQS+NJ+NNI",
            "MSA+NJ", "MSA+NJ+NNI", "MSA+ML", "MSA+RAXML", "GS",
        ]
        if m in df["pipeline_disp"].unique()
    ]
    lba_labels = [lbl for lbl, _, _ in DEFAULT_DATASETS if lbl in df["lba_label"].unique()]
    palette: Dict[str, str] = {
        "LBA 1x (Baseline)": "#1f77b4",
        "LBA 3x (T1 & T20)": "#ff7f0e",
        "LBA 5x (T1 & T20)": "#d62728",
    }

    # 1. nRF vs Distance grid (rows: methods, cols: sequence lengths)
    lengths = sorted(df["length"].unique())
    n_rows = len(methods)
    n_cols = len(lengths)
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(5.2 * n_cols, 3.0 * n_rows), sharex=True, sharey=True)
    if n_rows == 1:
        axes = np.expand_dims(axes, axis=0)
    if n_cols == 1:
        axes = np.expand_dims(axes, axis=1)

    for r_idx, method in enumerate(methods):
        for c_idx, length in enumerate(lengths):
            ax = axes[r_idx, c_idx]
            sub = df[(df["pipeline_disp"] == method) & (df["length"] == length)]
            if not sub.empty:
                sns.lineplot(
                    data=sub,
                    x="distance",
                    y="nrf_distance",
                    hue="lba_label",
                    hue_order=lba_labels,
                    palette=palette,
                    marker="o",
                    errorbar=("ci", 95),
                    ax=ax,
                    legend=(r_idx == 0 and c_idx == n_cols - 1),
                )
            if r_idx == 0:
                ax.set_title(f"Length = {length} aa", fontweight="bold")
            if c_idx == 0:
                ax.set_ylabel(f"{method}\nnRF", fontweight="bold")
            else:
                ax.set_ylabel("")
            if r_idx == n_rows - 1:
                ax.set_xlabel("Distance D (substitutions/site)", fontweight="bold")
            ax.set_ylim(-0.02, 0.70)
            ax.grid(True, linestyle="--", alpha=0.4)

    fig.suptitle("NELSI Yule (taxa=20): Impact of LBA (1x vs 3x vs 5x) on nRF", fontweight="bold", fontsize=13, y=1.005)
    plt.tight_layout()
    fig.savefig(outdir / "nelsi_lba_nrf_curves.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    # 2. Mean Delta-nRF penalty barplot across methods (LBA > 1.0)
    pen = summary[summary["lba_ratio"] > 1.0].dropna(subset=["delta_nrf_vs_1x"])
    if not pen.empty:
        fig2, ax2 = plt.subplots(figsize=(11, 5.5))
        sns.barplot(
            data=pen,
            x="pipeline_disp",
            y="delta_nrf_vs_1x",
            hue="lba_label",
            order=methods,
            palette={k: v for k, v in palette.items() if k in pen["lba_label"].unique()},
            errorbar=("ci", 95),
            ax=ax2,
        )
        ax2.axhline(0.0, color="black", linewidth=0.9, linestyle="--")
        ax2.set_title("Mean LBA Topological Penalty (Delta nRF = nRF_LBA - nRF_1x)", fontweight="bold")
        ax2.set_xlabel("Inference Pipeline", fontweight="bold")
        ax2.set_ylabel("Delta nRF (positive = worse under LBA)", fontweight="bold")
        ax2.tick_params(axis="x", rotation=25)
        ax2.grid(True, axis="y", linestyle="--", alpha=0.4)
        plt.tight_layout()
        fig2.savefig(outdir / "nelsi_lba_delta_nrf_bar.png", dpi=300, bbox_inches="tight")
        plt.close(fig2)


def parse_args() -> argparse.Namespace:
    """
    Parse command-line arguments.

    Returns
    -------
    argparse.Namespace
        Parsed CLI arguments.
    """
    parser = argparse.ArgumentParser(description="Analyze NELSI LBA (taxa=20) benchmark results.")
    parser.add_argument(
        "--outdir",
        type=Path,
        default=Path("results/results_nelsi_lba_comparison"),
        help="Output directory for LBA comparison plots and CSV.",
    )
    return parser.parse_args()


def main() -> None:
    """
    Main entry point for NELSI LBA comparison analysis.

    Returns
    -------
    None
    """
    args = parse_args()
    df = load_lba_datasets(DEFAULT_DATASETS)
    summary = summarize_lba_penalty(df)
    plot_lba_comparison(df, summary, args.outdir)
    print(f"Saved LBA analysis outputs to {args.outdir}")


if __name__ == "__main__":
    main()
