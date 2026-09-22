#!/usr/bin/env python3
"""
Model Selection Impact Analysis: True Model (WAG+G4) vs Mis-selected Models in MSA+ML.

This script compares topological accuracy (normalized Robinson-Foulds distance: nRF)
between replicates where ModelFinder correctly selected the true generative model
(WAG+G4) versus replicates where it selected alternative models (WAG+I+G4, WAG+R3,
VT, BLOSUM62, PMB) under Experiment 16 (INDELible v1.03 benchmark).

Units:
------
- Evolution Distance (distance): substitutions/site (unitless rate scale)
- Normalized Robinson-Foulds distance (nRF): unitless in [0.0, 1.0]
- Sequence Similarity Score (SSS, mean_w): unitless normalized bitscore in [0.0, 1.0]
- Sequence length: amino acid residues (aa)
"""

from typing import Dict, List, Tuple, Any, Optional
import os
import sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns

PROJECT_ROOT: Path = Path(__file__).resolve().parent.parent.parent
CSV_PATH: Path = PROJECT_ROOT / "results" / "results_paper_tree2" / "benchmark_summary.csv"
OUT_DIR: Path = PROJECT_ROOT / "results" / "results_paper_tree2"
OUT_PNG: Path = OUT_DIR / "model_selection_nrf_comparison.png"
OUT_CSV: Path = OUT_DIR / "model_selection_comparison_stats.csv"

# Global plotting aesthetics
plt.rcParams.update({
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "font.family": "sans-serif",
    "axes.edgecolor": "#333333",
    "axes.linewidth": 1.1,
    "grid.color": "#e5e5e5",
    "grid.linestyle": "--",
    "grid.alpha": 0.7,
    "savefig.dpi": 300,
    "figure.dpi": 300,
})


def categorize_model(model_name: str) -> str:
    """
    Categorizes the BIC-selected substitution model into distinct analytical classes.

    Parameters
    ----------
    model_name : str
        The raw model string selected by ModelFinder (e.g., 'WAG+G4', 'WAG+I+G4').

    Returns
    -------
    str
        Simplified category label:
        - 'WAG+G4 (True Model)'
        - 'WAG+I+G4 (Invariable Sites)'
        - 'WAG+R (FreeRate)'
        - 'Alternative Matrices (VT/BLOSUM/PMB)'
    """
    if not isinstance(model_name, str):
        return "Unknown"
    m = model_name.strip()
    if m == "WAG+G4":
        return "WAG+G4 (True)"
    elif "WAG" in m and "+I" in m:
        return "WAG+I+G4 (+I)"
    elif "WAG" in m and ("+R" in m or "+R2" in m or "+R3" in m):
        return "WAG+R (FreeRate)"
    elif "WAG" in m:
        return "WAG Other"
    elif any(mat in m for mat in ["VT", "BLOSUM", "PMB", "JTT", "LG", "DAYHOFF"]):
        return "Other Matrix (VT/BLOSUM/PMB)"
    return "Other"


def load_and_preprocess_data(csv_path: Path) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Loads benchmark results and preprocesses MSA+ML and MSA+RAXML subsets.

    Parameters
    ----------
    csv_path : Path
        Path to benchmark_summary.csv.

    Returns
    -------
    Tuple[pd.DataFrame, pd.DataFrame]
        Tuple of (df_ml, df_raxml) with added categorization columns.

    Raises
    ------
    FileNotFoundError
        If benchmark_summary.csv does not exist.
    """
    if not csv_path.exists():
        raise FileNotFoundError(f"Benchmark summary file not found: {csv_path}")

    df = pd.read_csv(csv_path)

    # Filter MSA+ML
    df_ml = df[df["pipeline"] == "MSA+ML"].copy()
    df_ml["is_true_model"] = df_ml["best_model_bic"] == "WAG+G4"
    df_ml["model_category"] = df_ml["best_model_bic"].apply(categorize_model)
    df_ml["model_group"] = np.where(df_ml["is_true_model"], "True Model (WAG+G4)", "Other Models (BIC Mis-selected)")

    # Filter MSA+RAXML for reference
    df_raxml = df[df["pipeline"] == "MSA+RAXML"].copy()

    return df_ml, df_raxml


def compute_group_statistics(df_ml: pd.DataFrame, df_raxml: pd.DataFrame) -> pd.DataFrame:
    """
    Computes summary statistics and statistical significance tests across distances.

    Parameters
    ----------
    df_ml : pd.DataFrame
        Preprocessed DataFrame for MSA+ML pipeline.
    df_raxml : pd.DataFrame
        Preprocessed DataFrame for MSA+RAXML pipeline.

    Returns
    -------
    pd.DataFrame
        Summary statistics table indexed by distance.
    """
    records: List[Dict[str, Any]] = []
    distances = sorted(df_ml["distance"].unique())

    for d in distances:
        sub_ml = df_ml[df_ml["distance"] == d]
        sub_raxml = df_raxml[df_raxml["distance"] == d]

        true_grp = sub_ml[sub_ml["is_true_model"]]["nrf_distance"]
        other_grp = sub_ml[~sub_ml["is_true_model"]]["nrf_distance"]
        raxml_grp = sub_raxml["nrf_distance"]

        n_true = len(true_grp)
        n_other = len(other_grp)
        n_raxml = len(raxml_grp)

        mean_true = true_grp.mean() if n_true > 0 else np.nan
        std_true = true_grp.std() if n_true > 0 else np.nan
        sem_true = true_grp.sem() if n_true > 0 else np.nan

        mean_other = other_grp.mean() if n_other > 0 else np.nan
        std_other = other_grp.std() if n_other > 0 else np.nan
        sem_other = other_grp.sem() if n_other > 0 else np.nan

        mean_raxml = raxml_grp.mean() if n_raxml > 0 else np.nan
        std_raxml = raxml_grp.std() if n_raxml > 0 else np.nan

        p_mwu = np.nan
        if n_true > 0 and n_other > 0:
            try:
                _, p_mwu = stats.mannwhitneyu(true_grp, other_grp, alternative="two-sided")
            except ValueError:
                p_mwu = 1.0

        records.append({
            "distance": d,
            "N_true": n_true,
            "mean_nrf_true": mean_true,
            "std_nrf_true": std_true,
            "sem_nrf_true": sem_true,
            "N_other": n_other,
            "mean_nrf_other": mean_other,
            "std_nrf_other": std_other,
            "sem_nrf_other": sem_other,
            "diff_mean_nrf (other - true)": mean_other - mean_true if not np.isnan(mean_other) and not np.isnan(mean_true) else np.nan,
            "p_value_mwu": p_mwu,
            "mean_nrf_raxml": mean_raxml,
            "std_nrf_raxml": std_raxml,
        })

    stat_df = pd.DataFrame(records)
    return stat_df


def generate_publication_figure(
    df_ml: pd.DataFrame,
    df_raxml: pd.DataFrame,
    stat_df: pd.DataFrame,
    out_png: Path
) -> None:
    """
    Renders a 4-panel publication-grade figure comparing nRF topological error.

    Parameters
    ----------
    df_ml : pd.DataFrame
        MSA+ML data.
    df_raxml : pd.DataFrame
        MSA+RAXML reference data.
    stat_df : pd.DataFrame
        Statistical summary table.
    out_png : Path
        Destination path for saving the PNG image.
    """
    fig = plt.figure(figsize=(15, 12))
    gs = fig.add_gridspec(2, 2, hspace=0.28, wspace=0.22)

    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_d = fig.add_subplot(gs[1, 1])

    # -------------------------------------------------------------------------
    # Panel A: Evolution Distance vs Mean nRF (True vs Other vs RAxML)
    # -------------------------------------------------------------------------
    distances = stat_df["distance"].values

    # True Model
    valid_t = stat_df["N_true"] > 0
    ax_a.errorbar(
        stat_df.loc[valid_t, "distance"],
        stat_df.loc[valid_t, "mean_nrf_true"],
        yerr=stat_df.loc[valid_t, "sem_nrf_true"],
        fmt="o-",
        color="#1f77b4",
        linewidth=2.5,
        markersize=8,
        capsize=4,
        capthick=1.5,
        label="MSA+ML: True Model (WAG+G4)",
        zorder=5
    )

    # Other Models
    valid_o = stat_df["N_other"] > 0
    ax_a.errorbar(
        stat_df.loc[valid_o, "distance"],
        stat_df.loc[valid_o, "mean_nrf_other"],
        yerr=stat_df.loc[valid_o, "sem_nrf_other"],
        fmt="s--",
        color="#d62728",
        linewidth=2.5,
        markersize=8,
        capsize=4,
        capthick=1.5,
        label="MSA+ML: Other Models (BIC Mis-selected)",
        zorder=4
    )

    # RAxML (Fixed PROTGAMMAIWAGX)
    ax_a.plot(
        stat_df["distance"],
        stat_df["mean_nrf_raxml"],
        "^:",
        color="#2ca02c",
        linewidth=2.0,
        markersize=7,
        label="MSA+RAXML: Fixed Model (WAG+G4+I)",
        zorder=3
    )

    # Replicate scatter points with jitter
    rng = np.random.default_rng(42)
    for is_true, c in [(True, "#1f77b4"), (False, "#d62728")]:
        sub = df_ml[df_ml["is_true_model"] == is_true]
        jitter = rng.uniform(-0.03, 0.03, size=len(sub))
        ax_a.scatter(
            sub["distance"] + jitter,
            sub["nrf_distance"],
            color=c,
            alpha=0.18,
            s=22,
            edgecolor="none",
            zorder=2
        )

    # Add statistical significance annotations above points
    for _, row in stat_df.iterrows():
        d_val = row["distance"]
        p_val = row["p_value_mwu"]
        if not np.isnan(p_val) and p_val < 0.05 and d_val >= 1.0:
            y_pos = max(row["mean_nrf_other"], row["mean_nrf_true"]) + 0.04
            star = "***" if p_val < 0.001 else "**" if p_val < 0.01 else "*"
            ax_a.annotate(
                f"{star}\np={p_val:.3f}",
                xy=(d_val, y_pos),
                xytext=(0, 2),
                textcoords="offset points",
                ha="center",
                fontsize=8.5,
                fontweight="bold",
                color="#800000"
            )

    ax_a.set_title("A. Topological Error (nRF) vs Distance D", fontsize=13, fontweight="bold", pad=10)
    ax_a.set_xlabel("Evolutionary Distance D (substitutions/site)", fontsize=11, fontweight="bold")
    ax_a.set_ylabel("Normalized RF Distance (nRF, lower is better)", fontsize=11, fontweight="bold")
    ax_a.set_ylim(-0.02, 0.70)
    ax_a.grid(True)
    ax_a.legend(frameon=True, facecolor="white", edgecolor="#cccccc", fontsize=9.5, loc="upper left")

    # -------------------------------------------------------------------------
    # Panel B: Model Selection Breakdown across Distance D (Stacked Bar)
    # -------------------------------------------------------------------------
    category_order = [
        "WAG+G4 (True)",
        "WAG+I+G4 (+I)",
        "WAG+R (FreeRate)",
        "Other Matrix (VT/BLOSUM/PMB)"
    ]
    palette_cat = {
        "WAG+G4 (True)": "#1f77b4",
        "WAG+I+G4 (+I)": "#ff7f0e",
        "WAG+R (FreeRate)": "#9467bd",
        "Other Matrix (VT/BLOSUM/PMB)": "#d62728"
    }

    # Group counts
    ctab = pd.crosstab(df_ml["distance"], df_ml["model_category"], normalize="index") * 100
    for col in category_order:
        if col not in ctab.columns:
            ctab[col] = 0.0
    ctab = ctab[category_order]

    x_indices = np.arange(len(ctab.index))
    bar_width = 0.55
    bottom = np.zeros(len(ctab.index))

    for col in category_order:
        values = ctab[col].values
        ax_b.bar(
            x_indices,
            values,
            bar_width,
            bottom=bottom,
            label=col,
            color=palette_cat[col],
            edgecolor="white",
            linewidth=1.0
        )
        # Label percentages inside bars if >= 8%
        for idx, (v, b) in enumerate(zip(values, bottom)):
            if v >= 9.0:
                ax_b.text(
                    idx,
                    b + v / 2,
                    f"{v:.0f}%",
                    ha="center",
                    va="center",
                    color="white",
                    fontweight="bold",
                    fontsize=8.5
                )
        bottom += values

    ax_b.set_title("B. ModelFinder Selected Models vs Distance D", fontsize=13, fontweight="bold", pad=10)
    ax_b.set_xlabel("Evolutionary Distance D (substitutions/site)", fontsize=11, fontweight="bold")
    ax_b.set_ylabel("Selection Frequency (%)", fontsize=11, fontweight="bold")
    ax_b.set_xticks(x_indices)
    ax_b.set_xticklabels([f"{d:.1f}" for d in ctab.index])
    ax_b.set_ylim(0, 100)
    ax_b.grid(axis="y")
    ax_b.legend(frameon=True, facecolor="white", edgecolor="#cccccc", fontsize=9.0, loc="lower left")

    # -------------------------------------------------------------------------
    # Panel C: Detailed Model Group Boxplot at High Divergence (D in [1.0, 1.5, 2.0])
    # -------------------------------------------------------------------------
    high_d_df = df_ml[df_ml["distance"].isin([1.0, 1.5, 2.0])].copy()

    # Define grouped boxplot
    sns.boxplot(
        data=high_d_df,
        x="distance",
        y="nrf_distance",
        hue="model_category",
        hue_order=category_order,
        palette=palette_cat,
        ax=ax_c,
        width=0.7,
        fliersize=0,
        boxprops=dict(alpha=0.85, linewidth=1.2)
    )

    # Overlay stripplot
    sns.stripplot(
        data=high_d_df,
        x="distance",
        y="nrf_distance",
        hue="model_category",
        hue_order=category_order,
        palette=palette_cat,
        dodge=True,
        jitter=0.2,
        alpha=0.6,
        size=4.5,
        ax=ax_c,
        legend=False
    )

    ax_c.set_title("C. nRF Error Distribution by Model Class (D ≥ 1.0)", fontsize=13, fontweight="bold", pad=10)
    ax_c.set_xlabel("Evolutionary Distance D (substitutions/site)", fontsize=11, fontweight="bold")
    ax_c.set_ylabel("Normalized RF Distance (nRF)", fontsize=11, fontweight="bold")
    ax_c.set_ylim(-0.02, 0.75)
    ax_c.grid(True)
    ax_c.legend(frameon=True, facecolor="white", edgecolor="#cccccc", fontsize=9.0, loc="upper left")

    # -------------------------------------------------------------------------
    # Panel D: Sequence Similarity Score (SSS, mean_w) vs nRF
    # -------------------------------------------------------------------------
    sns.scatterplot(
        data=df_ml,
        x="sss_mean",
        y="nrf_distance",
        hue="model_group",
        hue_order=["True Model (WAG+G4)", "Other Models (BIC Mis-selected)"],
        palette={"True Model (WAG+G4)": "#1f77b4", "Other Models (BIC Mis-selected)": "#d62728"},
        alpha=0.65,
        s=45,
        edgecolor="#333333",
        linewidth=0.5,
        ax=ax_d
    )

    ax_d.set_title("D. nRF Error vs Sequence Similarity Score (SSS: w)", fontsize=13, fontweight="bold", pad=10)
    ax_d.set_xlabel("Mean Sequence Similarity Score (SSS, w; log scale)", fontsize=11, fontweight="bold")
    ax_d.set_ylabel("Normalized RF Distance (nRF)", fontsize=11, fontweight="bold")
    ax_d.set_xscale("log")
    ax_d.set_ylim(-0.02, 0.70)
    ax_d.grid(True, which="both")
    ax_d.legend(frameon=True, facecolor="white", edgecolor="#cccccc", fontsize=9.0, loc="upper right")

    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[OK] Publication figure saved to: {out_png}")


def main() -> None:
    """Main execution entrypoint."""
    print(f"Loading data from: {CSV_PATH}")
    df_ml, df_raxml = load_and_preprocess_data(CSV_PATH)

    print("Computing statistical summary...")
    stat_df = compute_group_statistics(df_ml, df_raxml)
    stat_df.to_csv(OUT_CSV, index=False)
    print(f"[OK] Statistics saved to: {OUT_CSV}")

    # Print summary table to stdout
    print("\n" + "=" * 95)
    print(f"{'D':<5} | {'True N':<6} {'Mean nRF':<10} {'SEM':<8} | {'Other N':<7} {'Mean nRF':<10} {'SEM':<8} | {'Diff(O-T)':<10} | {'MWU p-value':<12}")
    print("-" * 95)
    for _, r in stat_df.iterrows():
        p_str = f"{r['p_value_mwu']:.4e}" if not np.isnan(r['p_value_mwu']) else "N/A"
        diff_str = f"{r['diff_mean_nrf (other - true)']:+.4f}" if not np.isnan(r['diff_mean_nrf (other - true)']) else "N/A"
        print(f"{r['distance']:<5.1f} | {int(r['N_true']):<6} {r['mean_nrf_true']:<10.4f} {r['sem_nrf_true']:<8.4f} | {int(r['N_other']):<7} {r['mean_nrf_other']:<10.4f} {r['sem_nrf_other']:<8.4f} | {diff_str:<10} | {p_str:<12}")
    print("=" * 95 + "\n")

    print("Generating publication figure...")
    generate_publication_figure(df_ml, df_raxml, stat_df, OUT_PNG)
    print("Execution completed successfully.")


if __name__ == "__main__":
    main()
