#!/usr/bin/env python3
"""
analyze_within_hog_tcs.py
=========================
Performs rigorous within-HOG comparisons among the 4 phylogenetic methods
(MSA+ML, MSA+NJ, PSA+NJ, GS) using MAD-rooted TCS scores across 400 real biological HOGs.

Evaluates:
  1. Overall & per-group win rates (HOGs where method achieved max TCS).
  2. Average rank within each HOG (1.0 = best, 4.0 = worst).
  3. Normalized TCS ratio within HOG: TCS / max(TCS in HOG).
  4. Pairwise Head-to-Head win / tie / loss records.
  5. Stratified statistics across SSS intervals (<0.2, 0.2-0.4, 0.4-0.6, 0.6-0.8, 0.8-1.0).

Outputs:
  - dataset/results/within_hog_win_rates.csv
  - dataset/results/sss_binned_tcs_stats.csv
  - dataset/results/plot_sss_binned_win_rates.png
  - dataset/results/plot_within_hog_normalized_tcs.png
"""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent.parent
DATASET_DIR = SCRIPT_DIR.parent
RESULTS_DIR = DATASET_DIR / "results"

plt.rcParams.update({
    "font.sans-serif": ["Arial", "DejaVu Sans", "Helvetica", "sans-serif"],
    "font.size": 11,
    "axes.labelsize": 12,
    "axes.titlesize": 13,
    "xtick.labelsize": 10.5,
    "ytick.labelsize": 10.5,
    "legend.fontsize": 10.5,
    "figure.titlesize": 15,
    "axes.grid": True,
    "grid.alpha": 0.3,
    "grid.linestyle": "--"
})

METHODS = [
    ("msa_ml", "MSA+ML (IQ-TREE 2)", "#2980b9", "o"),
    ("msa_nj", "MSA+NJ (RapidNJ)",   "#e67e22", "s"),
    ("psa_nj", "PSA+NJ (RapidNJ)",   "#27ae60", "^"),
    ("gs",     "GS (Graph Splitting)","#c0392b", "D")
]

METHOD_KEYS = [m[0] for m in METHODS]
METHOD_NAMES = {m[0]: m[1] for m in METHODS}
METHOD_COLORS = {m[0]: m[2] for m in METHODS}
METHOD_MARKERS = {m[0]: m[3] for m in METHODS}

GROUPS = [
    ("amn", "Amniota (10 taxa)"),
    ("euk", "Eukaryota (50 taxa)"),
    ("act", "Actinobacteria (50 taxa)"),
    ("bac", "Bacteria (50 taxa)")
]


def load_data() -> pd.DataFrame:
    csv_path = RESULTS_DIR / "dataset_tcs_summary.csv"
    if not csv_path.exists():
        raise FileNotFoundError(f"{csv_path} not found.")
    df = pd.read_csv(csv_path)

    # 1. Compute within-HOG ranks and normalized scores
    norm_scores = {m: [] for m in METHOD_KEYS}
    ranks = {m: [] for m in METHOD_KEYS}
    fractional_wins = {m: [] for m in METHOD_KEYS}

    for idx, row in df.iterrows():
        s_series = pd.Series({m: row[f"tcs_{m}"] for m in METHOD_KEYS})
        # rank descending: highest score = 1
        r = s_series.rank(ascending=False, method="average")
        max_val = s_series.max()

        # winners (handles ties)
        top_methods = s_series[s_series >= max_val - 1e-5].index.tolist()
        for m in METHOD_KEYS:
            ranks[m].append(r[m])
            norm_scores[m].append(row[f"tcs_{m}"] / max_val if max_val > 0 else 1.0)
            fractional_wins[m].append(1.0 / len(top_methods) if m in top_methods else 0.0)

    for m in METHOD_KEYS:
        df[f"rank_{m}"] = ranks[m]
        df[f"norm_{m}"] = norm_scores[m]
        df[f"win_{m}"] = fractional_wins[m]

    # SSS bins
    bins = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
    labels = ["< 0.2\n(超遠縁)", "0.2 - 0.4\n(遠縁)", "0.4 - 0.6\n(中等度)", "0.6 - 0.8\n(近縁)", "0.8 - 1.0\n(極近縁)"]
    labels_simple = ["< 0.2", "0.2 - 0.4", "0.4 - 0.6", "0.6 - 0.8", "0.8 - 1.0"]
    df["sss_bin"] = pd.cut(df["mean_sss"], bins=bins, labels=labels, include_lowest=True)
    df["sss_bin_simple"] = pd.cut(df["mean_sss"], bins=bins, labels=labels_simple, include_lowest=True)

    return df


def compute_pairwise_head_to_head(df: pd.DataFrame) -> pd.DataFrame:
    """Computes pairwise win/tie/loss across all HOGs."""
    pairs = [
        ("msa_ml", "gs"),
        ("msa_ml", "msa_nj"),
        ("msa_ml", "psa_nj"),
        ("psa_nj", "msa_nj"),
        ("gs", "msa_nj"),
        ("gs", "psa_nj")
    ]
    records = []
    for m1, m2 in pairs:
        diff = df[f"tcs_{m1}"] - df[f"tcs_{m2}"]
        w1 = (diff > 1e-4).sum()
        tie = (diff.abs() <= 1e-4).sum()
        w2 = (diff < -1e-4).sum()
        records.append({
            "Method 1": METHOD_NAMES[m1],
            "Method 2": METHOD_NAMES[m2],
            "M1 Wins": w1,
            "M1 Win %": round(w1 / len(df) * 100, 1),
            "Ties": tie,
            "Tie %": round(tie / len(df) * 100, 1),
            "M2 Wins": w2,
            "M2 Win %": round(w2 / len(df) * 100, 1),
            "Mean ΔTCS (M1-M2)": round(diff.mean(), 3)
        })
    return pd.DataFrame(records)


def compute_group_summary(df: pd.DataFrame) -> pd.DataFrame:
    """Computes win rates, mean rank, and mean normalized score per group."""
    records = []
    for g_key, g_title in GROUPS + [("all", "All Biological Groups (400 HOGs)")]:
        sub = df if g_key == "all" else df[df["dataset"] == g_key]
        n_hogs = len(sub)
        for m in METHOD_KEYS:
            wins = sub[f"win_{m}"].sum()
            records.append({
                "Group": g_key,
                "Group_Name": g_title,
                "Method_Key": m,
                "Method": METHOD_NAMES[m],
                "HOG_Count": n_hogs,
                "Wins": round(wins, 1),
                "Win_Rate_%": round(wins / n_hogs * 100, 1),
                "Mean_Rank": round(sub[f"rank_{m}"].mean(), 2),
                "Mean_Norm_TCS": round(sub[f"norm_{m}"].mean(), 4)
            })
    return pd.DataFrame(records)


def compute_sss_binned_summary(df: pd.DataFrame) -> pd.DataFrame:
    """Computes win rates, mean rank, and mean normalized score per SSS bin."""
    records = []
    for b in ["< 0.2", "0.2 - 0.4", "0.4 - 0.6", "0.6 - 0.8", "0.8 - 1.0"]:
        sub = df[df["sss_bin_simple"] == b]
        n_hogs = len(sub)
        if n_hogs == 0:
            continue
        for m in METHOD_KEYS:
            wins = sub[f"win_{m}"].sum()
            records.append({
                "SSS_Bin": b,
                "HOG_Count": n_hogs,
                "Method_Key": m,
                "Method": METHOD_NAMES[m],
                "Wins": round(wins, 1),
                "Win_Rate_%": round(wins / n_hogs * 100, 1),
                "Mean_Rank": round(sub[f"rank_{m}"].mean(), 2),
                "Mean_Norm_TCS": round(sub[f"norm_{m}"].mean(), 4)
            })
    return pd.DataFrame(records)


def plot_win_rates_by_sss(df_sss: pd.DataFrame, out_path: Path):
    """Plots grouped bar chart of win rates across SSS intervals."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))

    bins_order = ["< 0.2", "0.2 - 0.4", "0.4 - 0.6", "0.6 - 0.8", "0.8 - 1.0"]
    x = np.arange(len(bins_order))
    width = 0.2

    # Left: Win Rate %
    for i, (m_key, m_name, color, _) in enumerate(METHODS):
        sub_m = df_sss[df_sss["Method_Key"] == m_key].set_index("SSS_Bin").reindex(bins_order)
        bars = ax1.bar(x + (i - 1.5) * width, sub_m["Win_Rate_%"], width, label=m_name, color=color, alpha=0.88, edgecolor="black", linewidth=0.5)
        for bar in bars:
            h = bar.get_height()
            if h > 5:
                ax1.annotate(f"{h:.1f}%",
                             xy=(bar.get_x() + bar.get_width() / 2, h),
                             xytext=(0, 3), textcoords="offset points",
                             ha="center", va="bottom", fontsize=8.5, fontweight="bold")

    ax1.set_title("(A) Win Rate (% HOGs where Method is Best) by SSS Interval", fontweight="bold", fontsize=12.5)
    ax1.set_xlabel("Sequence Similarity Score (SSS: $\\bar{w}$)")
    ax1.set_ylabel("Win Rate (%) within Same HOG")
    ax1.set_xticks(x)
    ax1.set_xticklabels(["< 0.2\n(N=31)", "0.2 - 0.4\n(N=97)", "0.4 - 0.6\n(N=138)", "0.6 - 0.8\n(N=96)", "0.8 - 1.0\n(N=38)"])
    ax1.set_ylim(0, 52)
    ax1.legend(frameon=True, loc="upper right")

    # Right: Mean Relative TCS (TCS / Max_TCS in HOG)
    for i, (m_key, m_name, color, marker) in enumerate(METHODS):
        sub_m = df_sss[df_sss["Method_Key"] == m_key].set_index("SSS_Bin").reindex(bins_order)
        ax2.plot(x, sub_m["Mean_Norm_TCS"], marker=marker, label=m_name, color=color, linewidth=2.4, markersize=8)
        for idx_pt, val in enumerate(sub_m["Mean_Norm_TCS"]):
            ax2.annotate(f"{val:.3f}", xy=(x[idx_pt], val), xytext=(0, 6 if i % 2 == 0 else -14),
                         textcoords="offset points", ha="center", fontsize=8.5, color=color, fontweight="bold")

    ax2.set_title("(B) Mean Relative TCS (Normalized to Best=1.0) by SSS Interval", fontweight="bold", fontsize=12.5)
    ax2.set_xlabel("Sequence Similarity Score (SSS: $\\bar{w}$)")
    ax2.set_ylabel("Relative TCS (TCS / Max_TCS in HOG)")
    ax2.set_xticks(x)
    ax2.set_xticklabels(["< 0.2\n(N=31)", "0.2 - 0.4\n(N=97)", "0.4 - 0.6\n(N=138)", "0.6 - 0.8\n(N=96)", "0.8 - 1.0\n(N=38)"])
    ax2.set_ylim(0.65, 1.02)
    ax2.legend(frameon=True, loc="lower right")

    plt.suptitle("Within-HOG Performance of Phylogenetic Methods Stratified by Sequence Similarity (SSS)", fontweight="bold", y=0.98)
    plt.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out_path}")


def plot_within_hog_normalized_scatter(df: pd.DataFrame, out_path: Path):
    """
    Plots Normalized TCS (TCS / Max_TCS in HOG) vs SSS across 4 biological groups (2x2).
    This eliminates tree-size confounding across different HOGs!
    """
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    axes = axes.flatten()

    for idx, (d_key, title_name) in enumerate(GROUPS):
        ax = axes[idx]
        sub = df[df["dataset"] == d_key]

        for m_key, m_label, color, marker in METHODS:
            norm_col = f"norm_{m_key}"
            valid = sub.dropna(subset=["mean_sss", norm_col])

            ax.scatter(
                valid["mean_sss"],
                valid[norm_col],
                label=m_label,
                color=color,
                marker=marker,
                alpha=0.6,
                s=45,
                edgecolors="black",
                linewidth=0.3
            )

            # Trend line
            if len(valid) > 2:
                sns.regplot(
                    data=valid,
                    x="mean_sss",
                    y=norm_col,
                    scatter=False,
                    ax=ax,
                    color=color,
                    line_kws={"linewidth": 2.2, "alpha": 0.85},
                    ci=None
                )

        ax.set_title(title_name, fontweight="bold", fontsize=13)
        ax.set_xlabel("Sequence Similarity Score (SSS: $\\bar{w}$)")
        ax.set_ylabel("Relative TCS (TCS / HOG Max Score: 0 to 1)")
        ax.set_xlim(-0.02, 1.02)
        ax.set_ylim(-0.02, 1.08)
        ax.legend(frameon=True, loc="lower left", framealpha=0.9)

    plt.suptitle("Within-HOG Normalized Taxonomic Congruence (Relative to Best Method) vs. SSS", fontweight="bold", y=0.98)
    plt.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out_path}")


def main():
    print("=" * 80)
    print("  Within-HOG Performance Analysis (Controlling for Tree Size & SSS)")
    print("=" * 80)

    df = load_data()

    # 1. Group summary
    df_group = compute_group_summary(df)
    group_csv = RESULTS_DIR / "within_hog_win_rates.csv"
    df_group.to_csv(group_csv, index=False)
    print(f"\nSaved Group Summary: {group_csv}")

    # 2. Pairwise head-to-head
    df_h2h = compute_pairwise_head_to_head(df)
    h2h_csv = RESULTS_DIR / "pairwise_head_to_head.csv"
    df_h2h.to_csv(h2h_csv, index=False)
    print(f"Saved Pairwise Head-to-Head: {h2h_csv}")

    # 3. SSS binned summary
    df_sss = compute_sss_binned_summary(df)
    sss_csv = RESULTS_DIR / "sss_binned_tcs_stats.csv"
    df_sss.to_csv(sss_csv, index=False)
    print(f"Saved SSS Binned Summary: {sss_csv}")

    # Display tables in console
    print("\n=== (1) OVERALL & PER-GROUP WITHIN-HOG WIN RATES & RANKS ===")
    pivot_wins = df_group.pivot(index="Group_Name", columns="Method", values="Win_Rate_%")
    pivot_rank = df_group.pivot(index="Group_Name", columns="Method", values="Mean_Rank")
    print("\n[Win Rate (% of HOGs where method is Best)]")
    print(pivot_wins.to_string())
    print("\n[Mean Rank (1.0 = Best, 4.0 = Worst)]")
    print(pivot_rank.to_string())

    print("\n=== (2) PAIRWISE HEAD-TO-HEAD WITHIN SAME HOG (N=400 HOGs) ===")
    print(df_h2h.to_string(index=False))

    print("\n=== (3) WITHIN-HOG STATISTICS BY SSS INTERVAL ===")
    pivot_sss_wins = df_sss.pivot(index="SSS_Bin", columns="Method", values="Win_Rate_%")
    pivot_sss_norm = df_sss.pivot(index="SSS_Bin", columns="Method", values="Mean_Norm_TCS")
    print("\n[Win Rate (%) across SSS Intervals]")
    print(pivot_sss_wins.to_string())
    print("\n[Mean Relative Score (TCS / HOG Max) across SSS Intervals]")
    print(pivot_sss_norm.to_string())

    # 4. Plots
    print("\nGenerating Plots...")
    out_bars = RESULTS_DIR / "plot_sss_binned_win_rates.png"
    plot_win_rates_by_sss(df_sss, out_bars)

    out_norm_scatter = RESULTS_DIR / "plot_within_hog_normalized_tcs.png"
    plot_within_hog_normalized_scatter(df, out_norm_scatter)

    print("\nAll within-HOG analyses and plots completed successfully!")


if __name__ == "__main__":
    main()
