#!/usr/bin/env python3
"""
analyze_branch_lengths.py
=========================
Calculates branch length distributions, min/max values, and performs clustering
across phylogenetic tree estimation results (MSA+ML, MSA+NJ, PSA+NJ, and GS)
for the 4 biological groups (Amniota, Eukaryota, Actinobacteria, Bacteria).

Key outputs:
1. Overall summary table of branch lengths (min, max, mean, median, std, percentiles)
2. Biological group breakdown of branch lengths
3. Per-HOG branch length metrics table (CSV)
4. Clustering analysis of branch length distributions (Hierarchical Clustering & PCA)
5. Visualizations:
   - Branch length distributions (violin / box plots)
   - Dendrogram and PCA projection of branch length clusters
"""

import sys
import re
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
import scipy.cluster.hierarchy as sch
from scipy.spatial.distance import pdist, squareform
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import Phylo

SCRIPT_DIR = Path(__file__).resolve().parent
BASE_DIR = SCRIPT_DIR.parent
RESULTS_DIR = BASE_DIR / "results"
TREES_DIR = RESULTS_DIR / "trees"

# Set publication style
plt.rcParams.update({
    "font.sans-serif": ["Arial", "DejaVu Sans", "Helvetica", "sans-serif"],
    "font.size": 11,
    "axes.labelsize": 12,
    "axes.titlesize": 13,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 10,
    "figure.titlesize": 14,
    "lines.linewidth": 1.5,
    "axes.grid": True,
    "grid.alpha": 0.3,
    "grid.linestyle": "--"
})

GROUP_NAMES = {
    "act": "Actinobacteria",
    "amn": "Amniota",
    "bac": "Bacteria",
    "euk": "Eukaryota"
}

GROUP_COLORS = {
    "Amniota": "#2b5c8f",
    "Eukaryota": "#3282b8",
    "Actinobacteria": "#e07a5f",
    "Bacteria": "#d9534f"
}

METHOD_LABELS = {
    "msa_ml": "MSA+ML (IQ-TREE 2)",
    "msa_nj": "MSA+NJ (RapidNJ)",
    "psa_nj": "PSA+NJ (RapidNJ)",
    "gs": "GS (Graph Splitting)"
}


def extract_branch_lengths_from_tree(tree_path: Path):
    """Extracts all branch lengths from a Newick tree file."""
    try:
        tree = Phylo.read(str(tree_path), "newick")
        lengths = [c.branch_length for c in tree.find_clades() if c.branch_length is not None]
        return lengths
    except Exception as e:
        print(f"Error parsing {tree_path}: {e}", file=sys.stderr)
        return []


def collect_all_branch_lengths():
    """Traverse all tree files and collect branch lengths and per-HOG summaries."""
    records_hog = []
    all_lengths_by_method_group = {m: {g: [] for g in GROUP_NAMES.values()} for m in METHOD_LABELS.keys()}
    all_lengths_by_method = {m: [] for m in METHOD_LABELS.keys()}

    # Also load SSS and TCS if available to enrich per-HOG table
    tcs_file = RESULTS_DIR / "dataset_tcs_summary.csv"
    df_tcs = pd.read_csv(tcs_file) if tcs_file.exists() else None
    
    hog_sim_file = RESULTS_DIR / "hog_similarity_summary.csv"
    df_sim = pd.read_csv(hog_sim_file) if hog_sim_file.exists() else None

    methods = ["msa_ml", "msa_nj", "psa_nj", "gs"]

    # Gather HOG IDs from msa_ml directory
    ml_files = sorted((TREES_DIR / "msa_ml").glob("*.nwk"))
    
    for ml_file in ml_files:
        stem = ml_file.stem
        # e.g. act_E0906941
        parts = stem.split("_")
        dataset_code = parts[0]
        hog_id = "_".join(parts[1:])
        group_name = GROUP_NAMES.get(dataset_code, dataset_code)

        # Retrieve SSS and TCS
        sss_val = np.nan
        if df_sim is not None:
            match_sim = df_sim[(df_sim["dataset"] == dataset_code) & (df_sim["hog_id"] == hog_id)]
            if not match_sim.empty:
                sss_val = match_sim["mean_sss"].values[0]

        row_data = {
            "dataset": dataset_code,
            "dataset_name": group_name,
            "hog_id": hog_id,
            "mean_sss": sss_val
        }

        for method in methods:
            tpath = TREES_DIR / method / f"{stem}.nwk"
            if not tpath.exists():
                continue
            bls = extract_branch_lengths_from_tree(tpath)

            if len(bls) > 0:
                bl_arr = np.array(bls)
                all_lengths_by_method[method].extend(bls)
                all_lengths_by_method_group[method][group_name].extend(bls)

                row_data[f"{method}_num_branches"] = len(bl_arr)
                row_data[f"{method}_tree_len"] = np.sum(bl_arr)
                row_data[f"{method}_min"] = np.min(bl_arr)
                row_data[f"{method}_max"] = np.max(bl_arr)
                row_data[f"{method}_mean"] = np.mean(bl_arr)
                row_data[f"{method}_median"] = np.median(bl_arr)
                row_data[f"{method}_std"] = np.std(bl_arr)
                row_data[f"{method}_q25"] = np.percentile(bl_arr, 25)
                row_data[f"{method}_q75"] = np.percentile(bl_arr, 75)
                row_data[f"{method}_iqr"] = np.percentile(bl_arr, 75) - np.percentile(bl_arr, 25)
                row_data[f"{method}_skew"] = stats.skew(bl_arr) if len(bl_arr) > 2 else 0.0
                row_data[f"{method}_neg_count"] = int(np.sum(bl_arr < 0))
                row_data[f"{method}_zero_or_micro_count"] = int(np.sum(bl_arr <= 1e-5))
            else:
                row_data[f"{method}_num_branches"] = 0
                row_data[f"{method}_tree_len"] = np.nan
                row_data[f"{method}_min"] = np.nan
                row_data[f"{method}_max"] = np.nan
                row_data[f"{method}_mean"] = np.nan
                row_data[f"{method}_median"] = np.nan
                row_data[f"{method}_std"] = np.nan
                row_data[f"{method}_q25"] = np.nan
                row_data[f"{method}_q75"] = np.nan
                row_data[f"{method}_iqr"] = np.nan
                row_data[f"{method}_skew"] = np.nan
                row_data[f"{method}_neg_count"] = 0
                row_data[f"{method}_zero_or_micro_count"] = 0

        # Also add TCS scores if available
        if df_tcs is not None:
            match_tcs = df_tcs[(df_tcs["dataset"] == dataset_code) & (df_tcs["hog_id"] == hog_id)]
            if not match_tcs.empty:
                for m in methods:
                    col = f"tcs_{m}"
                    if col in match_tcs.columns:
                        row_data[col] = match_tcs[col].values[0]

        records_hog.append(row_data)

    df_hog_summary = pd.DataFrame(records_hog)
    return df_hog_summary, all_lengths_by_method, all_lengths_by_method_group


def compute_distribution_stats(arr):
    """Compute comprehensive summary stats for an array of branch lengths."""
    if len(arr) == 0:
        return {
            "count": 0, "min": np.nan, "max": np.nan, "mean": np.nan, "median": np.nan,
            "std": np.nan, "q05": np.nan, "q25": np.nan, "q75": np.nan, "q95": np.nan,
            "q99": np.nan, "neg_frac": np.nan, "micro_frac": np.nan
        }
    a = np.array(arr)
    return {
        "count": len(a),
        "min": float(np.min(a)),
        "max": float(np.max(a)),
        "mean": float(np.mean(a)),
        "median": float(np.median(a)),
        "std": float(np.std(a)),
        "q05": float(np.percentile(a, 5)),
        "q25": float(np.percentile(a, 25)),
        "q75": float(np.percentile(a, 75)),
        "q95": float(np.percentile(a, 95)),
        "q99": float(np.percentile(a, 99)),
        "neg_frac": float(np.mean(a < 0)),
        "micro_frac": float(np.mean(a <= 1e-5))
    }


def perform_clustering_analysis(df_hog: pd.DataFrame):
    """
    Performs clustering of branch length distributions across HOGs.
    Uses standardized distribution features (mean, median, std, min, max, q25, q75, skew).
    """
    # We focus on MSA+ML (the standard reference) and compare with MSA+NJ / PSA+NJ
    features_ml = [
        "msa_ml_mean", "msa_ml_median", "msa_ml_std",
        "msa_ml_min", "msa_ml_max", "msa_ml_q25", "msa_ml_q75", "msa_ml_skew"
    ]
    
    # Drop rows with NaN
    valid_df = df_hog.dropna(subset=features_ml).copy()
    X = valid_df[features_ml].values

    # Standardize features (Z-score)
    mean_X = np.mean(X, axis=0)
    std_X = np.std(X, axis=0)
    std_X[std_X == 0] = 1.0
    X_scaled = (X - mean_X) / std_X

    # Hierarchical Clustering (Ward's linkage on Euclidean distance)
    dists = pdist(X_scaled, metric="euclidean")
    linkage_matrix = sch.linkage(dists, method="ward")

    # Form flat clusters (k=3 or k=4)
    # Using k=4 clusters
    k = 4
    cluster_labels = sch.fcluster(linkage_matrix, t=k, criterion="maxclust")
    valid_df["cluster"] = cluster_labels

    # SVD for PCA
    U, S, Vt = np.linalg.svd(X_scaled, full_matrices=False)
    pca_proj = np.dot(X_scaled, Vt.T[:, :2])
    valid_df["PC1"] = pca_proj[:, 0]
    valid_df["PC2"] = pca_proj[:, 1]
    var_exp = (S[:2] ** 2) / np.sum(S ** 2) * 100

    return valid_df, linkage_matrix, var_exp


def plot_visualizations(df_hog: pd.DataFrame, valid_df: pd.DataFrame, linkage_matrix, var_exp,
                        all_lengths_by_method_group):
    """Generates visualization plots for distribution and clustering."""
    
    # 1. Branch Length Distribution Violin & Box Plots across biological groups
    fig, axes = plt.subplots(1, 3, figsize=(18, 6), sharey=True)
    methods = ["msa_ml", "msa_nj", "psa_nj"]
    
    for idx, m in enumerate(methods):
        ax = axes[idx]
        plot_data = []
        for g_name, lens in all_lengths_by_method_group[m].items():
            for val in lens:
                plot_data.append({"Group": g_name, "BranchLength": val})
        
        df_p = pd.DataFrame(plot_data)
        
        # We clip extreme display values for visualization readability (e.g. 99.5th percentile)
        upper_limit = np.percentile(df_p["BranchLength"], 99.0) * 1.2
        
        sns.violinplot(
            data=df_p,
            x="Group",
            y="BranchLength",
            hue="Group",
            order=["Amniota", "Eukaryota", "Actinobacteria", "Bacteria"],
            palette=GROUP_COLORS,
            legend=False,
            ax=ax,
            inner="quartile",
            cut=0
        )
        ax.set_title(METHOD_LABELS[m], fontsize=13, fontweight="bold")
        ax.set_xlabel("Biological Group", fontsize=11)
        ax.set_ylabel("Branch Length (substitutions/site)" if idx == 0 else "")
        ax.set_ylim(-0.05, min(upper_limit, 2.5))
        
    plt.suptitle("Branch Length Distributions across Biological Groups (by Method)", fontsize=15, y=1.02)
    plt.tight_layout()
    plt.savefig(RESULTS_DIR / "branch_length_distributions_by_group.png", dpi=300, bbox_inches="tight")
    plt.close()

    # 2. Clustering Visualization (Dendrogram + PCA scatter plot)
    fig, (ax_dendro, ax_pca) = plt.subplots(1, 2, figsize=(18, 7))

    # Dendrogram
    sch.dendrogram(
        linkage_matrix,
        truncate_mode="lastp",
        p=25,
        leaf_rotation=90,
        leaf_font_size=10,
        show_contracted=True,
        ax=ax_dendro
    )
    ax_dendro.set_title("Hierarchical Clustering Dendrogram (Ward Linkage, Top 25 Merges)", fontsize=13, fontweight="bold")
    ax_dendro.set_xlabel("HOG Clusters (Contracted Leaf Nodes)")
    ax_dendro.set_ylabel("Euclidean Distance (Scaled Features)")

    # PCA Scatter
    cluster_palette = {
        1: "#e41a1c",  # red
        2: "#377eb8",  # blue
        3: "#4daf4a",  # green
        4: "#984ea3"   # purple
    }
    
    for cl_id in sorted(valid_df["cluster"].unique()):
        sub = valid_df[valid_df["cluster"] == cl_id]
        ax_pca.scatter(
            sub["PC1"], sub["PC2"],
            label=f"Cluster {cl_id} (n={len(sub)})",
            color=cluster_palette.get(cl_id, "#555"),
            alpha=0.7,
            s=40,
            edgecolors="none"
        )
        
    ax_pca.set_title(f"PCA Projection of Branch Length Features\n(PC1: {var_exp[0]:.1f}%, PC2: {var_exp[1]:.1f}%)", fontsize=13, fontweight="bold")
    ax_pca.set_xlabel("Principal Component 1 (Overall Tree Rate / Length Scale)")
    ax_pca.set_ylabel("Principal Component 2 (Variance & Skewness / Rate Heterogeneity)")
    ax_pca.legend(title="Branch Length Clusters", loc="best")

    plt.tight_layout()
    plt.savefig(RESULTS_DIR / "branch_length_clustering.png", dpi=300, bbox_inches="tight")
    plt.close()
    
    print("Visualizations saved successfully.")


def main():
    print("Collecting branch lengths from all tree files...")
    df_hog_summary, all_lengths_by_method, all_lengths_by_method_group = collect_all_branch_lengths()

    # Save per-HOG summary CSV
    out_hog_csv = RESULTS_DIR / "branch_length_hog_summary.csv"
    df_hog_summary.to_csv(out_hog_csv, index=False)
    print(f"Saved per-HOG branch length summary to {out_hog_csv}")

    # Compute overall summary table
    overall_rows = []
    for method, label in METHOD_LABELS.items():
        st = compute_distribution_stats(all_lengths_by_method[method])
        st["method"] = method
        st["method_label"] = label
        overall_rows.append(st)
    
    df_overall = pd.DataFrame(overall_rows)
    out_overall_csv = RESULTS_DIR / "branch_length_overall_summary.csv"
    df_overall.to_csv(out_overall_csv, index=False)
    print(f"Saved overall branch length summary to {out_overall_csv}")

    # Compute group breakdown table
    group_rows = []
    for method, label in METHOD_LABELS.items():
        for group in ["Amniota", "Eukaryota", "Actinobacteria", "Bacteria"]:
            st = compute_distribution_stats(all_lengths_by_method_group[method][group])
            st["method"] = method
            st["method_label"] = label
            st["group"] = group
            group_rows.append(st)

    df_group = pd.DataFrame(group_rows)
    out_group_csv = RESULTS_DIR / "branch_length_group_summary.csv"
    df_group.to_csv(out_group_csv, index=False)
    print(f"Saved group branch length summary to {out_group_csv}")

    # Perform clustering analysis
    print("Performing branch length distribution clustering...")
    valid_df, linkage_matrix, var_exp = perform_clustering_analysis(df_hog_summary)
    
    out_cluster_csv = RESULTS_DIR / "branch_length_clusters.csv"
    valid_df.to_csv(out_cluster_csv, index=False)
    print(f"Saved clustering results to {out_cluster_csv}")

    # Plot visualizations
    print("Generating visualizations...")
    plot_visualizations(df_hog_summary, valid_df, linkage_matrix, var_exp, all_lengths_by_method_group)

    print("\n--- SUMMARY OF RESULTS ---")
    print(df_overall[["method_label", "count", "min", "max", "mean", "median", "std", "q05", "q95", "neg_frac"]].to_string(index=False))


if __name__ == "__main__":
    main()
