#!/usr/bin/env python3
"""
calculate_dataset_tcs.py
========================
Calculates the Taxonomy Congruence Score (TCS / Taxonomy Overlap Score)
for phylogenetic trees reconstructed by the 4 methods (MSA+ML, MSA+NJ, PSA+NJ, GS)
across the 400 real biological HOGs (Amniota, Eukaryota, Bacteria, Actinobacteria).

Adheres strictly to the TCS algorithm in DessimozLab/ampliphy-analysis (Moi et al., 2023; Kim et al.):
  - Traverses tree hierarchy to compute shared taxonomic lineage overlap at each internal node.
  - Scores taxonomic consistency against the NCBI taxonomy definitions in dataset/*_taxa.tsv.
  - Normalized formula: (root.tax_score - root.leaf_acc * len(overlap)) / root.leaf_size

Generates:
  1. dataset/results/dataset_tcs_summary.csv
  2. Publication-quality scatterplots (300 DPI):
     - dataset/results/scatter_sss_vs_tcs_by_group.png (2x2 grid for amn, euk, act, bac)
     - dataset/results/scatter_sss_vs_tcs_amn.png
     - dataset/results/scatter_sss_vs_tcs_euk.png
     - dataset/results/scatter_sss_vs_tcs_act.png
     - dataset/results/scatter_sss_vs_tcs_bac.png
"""

import sys
import time
from pathlib import Path
from typing import Dict, Set, List, Any

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import toytree

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent.parent
DATASET_DIR = SCRIPT_DIR.parent
RESULTS_DIR = DATASET_DIR / "results"
TREES_DIR = RESULTS_DIR / "trees"

plt.rcParams.update({
    "font.sans-serif": ["Arial", "DejaVu Sans", "Helvetica", "sans-serif"],
    "font.size": 11,
    "axes.labelsize": 12,
    "axes.titlesize": 13,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 10.5,
    "figure.titlesize": 15,
    "lines.linewidth": 1.8,
    "axes.grid": True,
    "grid.alpha": 0.3,
    "grid.linestyle": "--"
})

GROUPS = [
    ("amn", "Amniota (10 taxa)"),
    ("euk", "Eukaryota (50 taxa)"),
    ("act", "Actinobacteria (50 taxa)"),
    ("bac", "Bacteria (50 taxa)")
]

METHODS = [
    ("msa_ml", "MSA+ML (IQ-TREE 2)", "#2980b9", "o"),
    ("msa_nj", "MSA+NJ (RapidNJ)",   "#e67e22", "s"),
    ("psa_nj", "PSA+NJ (RapidNJ)",   "#27ae60", "^"),
    ("gs",     "GS (Graph Splitting)","#c0392b", "D")
]


def load_taxa_lineages(dataset_dir: Path) -> Dict[str, Dict[str, Set[str]]]:
    """
    Loads taxon mnemonic codes and their lineage sets from dataset/*_taxa.tsv.
    Returns {group: {code: set(taxonomic_ranks)}}.
    """
    taxa_maps = {}
    for g_key in ["amn", "euk", "bac", "act"]:
        tsv_path = dataset_dir / f"{g_key}_taxa.tsv"
        if not tsv_path.exists():
            print(f"Warning: {tsv_path} not found.")
            continue

        g_map = {}
        with open(tsv_path, "r", encoding="utf-8") as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 3:
                    code = parts[0].strip()
                    lineage_str = parts[2].strip()
                    g_map[code] = set(x.strip() for x in lineage_str.split(",") if x.strip())
        taxa_maps[g_key] = g_map
    return taxa_maps


from loguru import logger
logger.disable("toytree")


def compute_tree_tcs(tree_newick: str, taxa_map: Dict[str, Set[str]], apply_rooting: bool = True, method: str = "") -> float:
    """
    Computes the Taxonomy Congruence Score (TCS) on a Newick tree
    matching DessimozLab/ampliphy-analysis (tcs_batch.py and scripts/run_tree).

    If apply_rooting is True:
      - For methods with branch lengths (msa_ml, msa_nj, psa_nj), applies MAD rooting
        (Minimal Ancestor Deviation) using toytree.mod.root_on_minimal_ancestor_deviation(),
        falling back to midpoint rooting if needed.
      - For GS (Graph Splitting), the tree is inherently rooted via recursive top-down
        spectral bisection.
    """
    if not tree_newick or not tree_newick.strip():
        return np.nan

    try:
        t = toytree.tree(tree_newick.strip())
    except Exception:
        return np.nan

    if apply_rooting and method != "gs":
        try:
            t = t.mod.root_on_minimal_ancestor_deviation()
        except Exception:
            try:
                t = t.mod.root_on_midpoint()
            except Exception:
                pass

    # Assign lineage set to each leaf
    for n in t.treenode.iter_leaves():
        clean_name = n.name.strip("'\"")
        code = clean_name[:5]
        n.lineage = taxa_map.get(code, set())

    def getTaxOverlap(node):
        if node.is_leaf():
            node.leaf_size = 1
            node.leaf_acc = 0
            node.tax_score = 0
            return node.lineage
        else:
            sets = []
            leaf_size = 0
            leaf_acc = 0
            tax_score = 0
            for c in node.children:
                sets.append(getTaxOverlap(c))
                leaf_size += c.leaf_size
                leaf_acc += c.leaf_acc
                tax_score += c.tax_score

            sets = [s for s in sets if s]
            if len(sets) > 0:
                nset = sets[0].copy()
                for cset in sets[1:]:
                    nset = nset.intersection(cset)
                node.tax_score = tax_score + len(nset) * leaf_size
            else:
                nset = set()
                node.tax_score = tax_score

            node.leaf_size = leaf_size
            node.leaf_acc = leaf_acc + leaf_size
            return nset

    root = t.treenode
    overlap = getTaxOverlap(root)
    if root.leaf_size <= 0:
        return 0.0

    adjscore = (root.tax_score - root.leaf_acc * len(overlap)) / root.leaf_size
    return float(adjscore)


def process_all_trees(taxa_maps: Dict[str, Dict[str, Set[str]]]) -> pd.DataFrame:
    """Calculates TCS scores (both MAD-rooted and unrooted raw) for all 400 HOGs across 4 methods."""
    records = []

    for d_key, _ in GROUPS:
        taxa_map = taxa_maps.get(d_key, {})
        hogs_file = DATASET_DIR / f"{d_key}_hogs.tsv"
        with open(hogs_file, "r") as f:
            hog_ids = [line.strip().split("\t")[0].replace("HOG_", "").replace("HOG:", "") for line in f if line.strip()]

        for h_id in hog_ids:
            row = {
                "dataset": d_key,
                "hog_id": h_id
            }
            for m_key, _, _, _ in METHODS:
                tpath = TREES_DIR / m_key / f"{d_key}_{h_id}.nwk"
                if tpath.exists() and tpath.stat().st_size > 0:
                    with open(tpath, "r") as f:
                        tree_str = f.read()
                    score_rooted = compute_tree_tcs(tree_str, taxa_map, apply_rooting=True, method=m_key)
                    score_raw = compute_tree_tcs(tree_str, taxa_map, apply_rooting=False, method=m_key)
                    row[f"tcs_{m_key}"] = round(score_rooted, 4)
                    row[f"tcs_raw_{m_key}"] = round(score_raw, 4)
                else:
                    row[f"tcs_{m_key}"] = np.nan
                    row[f"tcs_raw_{m_key}"] = np.nan
            records.append(row)

    df_tcs = pd.DataFrame(records)

    # Merge with SSS and identity
    sim_file = RESULTS_DIR / "hog_similarity_summary.csv"
    if sim_file.exists():
        df_sim = pd.read_csv(sim_file)
        cols_merge = ["dataset", "hog_id", "mean_sss", "median_sss", "mean_identity"]
        df_merged = pd.merge(df_tcs, df_sim[cols_merge], on=["dataset", "hog_id"], how="left")
    else:
        df_merged = df_tcs

    return df_merged


def plot_combined_scatter(df: pd.DataFrame, out_path: Path):
    """
    Plots a 2x2 multi-panel figure where each panel is a biological group (amn, euk, act, bac):
      - X-axis: Sequence Similarity Score (SSS)
      - Y-axis: Taxonomy Congruence Score (TCS)
      - Hue: Method (MSA+ML, MSA+NJ, PSA+NJ, GS)
    """
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    axes = axes.flatten()

    for idx, (d_key, title_name) in enumerate(GROUPS):
        ax = axes[idx]
        sub = df[df["dataset"] == d_key]

        for m_key, m_label, color, marker in METHODS:
            col = f"tcs_{m_key}"
            valid = sub.dropna(subset=["mean_sss", col])

            ax.scatter(
                valid["mean_sss"],
                valid[col],
                label=m_label,
                color=color,
                marker=marker,
                alpha=0.65,
                s=48,
                edgecolors="black",
                linewidth=0.4
            )

            # Trend line
            if len(valid) > 2:
                sns.regplot(
                    data=valid,
                    x="mean_sss",
                    y=col,
                    scatter=False,
                    ax=ax,
                    color=color,
                    line_kws={"linewidth": 2.2, "alpha": 0.85},
                    ci=None
                )

        ax.set_title(title_name, fontweight="bold", fontsize=13)
        ax.set_xlabel("Sequence Similarity Score (SSS: $\\bar{w}$)")
        ax.set_ylabel("Taxonomy Congruence Score (TCS)")
        ax.set_xlim(-0.02, 1.02)
        ax.legend(frameon=True, loc="upper left", framealpha=0.9)

    plt.suptitle("Taxonomic Congruence (TCS) vs. Sequence Similarity (SSS) Across Biological Groups", fontweight="bold", y=0.98)
    plt.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved combined plot: {out_path}")


def plot_individual_scatters(df: pd.DataFrame):
    """Generates individual high-resolution scatterplots for each biological group."""
    for d_key, title_name in GROUPS:
        fig, ax = plt.subplots(figsize=(9, 6.5))
        sub = df[df["dataset"] == d_key]

        for m_key, m_label, color, marker in METHODS:
            col = f"tcs_{m_key}"
            valid = sub.dropna(subset=["mean_sss", col])

            ax.scatter(
                valid["mean_sss"],
                valid[col],
                label=m_label,
                color=color,
                marker=marker,
                alpha=0.7,
                s=55,
                edgecolors="black",
                linewidth=0.4
            )

            if len(valid) > 2:
                sns.regplot(
                    data=valid,
                    x="mean_sss",
                    y=col,
                    scatter=False,
                    ax=ax,
                    color=color,
                    line_kws={"linewidth": 2.4, "alpha": 0.85},
                    ci=None
                )

        ax.set_title(f"{title_name}: SSS vs. TCS (Taxonomy Congruence Score)", fontweight="bold", fontsize=13)
        ax.set_xlabel("Sequence Similarity Score (SSS: $\\bar{w}$)", fontsize=12)
        ax.set_ylabel("Taxonomy Congruence Score (TCS)", fontsize=12)
        ax.set_xlim(-0.02, 1.02)
        ax.legend(frameon=True, loc="upper left", framealpha=0.9)
        plt.tight_layout()

        out_name = RESULTS_DIR / f"scatter_sss_vs_tcs_{d_key}.png"
        fig.savefig(out_name, dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(f"Saved individual plot: {out_name}")


def main():
    print("=" * 80)
    print("  Calculating Taxonomy Congruence Score (TCS) — AmpliPhy Methodology")
    print(f"  Dataset Dir : {DATASET_DIR}")
    print(f"  Trees Dir   : {TREES_DIR}")
    print(f"  Results Dir : {RESULTS_DIR}")
    print("=" * 80)

    t0 = time.time()
    taxa_maps = load_taxa_lineages(DATASET_DIR)
    df_results = process_all_trees(taxa_maps)

    csv_path = RESULTS_DIR / "dataset_tcs_summary.csv"
    df_results.to_csv(csv_path, index=False)
    print(f"\nSaved TCS Summary: {csv_path} ({len(df_results)} rows in {time.time() - t0:.2f}s)")

    # Print mean TCS summary by group
    cols_tcs = [f"tcs_{m[0]}" for m in METHODS]
    cols_tcs_raw = [f"tcs_raw_{m[0]}" for m in METHODS]

    print("\n=== Unrooted (Raw Outer-Bracket Root) Mean TCS ===")
    print(df_results.groupby("dataset")[cols_tcs_raw].mean().round(3).to_string())

    print("\n=== MAD-Rooted (AmpliPhy Pipeline Standard) Mean TCS (Higher is Better) ===")
    print(df_results.groupby("dataset")[cols_tcs].mean().round(3).to_string())

    # Plot
    print("\nGenerating Scatterplots (X=SSS, Y=TCS, Hue=Method)...")
    out_combined = RESULTS_DIR / "scatter_sss_vs_tcs_by_group.png"
    plot_combined_scatter(df_results, out_combined)
    plot_individual_scatters(df_results)
    print("\nAll TCS calculations and scatterplots completed successfully!")


if __name__ == "__main__":
    main()
