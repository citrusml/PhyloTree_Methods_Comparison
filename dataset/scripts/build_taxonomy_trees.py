#!/usr/bin/env python3
"""
build_taxonomy_trees.py
========================
Builds ground-truth taxonomy topology Newick trees for all 400 HOGs
based on the NCBI taxonomy lineages in dataset/*_taxa.tsv and the exact
sequence IDs in dataset/data/<group>/<hog_id>.fa.

Also collects all reconstructed trees (MSA+ML, MSA+NJ, PSA+NJ, GS)
and outputs a unified, highly-compressed JSON/JS data file for the
dual tree viewer (left: estimated tree, right: ground-truth topology).
"""

import os
import re
import csv
import json
from pathlib import Path
from typing import Dict, List, Set, Any, Tuple

SCRIPT_DIR = Path(__file__).resolve().parent
DATASET_DIR = SCRIPT_DIR.parent
DATA_DIR = DATASET_DIR / "data"
RESULTS_DIR = DATASET_DIR / "results"
TREES_DIR = RESULTS_DIR / "trees"

GROUPS = [
    ("amn", "Amniota (10 taxa)", "amniota"),
    ("euk", "Eukaryota (50 taxa)", "eukaryota"),
    ("act", "Actinobacteria (50 taxa)", "actinomycete"),
    ("bac", "Bacteria (50 taxa)", "bacteria")
]

METHODS = ["msa_ml", "msa_nj", "psa_nj", "gs"]


class TaxonNode:
    """Trie node representing an NCBI taxonomy rank or taxon."""
    def __init__(self, node_id: str):
        self.node_id = node_id
        self.children: Dict[str, "TaxonNode"] = {}
        self.leaves: List[str] = []

    def add_path(self, path: List[str], leaf_name: str):
        curr = self
        for nid in path:
            if nid not in curr.children:
                curr.children[nid] = TaxonNode(nid)
            curr = curr.children[nid]
        curr.leaves.append(leaf_name)

    def collapse_unbranched(self) -> "TaxonNode":
        """Collapses intermediate internal nodes with single child and no leaves."""
        new_children = {}
        for nid, child in self.children.items():
            collapsed_child = child.collapse_unbranched()
            new_children[nid] = collapsed_child

        self.children = new_children

        # If this node has no leaves and exactly one child, bypass it
        if not self.leaves and len(self.children) == 1:
            return next(iter(self.children.values()))
        return self

    def to_newick(self) -> str:
        """Serializes this subtree to a clean Newick string."""
        parts = []
        for nid in sorted(self.children.keys()):
            child = self.children[nid]
            child_nwk = child.to_newick()
            if child_nwk:
                parts.append(child_nwk)

        for leaf in sorted(self.leaves):
            parts.append(leaf)

        if not parts:
            return ""
        if len(parts) == 1:
            return parts[0]
        return f"({','.join(parts)})"


def load_taxa_lineages(dataset_dir: Path) -> Dict[str, Dict[str, List[str]]]:
    """
    Loads taxonomy lineage path (list of node IDs) for each species code in each group.
    Returns {group: {species_code: [tax_id_1, tax_id_2, ...]}}.
    """
    lineage_maps = {}
    for g_key, _, _ in GROUPS:
        tsv_path = dataset_dir / f"{g_key}_taxa.tsv"
        if not tsv_path.exists():
            continue
        g_map = {}
        with open(tsv_path, "r", encoding="utf-8") as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 3:
                    code = parts[0].strip()
                    lineage_str = parts[2].strip()
                    node_ids = []
                    for item in lineage_str.split(","):
                        item = item.strip()
                        if not item:
                            continue
                        m = re.match(r"^(\d+)", item)
                        if m:
                            node_ids.append(m.group(1))
                    g_map[code] = node_ids
        lineage_maps[g_key] = g_map
    return lineage_maps


def load_fasta_leaves(fasta_path: Path) -> List[str]:
    """Extracts exact leaf names directly from the HOG FASTA file headers."""
    if not fasta_path.exists():
        return []
    leaves = []
    with open(fasta_path, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith(">"):
                header = line[1:].strip().split()[0]
                clean = header.strip("'\"")
                if clean:
                    leaves.append(clean)
    return sorted(list(set(leaves)))


def sanitize_newick(nwk: str) -> str:
    """Removes enclosing single/double quotes around taxon names in Newick."""
    if not nwk:
        return ""
    # RapidNJ outputs 'SPECIES123':0.123
    clean = re.sub(r"['\"]([A-Za-z0-9_]+)['\"]", r"\1", nwk)
    return clean.strip()


def build_ground_truth_tree(leaves: List[str], group_lineages: Dict[str, List[str]]) -> str:
    """
    Builds a ground truth Newick tree for the exact leaf sequences
    using their species NCBI taxonomy lineage.
    """
    root = TaxonNode("root")
    for leaf in leaves:
        sp_code = leaf[:5]
        lineage = group_lineages.get(sp_code, ["unknown", sp_code])
        root.add_path(lineage, leaf)

    collapsed_root = root.collapse_unbranched()
    nwk = collapsed_root.to_newick()
    if not nwk.endswith(";"):
        nwk += ";"
    return nwk


def load_metadata() -> Tuple[Dict[str, Dict[str, Any]], Dict[str, Dict[str, Any]], Dict[str, Dict[str, Any]]]:
    """Loads SSS/Identity summary, TCS summary, and Branch Length Cluster summary tables."""
    meta_sss = {}
    meta_tcs = {}
    meta_cluster = {}

    sss_csv = RESULTS_DIR / "hog_similarity_summary.csv"
    if sss_csv.exists():
        with open(sss_csv, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            for row in reader:
                key = f"{row['dataset']}_{row['hog_id']}"
                meta_sss[key] = {
                    "mean_sss": float(row["mean_sss"]),
                    "median_sss": float(row["median_sss"]),
                    "mean_identity": float(row["mean_identity"]) * 100.0,
                    "num_taxa": int(row["num_taxa"])
                }

    tcs_csv = RESULTS_DIR / "dataset_tcs_summary.csv"
    if tcs_csv.exists():
        with open(tcs_csv, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            for row in reader:
                key = f"{row['dataset']}_{row['hog_id']}"
                meta_tcs[key] = {
                    "msa_ml": float(row.get("tcs_msa_ml", 0.0)),
                    "msa_nj": float(row.get("tcs_msa_nj", 0.0)),
                    "psa_nj": float(row.get("tcs_psa_nj", 0.0)),
                    "gs": float(row.get("tcs_gs", 0.0))
                }

    cluster_csv = RESULTS_DIR / "branch_length_clusters.csv"
    if cluster_csv.exists():
        with open(cluster_csv, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            for row in reader:
                key = f"{row['dataset']}_{row['hog_id']}"
                meta_cluster[key] = {
                    "cluster": int(row.get("cluster", 0)),
                    "msa_ml_mean": float(row.get("msa_ml_mean", 0.0)),
                    "msa_ml_max": float(row.get("msa_ml_max", 0.0)),
                    "msa_ml_skew": float(row.get("msa_ml_skew", 0.0)),
                    "msa_ml_tree_len": float(row.get("msa_ml_tree_len", 0.0))
                }

    return meta_sss, meta_tcs, meta_cluster


def collect_all_data() -> Dict[str, Any]:
    """Scans all trees and metadata, creating the complete viewer dataset."""
    lineages = load_taxa_lineages(DATASET_DIR)
    meta_sss, meta_tcs, meta_cluster = load_metadata()

    dataset_dict = {
        "groups": {k: name for k, name, _ in GROUPS},
        "methods": {
            "msa_ml": "MSA+ML (IQ-TREE 2)",
            "msa_nj": "MSA+NJ (RapidNJ)",
            "psa_nj": "PSA+NJ (RapidNJ)",
            "gs":     "GS (Graph Splitting)"
        },
        "clusters": {
            "1": {
                "id": 1,
                "name": "Cluster 1: 超遠縁・全体長枝型",
                "badge": "C1: 超遠縁",
                "color": "#ef4444",
                "bg": "rgba(239, 68, 68, 0.18)",
                "description": "全枝が一様に長い置換飽和領域。中央値TCSでNJ法がMLを逆転。",
                "top_method": "MSA+ML (35.5%) / NJ連合 (48.4%)"
            },
            "2": {
                "id": 2,
                "name": "Cluster 2: 標準遠縁・伸長型",
                "badge": "C2: 標準遠縁",
                "color": "#f59e0b",
                "bg": "rgba(245, 158, 11, 0.18)",
                "description": "原核生物中心。枝が均等伸長。GS法が躍進(30.6%)、PSA+NJがMSA+NJを圧倒。",
                "top_method": "ML (30.8%) vs GS (30.6%) / 順位首位: PSA+NJ"
            },
            "3": {
                "id": 3,
                "name": "Cluster 3: 保存的・短枝型",
                "badge": "C3: 高保存",
                "color": "#10b981",
                "bg": "rgba(16, 185, 129, 0.18)",
                "description": "羊膜類・真核生物中心。機能制約が強く高保存。全手法が高精度で安定。",
                "top_method": "MSA+ML (37.2%)"
            },
            "4": {
                "id": 4,
                "name": "Cluster 4: 局所長枝・LBAリスク型",
                "badge": "C4: LBA注意",
                "color": "#a855f7",
                "bg": "rgba(168, 85, 247, 0.18)",
                "description": "最大枝長3.14, 歪度4.23と極端な長枝を持つ。ML法がLBAを回避して他を圧倒。",
                "top_method": "MSA+ML (勝率 40.4% / 平均順位 2.12位 独走)"
            }
        },
        "hogs": {}
    }

    for g_key, _, g_dir in GROUPS:
        group_lineage = lineages.get(g_key, {})
        fasta_folder = DATA_DIR / g_dir

        ml_dir = TREES_DIR / "msa_ml"
        hog_files = sorted(ml_dir.glob(f"{g_key}_*.nwk"))

        for hf in hog_files:
            hog_id = hf.stem[len(g_key) + 1:]  # strip 'act_' -> 'E0906941'
            unique_key = f"{g_key}_{hog_id}"

            # 1. Read exact leaves from source FASTA file
            fasta_path = fasta_folder / f"{hog_id}.fa"
            if not fasta_path.exists():
                fasta_path = fasta_folder / f"HOG_{hog_id}.fa"
            fasta_leaves = load_fasta_leaves(fasta_path)

            # 2. Read and sanitize all 4 method trees
            method_trees = {}
            for m in METHODS:
                m_path = TREES_DIR / m / f"{unique_key}.nwk"
                if m_path.exists():
                    with open(m_path, "r", encoding="utf-8") as tf:
                        raw_nwk = tf.read().strip()
                        method_trees[m] = sanitize_newick(raw_nwk)
                else:
                    method_trees[m] = ""

            # 3. Build Ground Truth taxonomy topology tree using exact FASTA leaves
            gt_tree = build_ground_truth_tree(fasta_leaves, group_lineage)

            # Metadata
            s_info = meta_sss.get(unique_key, {})
            t_info = meta_tcs.get(unique_key, {})
            c_info = meta_cluster.get(unique_key, {})

            dataset_dict["hogs"][unique_key] = {
                "group": g_key,
                "hog_id": hog_id,
                "num_leaves": len(fasta_leaves),
                "mean_sss": round(s_info.get("mean_sss", 0.0), 4),
                "median_sss": round(s_info.get("median_sss", 0.0), 4),
                "mean_identity": round(s_info.get("mean_identity", 0.0), 2),
                "cluster": c_info.get("cluster", 0),
                "branch_stats": {
                    "ml_mean": round(c_info.get("msa_ml_mean", 0.0), 3),
                    "ml_max": round(c_info.get("msa_ml_max", 0.0), 3),
                    "ml_skew": round(c_info.get("msa_ml_skew", 0.0), 2),
                    "ml_tree_len": round(c_info.get("msa_ml_tree_len", 0.0), 2)
                },
                "tcs": {
                    "msa_ml": round(t_info.get("msa_ml", 0.0), 2),
                    "msa_nj": round(t_info.get("msa_nj", 0.0), 2),
                    "psa_nj": round(t_info.get("psa_nj", 0.0), 2),
                    "gs":     round(t_info.get("gs", 0.0), 2)
                },
                "ground_truth_tree": gt_tree,
                "trees": method_trees
            }

    return dataset_dict


def main():
    print("=== Building Ground Truth Taxonomy Trees & Exporting Viewer Data ===")
    data = collect_all_data()
    total_hogs = len(data["hogs"])
    print(f"Total HOGs processed: {total_hogs}")

    json_path = RESULTS_DIR / "dual_tree_data.json"
    js_path = RESULTS_DIR / "dual_tree_data.js"

    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(data, f, separators=(',', ':'))
    json_size_mb = json_path.stat().st_size / (1024 * 1024)
    print(f"Exported JSON: {json_path} ({json_size_mb:.2f} MB)")

    with open(js_path, "w", encoding="utf-8") as f:
        f.write("/* Auto-generated Dual Tree Dataset */\nwindow.DUAL_TREE_DATA = ")
        json.dump(data, f, separators=(',', ':'))
        f.write(";\n")
    js_size_mb = js_path.stat().st_size / (1024 * 1024)
    print(f"Exported JS:   {js_path} ({js_size_mb:.2f} MB)")

    print("Success! Data aggregation complete.")


if __name__ == "__main__":
    main()
