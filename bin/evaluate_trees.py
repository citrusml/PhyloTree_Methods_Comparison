#!/usr/bin/env python3
"""
Evaluate Trees Script: evaluate_trees.py

Calculates tree topology evaluation metrics between True Tree and Estimated Tree
according to Matsui & Iwasaki (2020, Systematic Biology, 69(2):265-279):
1. RCRB (Ratio of Correctly Recovered Branches / Correctly Recovered Branch Ratio):
       RCRB = 1 - (RF + O) / (2 * (M - 3)) = correct_branches / (M - 3)
2. RF (Robinson-Foulds unrooted distance) and Normalized RF (nRF).
3. Precision: correct_branches / est_branches.
4. DS (Shortest path branch distance between terminal taxa, e.g. T1 and TN, for LBA analysis).

Appends result row to CSV log.
"""

from typing import Optional
import sys
import os
import json
import argparse
import pandas as pd
import dendropy
from dendropy.calculate import treecompare


def load_unrooted_tree(tree_file: str, taxon_namespace: dendropy.TaxonNamespace) -> dendropy.Tree:
    """Loads a Newick tree, unroots it, and aligns taxon namespace."""
    tree = dendropy.Tree.get(
        path=tree_file,
        schema="newick",
        taxon_namespace=taxon_namespace,
        preserve_underscores=True
    )
    tree.is_rooted = False
    tree.deroot()
    return tree


def compute_shortest_path_edges(
    tree: dendropy.Tree,
    taxon1: str,
    taxon2: str
) -> Optional[int]:
    """
    Computes topological shortest path distance (number of edges/branches)
    between two taxa in an unrooted tree via Breadth-First Search (BFS).
    Used for LBA DS metric (shortest path between nodes 1 and M).
    """
    node1 = tree.find_node_with_taxon_label(taxon1)
    node2 = tree.find_node_with_taxon_label(taxon2)
    if not node1 or not node2:
        return None

    visited = {node1: 0}
    queue = [node1]
    while queue:
        cur = queue.pop(0)
        dist = visited[cur]
        if cur == node2:
            return dist
        neighbors = []
        if cur.parent_node:
            neighbors.append(cur.parent_node)
        for child in cur.child_nodes():
            neighbors.append(child)
        for nbr in neighbors:
            if nbr not in visited:
                visited[nbr] = dist + 1
                queue.append(nbr)
    return None


def main() -> None:
    parser = argparse.ArgumentParser(description="Evaluate Tree Topology Accuracy")
    parser.add_argument("--truetree", required=True, help="Path to True Newick Tree file")
    parser.add_argument("--esttree", required=True, help="Path to Estimated Newick Tree file")
    parser.add_argument(
        "--pipeline",
        required=True,
        help="Pipeline name (e.g. PWA+NJ, MSA+NJ, MSA+ML, GS, etc.)"
    )
    parser.add_argument("--distance", type=float, required=True, help="Evolutionary distance D")
    parser.add_argument("--length", type=int, required=True, help="Sequence length L")
    parser.add_argument("--replicate", type=int, required=True, help="Replicate ID")
    parser.add_argument("--alpha", type=float, default=1.0, help="Gamma shape parameter alpha")
    parser.add_argument(
        "--ics_prop",
        type=float,
        default=0.0,
        help="ICS proportion (simulation invariant category site ratio)"
    )
    parser.add_argument("--json", help="Path to ML metadata JSON file")
    parser.add_argument("--outcsv", required=True, help="Output summary CSV file")
    args = parser.parse_args()

    # Load shared taxon namespace from true tree first
    tns = dendropy.TaxonNamespace()
    true_tree = load_unrooted_tree(args.truetree, tns)
    est_tree = load_unrooted_tree(args.esttree, tns)

    # Encode splits
    true_tree.encode_bipartitions()
    est_tree.encode_bipartitions()

    # Compute unrooted Robinson-Foulds distance
    rf = treecompare.symmetric_difference(true_tree, est_tree)
    num_taxa = len(tns)
    max_rf = 2 * (num_taxa - 3) if num_taxa > 3 else 1
    nrf = rf / max_rf if max_rf > 0 else 0.0

    # Extract internal non-trivial bipartitions (splits)
    bp_true = set(bp for bp in true_tree.bipartition_encoding if not bp.is_trivial())
    bp_est = set(bp for bp in est_tree.bipartition_encoding if not bp.is_trivial())

    true_branches = len(bp_true)
    est_branches = len(bp_est)
    branch_diff_o = abs(true_branches - est_branches)
    correct_branches = len(bp_true.intersection(bp_est))

    # Correctly Recovered Branch Ratio (RCRB) matching Matsui & Iwasaki (2020) Page 5:
    # RCRB = 1 - (RF + O) / (2 * (M - 3))
    denom = 2 * (num_taxa - 3)
    if denom > 0:
        rcrb = 1.0 - (rf + branch_diff_o) / denom
        rcrb = max(0.0, min(1.0, rcrb))
    else:
        rcrb = 1.0

    precision = (correct_branches / est_branches) if est_branches > 0 else 0.0

    # LBA DS metric: shortest path between T1 and TN (or first and last taxa)
    taxon_labels = sorted([t.label for t in tns if t.label is not None])
    ds_distance = None
    ds_true = None
    if len(taxon_labels) >= 2:
        first_taxon = taxon_labels[0]
        last_taxon = taxon_labels[-1]
        ds_distance = compute_shortest_path_edges(est_tree, first_taxon, last_taxon)
        ds_true = compute_shortest_path_edges(true_tree, first_taxon, last_taxon)

    # Load ML metadata if available
    best_model = "N/A"
    gamma_alpha = "N/A"
    if args.json and os.path.exists(args.json):
        try:
            with open(args.json, "r") as f:
                meta = json.load(f)
                best_model = meta.get("best_model_bic", meta.get("model", "N/A"))
                gamma_alpha = str(meta.get("gamma_alpha", "N/A"))
        except Exception:
            pass

    record = {
        "alpha": args.alpha,
        "ics_prop": args.ics_prop,
        "distance": args.distance,
        "length": args.length,
        "replicate": args.replicate,
        "pipeline": args.pipeline,
        "num_taxa": num_taxa,
        "rf_distance": rf,
        "nrf_distance": round(nrf, 6),
        "rcrb": round(rcrb, 6),
        "correct_branches": correct_branches,
        "true_branches": true_branches,
        "est_branches": est_branches,
        "branch_diff_o": branch_diff_o,
        "precision": round(precision, 6),
        "ds_distance": ds_distance if ds_distance is not None else "N/A",
        "ds_true": ds_true if ds_true is not None else "N/A",
        "best_model_bic": best_model,
        "gamma_alpha": gamma_alpha
    }

    df = pd.DataFrame([record])
    header = not os.path.exists(args.outcsv)
    df.to_csv(args.outcsv, mode="a", index=False, header=header)

    print(
        f"Evaluated {args.pipeline} (D={args.distance}, L={args.length}, rep={args.replicate}): "
        f"RCRB={rcrb:.4f}, nRF={nrf:.4f}, RF={rf} (correct={correct_branches}/{true_branches})"
    )


if __name__ == "__main__":
    main()
