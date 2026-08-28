#!/usr/bin/env python3
"""
Evaluate Trees Script
Calculates Robinson-Foulds (RF) distance and Normalized RF (nRF) distance
between True Tree and Estimated Tree using DendroPy.
Appends result row to CSV log.
"""

import sys
import os
import json
import argparse
import pandas as pd
import dendropy
from dendropy.calculate import treecompare

def load_unrooted_tree(tree_file, taxon_namespace):
    """Loads a Newick tree, unroots it, and aligns taxon namespace."""
    tree = dendropy.Tree.get(path=tree_file, schema="newick", taxon_namespace=taxon_namespace, preserve_underscores=True)
    tree.is_rooted = False
    tree.deroot()
    return tree

def main():
    parser = argparse.ArgumentParser(description="Evaluate Tree Topology Accuracy")
    parser.add_argument("--truetree", required=True, help="Path to True Newick Tree file")
    parser.add_argument("--esttree", required=True, help="Path to Estimated Newick Tree file")
    parser.add_argument("--pipeline", required=True, help="Pipeline name (e.g. PWA+NJ, MSA+NJ, MSA+ML, PWA+FastME, MSA+FastME)")
    parser.add_argument("--distance", type=float, required=True, help="Evolutionary distance D")
    parser.add_argument("--length", type=int, required=True, help="Sequence length L")
    parser.add_argument("--replicate", type=int, required=True, help="Replicate ID")
    parser.add_argument("--alpha", type=float, default=1.0, help="Gamma shape parameter alpha (simulation true alpha)")
    parser.add_argument("--ics_prop", type=float, default=0.0, help="ICS proportion (simulation invariant category site ratio)")
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
    max_rf = 2 * (num_taxa - 3)
    nrf = rf / max_rf if max_rf > 0 else 0.0

    # Load ML metadata if available
    best_model = "N/A"
    gamma_alpha = "N/A"
    if args.json and os.path.exists(args.json):
        try:
            with open(args.json, "r") as f:
                meta = json.load(f)
                best_model = meta.get("best_model_bic", "N/A")
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
        "best_model_bic": best_model,
        "gamma_alpha": gamma_alpha
    }

    df = pd.DataFrame([record])
    header = not os.path.exists(args.outcsv)
    df.to_csv(args.outcsv, mode="a", index=False, header=header)

    print(f"Evaluated {args.pipeline} (D={args.distance}, L={args.length}, ics_prop={args.ics_prop}, rep={args.replicate}): RF={rf}, nRF={nrf:.4f}")

if __name__ == "__main__":
    main()
