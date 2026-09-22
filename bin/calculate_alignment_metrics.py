#!/usr/bin/env python3
"""
Alignment Metrics Script: calculate_alignment_metrics.py

Calculates alignment accuracy metrics (Sum-of-Pairs / SP score, Modeler score),
structural properties (length, gap ratio), and optional true patristic distance matrix
comparing True MSA (from AliSim/INDELible) and Estimated MSA (from MAFFT).
"""

from typing import Dict, List, Tuple, Set, Optional
import os
import sys
import argparse
import pandas as pd
from Bio import SeqIO
import dendropy


def compute_sp_score(
    true_msa_path: str,
    est_msa_path: str
) -> Tuple[Optional[float], Optional[float], int, int, int]:
    """
    Computes Sum-of-Pairs (SP) score (Developer's score f_D) and Modeler score (CS)
    based on homologous residue pairs between True MSA and Estimated MSA.

    Parameters
    ----------
    true_msa_path : str
        Path to reference True MSA (FASTA or PHYLIP format).
    est_msa_path : str
        Path to estimated MSA (FASTA format).

    Returns
    -------
    Tuple[Optional[float], Optional[float], int, int, int]
        (sp_score, modeler_score, correct_pairs, total_true_pairs, total_est_pairs)
    """
    format_true = "fasta" if str(true_msa_path).endswith((".fa", ".fas", ".fasta")) else "phylip"
    t_recs = {r.id: str(r.seq).upper() for r in SeqIO.parse(true_msa_path, format_true)}
    e_recs = {r.id: str(r.seq).upper() for r in SeqIO.parse(est_msa_path, "fasta")}

    taxa = sorted(set(t_recs.keys()).intersection(set(e_recs.keys())))
    if len(taxa) < 2:
        return None, None, 0, 0, 0

    # Extract homologous residue pairs from True MSA
    true_pairs: Set[Tuple[str, int, str, int]] = set()
    for i in range(len(taxa)):
        for j in range(i + 1, len(taxa)):
            t1, t2 = taxa[i], taxa[j]
            seq1, seq2 = t_recs[t1], t_recs[t2]
            pos1, pos2 = 0, 0
            for c1, c2 in zip(seq1, seq2):
                if c1 != "-" and c2 != "-":
                    true_pairs.add((t1, pos1, t2, pos2))
                if c1 != "-":
                    pos1 += 1
                if c2 != "-":
                    pos2 += 1

    total_true = len(true_pairs)
    if total_true == 0:
        return 1.0, 1.0, 0, 0, 0

    # Count matching homologous residue pairs in Estimated MSA
    correct = 0
    total_est = 0
    for i in range(len(taxa)):
        for j in range(i + 1, len(taxa)):
            t1, t2 = taxa[i], taxa[j]
            seq1, seq2 = e_recs[t1], e_recs[t2]
            pos1, pos2 = 0, 0
            for c1, c2 in zip(seq1, seq2):
                if c1 != "-" and c2 != "-":
                    total_est += 1
                    if (t1, pos1, t2, pos2) in true_pairs:
                        correct += 1
                if c1 != "-":
                    pos1 += 1
                if c2 != "-":
                    pos2 += 1

    sp_score = float(correct) / total_true if total_true > 0 else 1.0
    modeler_score = float(correct) / total_est if total_est > 0 else 1.0

    return sp_score, modeler_score, correct, total_true, total_est


def analyze_msa_properties(msa_path: str) -> Dict[str, float]:
    """
    Computes basic alignment properties (length, gap ratio).

    Parameters
    ----------
    msa_path : str
        Path to alignment file.

    Returns
    -------
    Dict[str, float]
        Dictionary of properties: aln_length, gap_ratio, num_taxa.
    """
    fmt = "fasta" if str(msa_path).endswith((".fa", ".fas", ".fasta")) else "phylip"
    records = list(SeqIO.parse(msa_path, fmt))
    if not records:
        return {"aln_length": 0, "gap_ratio": 0.0, "num_taxa": 0}

    num_seqs = len(records)
    aln_len = len(records[0].seq)
    total_cells = aln_len * num_seqs
    total_gaps = sum(str(r.seq).count("-") for r in records)
    gap_ratio = float(total_gaps) / total_cells if total_cells > 0 else 0.0

    return {
        "aln_length": float(aln_len),
        "gap_ratio": round(gap_ratio, 6),
        "num_taxa": float(num_seqs),
    }


def export_true_patristic_matrix(
    tree_path: str,
    out_matrix_path: str,
    ordered_taxa: Optional[List[str]] = None
) -> None:
    """
    Computes exact patristic distance matrix from Newick tree and writes in PHYLIP format.

    Parameters
    ----------
    tree_path : str
        Path to Newick tree.
    out_matrix_path : str
        Path to output PHYLIP distance matrix.
    ordered_taxa : Optional[List[str]], default=None
        Specified order of taxon names.
    """
    tns = dendropy.TaxonNamespace()
    tree = dendropy.Tree.get(
        path=tree_path,
        schema="newick",
        taxon_namespace=tns,
        preserve_underscores=True
    )
    pdm = tree.phylogenetic_distance_matrix()

    if ordered_taxa:
        taxa = [tns.get_taxon(name) for name in ordered_taxa if tns.has_taxon_label(name)]
    else:
        taxa = sorted([t for t in tns if t.label is not None], key=lambda t: t.label)

    n = len(taxa)
    lines = [f"   {n}\n"]
    for t1 in taxa:
        row_str = f"{t1.label:<10}"
        for t2 in taxa:
            dist = pdm(t1, t2)
            row_str += f"  {dist:.6f}"
        lines.append(row_str + "\n")

    os.makedirs(os.path.dirname(os.path.abspath(out_matrix_path)), exist_ok=True)
    with open(out_matrix_path, "w") as f:
        f.writelines(lines)


def main() -> None:
    parser = argparse.ArgumentParser(description="Calculate alignment SP score and structural properties")
    parser.add_argument("--truemsa", required=True, help="Path to True reference MSA (FASTA)")
    parser.add_argument("--estmsa", required=True, help="Path to Estimated MSA (FASTA)")
    parser.add_argument("--truetree", help="Optional path to True Newick Tree")
    parser.add_argument("--outmatrix", help="Optional output True Patristic Matrix (PHYLIP)")
    parser.add_argument("--distance", type=float, required=True, help="Evolutionary distance D")
    parser.add_argument("--length", type=int, required=True, help="Sequence length L")
    parser.add_argument("--replicate", type=int, required=True, help="Replicate ID")
    parser.add_argument("--outcsv", required=True, help="Output CSV path for chunk results")
    args = parser.parse_args()

    # 1. Compute SP score
    sp_score, modeler_score, correct_pairs, total_true, total_est = compute_sp_score(
        args.truemsa,
        args.estmsa
    )

    # 2. Structural properties
    t_props = analyze_msa_properties(args.truemsa)
    e_props = analyze_msa_properties(args.estmsa)

    len_ratio = None
    if t_props.get("aln_length", 0) > 0:
        len_ratio = round(e_props.get("aln_length", 0) / t_props["aln_length"], 6)

    # 3. Optional True Patristic Distance Matrix export
    if args.truetree and args.outmatrix:
        export_true_patristic_matrix(args.truetree, args.outmatrix)

    # 4. Record row
    row = {
        "distance": args.distance,
        "length": args.length,
        "replicate": args.replicate,
        "sp_score": round(sp_score, 6) if sp_score is not None else None,
        "modeler_score": round(modeler_score, 6) if modeler_score is not None else None,
        "correct_residue_pairs": correct_pairs,
        "true_residue_pairs": total_true,
        "est_residue_pairs": total_est,
        "true_msa_length": int(t_props.get("aln_length", 0)),
        "est_msa_length": int(e_props.get("aln_length", 0)),
        "msa_length_ratio": len_ratio,
        "true_gap_ratio": t_props.get("gap_ratio"),
        "est_gap_ratio": e_props.get("gap_ratio"),
    }

    df_new = pd.DataFrame([row])
    os.makedirs(os.path.dirname(os.path.abspath(args.outcsv)), exist_ok=True)
    if os.path.exists(args.outcsv):
        df_new.to_csv(args.outcsv, mode="a", header=False, index=False)
    else:
        df_new.to_csv(args.outcsv, mode="w", header=True, index=False)

    print(f"[calculate_alignment_metrics] D={args.distance} L={args.length} rep={args.replicate}: SP={sp_score:.4f} -> {args.outcsv}")


if __name__ == "__main__":
    main()
