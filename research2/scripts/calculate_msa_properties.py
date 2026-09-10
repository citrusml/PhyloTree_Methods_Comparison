#!/usr/bin/env python3
"""
calculate_msa_properties.py

Calculates comprehensive structural properties of True MSA and MAFFT MSA,
alignment accuracy (Developer's score f_D), sequence length dispersion,
and phylogenetic reconstruction accuracy (RCRB = 1 - nRF) across methods.
Used in research2 comparative benchmark (AliSim vs INDELible).
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd
from Bio import SeqIO
import dendropy
from dendropy.calculate import treecompare


def compute_fd_score(true_msa_path, est_msa_path):
    """
    Computes Developer's score f_D = c / r:
    c = correctly aligned homologous residue pairs in estimated MSA
    r = total aligned residue pairs in reference (true) MSA
    """
    try:
        format_true = "fasta" if str(true_msa_path).endswith((".fa", ".fas", ".fasta")) else "phylip"
        t_recs = {r.id: str(r.seq).upper() for r in SeqIO.parse(true_msa_path, format_true)}
        e_recs = {r.id: str(r.seq).upper() for r in SeqIO.parse(est_msa_path, "fasta")}
        
        taxa = sorted(set(t_recs.keys()).intersection(set(e_recs.keys())))
        if len(taxa) < 2:
            return None

        # Build pair maps for true MSA
        true_pairs = set()
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

        if not true_pairs:
            return 0.0

        # Check in est MSA
        correct = 0
        for i in range(len(taxa)):
            for j in range(i + 1, len(taxa)):
                t1, t2 = taxa[i], taxa[j]
                seq1, seq2 = e_recs[t1], e_recs[t2]
                pos1, pos2 = 0, 0
                for c1, c2 in zip(seq1, seq2):
                    if c1 != "-" and c2 != "-":
                        if (t1, pos1, t2, pos2) in true_pairs:
                            correct += 1
                    if c1 != "-":
                        pos1 += 1
                    if c2 != "-":
                        pos2 += 1

        return float(correct) / len(true_pairs)
    except Exception as e:
        sys.stderr.write(f"Warning: f_D calculation error: {e}\n")
        return None


def calculate_tree_metrics(true_tree_file, est_tree_file):
    """
    Computes RCRB (1 - nRF) and nRF between true tree and estimated tree.
    """
    if not os.path.exists(est_tree_file) or os.path.getsize(est_tree_file) == 0:
        return None, None
    try:
        tns = dendropy.TaxonNamespace()
        true_tr = dendropy.Tree.get(path=str(true_tree_file), schema="newick", taxon_namespace=tns, preserve_underscores=True)
        est_tr = dendropy.Tree.get(path=str(est_tree_file), schema="newick", taxon_namespace=tns, preserve_underscores=True)
        true_tr.is_rooted = False
        est_tr.is_rooted = False
        true_tr.deroot()
        est_tr.deroot()
        true_tr.encode_bipartitions()
        est_tr.encode_bipartitions()

        rf = treecompare.symmetric_difference(true_tr, est_tr)
        m = len(tns)
        max_rf = 2 * (m - 3)
        nrf = rf / max_rf if max_rf > 0 else 0.0
        rcrb = 1.0 - nrf
        return round(rcrb, 6), round(nrf, 6)
    except Exception as e:
        sys.stderr.write(f"Warning: Tree metric calculation error for {est_tree_file}: {e}\n")
        return None, None


def analyze_alignment_structure(msa_path):
    """
    Analyzes alignment length, gap ratio, indel block statistics.
    """
    fmt = "fasta" if str(msa_path).endswith((".fa", ".fas", ".fasta")) else "phylip"
    records = list(SeqIO.parse(msa_path, fmt))
    if not records:
        return {}
    
    num_seqs = len(records)
    aln_len = len(records[0].seq)
    total_cells = aln_len * num_seqs
    
    seq_lens = []
    total_gaps = 0
    indel_block_lengths = []
    
    for r in records:
        s = str(r.seq).upper()
        gaps = s.count("-")
        total_gaps += gaps
        seq_lens.append(len(s) - gaps)
        
        # Analyze consecutive gap blocks
        in_gap = False
        cur_gap_len = 0
        for ch in s:
            if ch == "-":
                in_gap = True
                cur_gap_len += 1
            else:
                if in_gap:
                    indel_block_lengths.append(cur_gap_len)
                    in_gap = False
                    cur_gap_len = 0
        if in_gap:
            indel_block_lengths.append(cur_gap_len)

    gap_ratio = float(total_gaps) / total_cells if total_cells > 0 else 0.0
    num_indel_blocks = len(indel_block_lengths)
    mean_indel_len = float(np.mean(indel_block_lengths)) if indel_block_lengths else 0.0

    return {
        "num_taxa": num_seqs,
        "aln_length": aln_len,
        "gap_ratio": round(gap_ratio, 6),
        "gap_pct": round(gap_ratio * 100.0, 2),
        "mean_seq_len": round(float(np.mean(seq_lens)), 2),
        "std_seq_len": round(float(np.std(seq_lens)), 2),
        "min_seq_len": int(np.min(seq_lens)),
        "max_seq_len": int(np.max(seq_lens)),
        "num_indel_blocks": num_indel_blocks,
        "mean_indel_block_len": round(mean_indel_len, 2),
    }


def main():
    parser = argparse.ArgumentParser(description="Calculate comprehensive MSA properties and tree metrics.")
    parser.add_argument("--cond_name", required=True, help="Condition name (e.g. alisim_rate0.10)")
    parser.add_argument("--simulator", required=True, help="Simulator (alisim or indelible)")
    parser.add_argument("--ins_rate", type=float, required=True, help="Insertion rate")
    parser.add_argument("--del_rate", type=float, required=True, help="Deletion rate")
    parser.add_argument("--dist", type=float, required=True, help="Evolutionary distance (scale)")
    parser.add_argument("--length", type=int, required=True, help="Root sequence length")
    parser.add_argument("--rep", type=int, required=True, help="Replicate ID")
    parser.add_argument("--true_tree", required=True, help="Path to true tree newick")
    parser.add_argument("--true_msa", required=True, help="Path to True MSA FASTA")
    parser.add_argument("--unaligned", required=True, help="Path to unaligned sequences FASTA")
    parser.add_argument("--mafft_msa", required=False, default=None, help="Path to MAFFT MSA FASTA")
    parser.add_argument("--tree_pwa_nj", required=False, default=None, help="PWA+NJ estimated tree")
    parser.add_argument("--tree_msa_nj", required=False, default=None, help="MSA+NJ estimated tree")
    parser.add_argument("--tree_msa_ml", required=False, default=None, help="MSA+ML estimated tree")
    parser.add_argument("--tree_true_msa_nj", required=False, default=None, help="TRUE_MSA+NJ estimated tree")
    parser.add_argument("--tree_true_msa_ml", required=False, default=None, help="TRUE_MSA+ML estimated tree")
    parser.add_argument("--sss_csv", required=False, default=None, help="Path to SSS CSV file")
    parser.add_argument("--outfile", required=True, help="Output CSV path")
    args = parser.parse_args()

    # 1. Analyze True MSA structure
    t_props = analyze_alignment_structure(args.true_msa)
    
    # 2. Analyze MAFFT MSA structure
    m_props = analyze_alignment_structure(args.mafft_msa) if args.mafft_msa and os.path.exists(args.mafft_msa) else {}

    # 3. Compute Developer's score f_D
    fd_score = None
    if args.mafft_msa and os.path.exists(args.mafft_msa):
        fd_score = compute_fd_score(args.true_msa, args.mafft_msa)

    # 4. Read SSS if provided
    sss_mean = None
    if args.sss_csv and os.path.exists(args.sss_csv):
        try:
            df_sss = pd.read_csv(args.sss_csv)
            if "sss_mean" in df_sss.columns:
                sss_mean = float(df_sss["sss_mean"].iloc[0])
            elif "mean_sss" in df_sss.columns:
                sss_mean = float(df_sss["mean_sss"].iloc[0])
        except Exception as e:
            sys.stderr.write(f"Warning reading SSS CSV: {e}\n")

    # 5. Calculate tree metrics
    rcrb_pwa_nj, nrf_pwa_nj = calculate_tree_metrics(args.true_tree, args.tree_pwa_nj) if args.tree_pwa_nj else (None, None)
    rcrb_msa_nj, nrf_msa_nj = calculate_tree_metrics(args.true_tree, args.tree_msa_nj) if args.tree_msa_nj else (None, None)
    rcrb_msa_ml, nrf_msa_ml = calculate_tree_metrics(args.true_tree, args.tree_msa_ml) if args.tree_msa_ml else (None, None)
    rcrb_true_nj, nrf_true_nj = calculate_tree_metrics(args.true_tree, args.tree_true_msa_nj) if args.tree_true_msa_nj else (None, None)
    rcrb_true_ml, nrf_true_ml = calculate_tree_metrics(args.true_tree, args.tree_true_msa_ml) if args.tree_true_msa_ml else (None, None)

    # Alignment Loss for ML
    aln_loss_ml = None
    if rcrb_true_ml is not None and rcrb_msa_ml is not None:
        aln_loss_ml = round(rcrb_true_ml - rcrb_msa_ml, 6)

    # Compile result row
    row = {
        "cond_name": args.cond_name,
        "simulator": args.simulator,
        "insert_rate": args.ins_rate,
        "delete_rate": args.del_rate,
        "distance": args.dist,
        "root_length": args.length,
        "replicate": args.rep,
        
        # True MSA properties
        "true_msa_length": t_props.get("aln_length"),
        "true_msa_expansion": round(t_props.get("aln_length", 0) / args.length, 3) if t_props.get("aln_length") else None,
        "true_gap_ratio": t_props.get("gap_ratio"),
        "true_gap_pct": t_props.get("gap_pct"),
        "mean_seq_len": t_props.get("mean_seq_len"),
        "std_seq_len": t_props.get("std_seq_len"),
        "min_seq_len": t_props.get("min_seq_len"),
        "max_seq_len": t_props.get("max_seq_len"),
        "true_num_indels": t_props.get("num_indel_blocks"),
        "true_mean_indel_len": t_props.get("mean_indel_block_len"),
        
        # MAFFT MSA properties
        "mafft_msa_length": m_props.get("aln_length"),
        "mafft_gap_ratio": m_props.get("gap_ratio"),
        "mafft_gap_pct": m_props.get("gap_pct"),
        "mafft_fd_score": round(fd_score, 6) if fd_score is not None else None,
        
        # Similarity
        "sss_mean": sss_mean,
        
        # Tree Metrics (RCRB = 1 - nRF)
        "rcrb_pwa_nj": rcrb_pwa_nj,
        "nrf_pwa_nj": nrf_pwa_nj,
        "rcrb_msa_nj": rcrb_msa_nj,
        "nrf_msa_nj": nrf_msa_nj,
        "rcrb_msa_ml": rcrb_msa_ml,
        "nrf_msa_ml": nrf_msa_ml,
        "rcrb_true_msa_nj": rcrb_true_nj,
        "nrf_true_msa_nj": nrf_true_nj,
        "rcrb_true_msa_ml": rcrb_true_ml,
        "nrf_true_msa_ml": nrf_true_ml,
        "alignment_loss_ml": aln_loss_ml,
    }

    df = pd.DataFrame([row])
    os.makedirs(os.path.dirname(os.path.abspath(args.outfile)), exist_ok=True)
    df.to_csv(args.outfile, index=False)
    print(f"Properties written to {args.outfile}")


if __name__ == "__main__":
    main()
