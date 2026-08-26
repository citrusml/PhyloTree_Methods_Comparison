#!/usr/bin/env python3
"""
Simulate Data Script
Generates true phylogenetic trees using IQ-TREE AliSim's RANDOM Birth-Death model
and simulates protein sequence alignments with substitutions (LG+G4 with customizable alpha)
and indels (customizable insert_rate and delete_rate).
"""

import sys
import os
import math
import random
import argparse
import subprocess
import shutil
from Bio import SeqIO
import dendropy

def generate_partition_nexus_file(filepath, length, ics_prop, base_model="LG+G4", alpha=1.0, seed=None):
    """
    Generates a NEXUS partition file for AliSim partitioning the sequence into
    regular LG and ICS partitions. Randomly samples round(length * ics_prop) sites as ICS.
    """
    rng = random.Random(seed)
    k_ics = int(round(length * ics_prop))
    k_ics = max(1, min(length - 1, k_ics)) if 0.0 < ics_prop < 1.0 else (length if ics_prop >= 1.0 else 0)

    all_sites = list(range(1, length + 1))
    ics_sites = sorted(rng.sample(all_sites, k_ics)) if k_ics > 0 else []
    lg_sites = sorted([s for s in all_sites if s not in set(ics_sites)])

    lg_model_str = format_model_string(base_model, alpha)
    # ICS model with corresponding Gamma rate heterogeneity
    ics_model_str = format_model_string("ICS+G4", alpha) if "+G" in base_model else "ICS"

    lines = ["#nexus", "begin sets;"]
    charpartitions = []

    if lg_sites:
        lines.append(f"    charset part_lg = {' '.join(map(str, lg_sites))};")
        charpartitions.append(f"{lg_model_str}:part_lg")
    if ics_sites:
        lines.append(f"    charset part_ics = {' '.join(map(str, ics_sites))};")
        charpartitions.append(f"{ics_model_str}:part_ics")

    lines.append(f"    charpartition mypart = {', '.join(charpartitions)};")
    lines.append("end;")
    lines.append("")

    with open(filepath, "w") as f:
        f.write("\n".join(lines))

def export_true_patristic_matrix(tree_file, out_matrix_file, ordered_labels=None):
    """
    Computes exact true patristic distance matrix from the true Newick tree using DendroPy
    and writes it as a PHYLIP distance matrix.
    """
    tns = dendropy.TaxonNamespace()
    tree = dendropy.Tree.get(path=tree_file, schema="newick", taxon_namespace=tns, preserve_underscores=True)
    pdm = tree.phylogenetic_distance_matrix()

    taxa_map = {t.label: t for t in tns}
    if ordered_labels:
        taxa_list = [taxa_map[lbl] for lbl in ordered_labels if lbl in taxa_map]
    else:
        taxa_list = sorted(list(tns), key=lambda x: x.label)

    N = len(taxa_list)
    with open(out_matrix_file, "w") as f:
        f.write(f"   {N}\n")
        for t1 in taxa_list:
            row_str = f"{t1.label:<10}"
            for t2 in taxa_list:
                d = pdm(t1, t2)
                row_str += f"  {d:.6f}"
            f.write(row_str + "\n")

def format_model_string(model, alpha=None):
    """
    Formats the substitution model string for AliSim.
    If the model string already specifies parameters (contains '{' and '}'), returns as is.
    If the model contains '+G' or '+G4' (Gamma rate heterogeneity) and alpha is provided, formats as 'MODEL+G4{alpha}'.
    If the model does not have '+G' (e.g., 'LG', 'WAG'), returns model string without alpha.
    """
    if not model:
        model = "LG+G4"

    # Already has parameter specification (e.g. "LG+G4{0.5}")
    if "{" in model and "}" in model:
        return model

    # If Gamma rate variation is in model, attach alpha if provided
    if "+G" in model:
        if alpha is not None:
            import re
            return re.sub(r'(\+G\d*)(\{[^}]*\})?', rf'\1{{{alpha}}}', model)
        return model
    else:
        # Non-gamma model (e.g., "LG", "WAG", "JC")
        return model

def find_iqtree_cmd():
    """Locates iqtree executable (iqtree3, iqtree2, or iqtree)."""
    candidates = ["iqtree3", "iqtree2", "iqtree"]
    # Check directory containing current python interpreter (e.g., conda / micromamba env)
    py_dir = os.path.dirname(sys.executable)
    for bin_name in candidates:
        full_path = os.path.join(py_dir, bin_name)
        if os.path.isfile(full_path) and os.access(full_path, os.X_OK):
            return full_path

    for bin_name in candidates:
        found = shutil.which(bin_name)
        if found:
            return found
        try:
            res = subprocess.run([bin_name, "--version"], capture_output=True, text=True)
            if res.returncode == 0:
                return bin_name
        except FileNotFoundError:
            pass
    return None

def simulate_tree_and_sequences_with_alisim(
    out_tree,
    out_fasta,
    out_true_msa=None,
    out_true_matrix=None,
    num_taxa=32,
    distance=1.0,
    length=300,
    birth_rate=0.1,
    death_rate=0.05,
    insert_rate=0.05,
    delete_rate=0.10,
    model="LG+G4",
    alpha=1.0,
    ics_prop=0.0,
    ics_model_file=None,
    seed=None
):
    """
    Simulates a random Birth-Death phylogenetic tree, sequences, and true MSA using IQ-TREE AliSim.
    Supports Invariant Category Sites (ICS) when ics_prop > 0.0 via AliSim partition models (-q)
    and pre-generated custom NEXUS model definition (--mdef).
    """
    iqtree_cmd = find_iqtree_cmd()
    if not iqtree_cmd:
        raise RuntimeError("Error: AliSim (iqtree3 / iqtree2 / iqtree) is not found in PATH. Please ensure IQ-TREE is installed and accessible.")

    prefix = out_fasta + ".alisim"
    tree_random_spec = f"RANDOM{{bd{{{birth_rate}/{death_rate}}}/{num_taxa}}}"

    if ics_prop > 0.0:
        # Locate pre-generated ICS model NEXUS file
        resolved_model_file = None
        candidates = [
            ics_model_file,
            "models/ics_model.nex",
            os.path.join(os.path.dirname(__file__), "../models/ics_model.nex"),
            os.path.join(os.getcwd(), "models/ics_model.nex")
        ]
        for cand in candidates:
            if cand and os.path.exists(cand) and os.path.getsize(cand) > 0:
                resolved_model_file = os.path.abspath(cand)
                break

        if not resolved_model_file:
            raise RuntimeError("Error: Invariant Category Sites (ICS) NEXUS model file not found. Please generate it using 'python3 bin/generate_ics_model.py' first.")

        k_ics = int(round(length * ics_prop))
        k_ics = max(1, min(length - 1, k_ics)) if 0.0 < ics_prop < 1.0 else (length if ics_prop >= 1.0 else 0)
        k_lg = length - k_ics

        if k_ics == length:
            # All sites are ICS
            model_ics_str = format_model_string("ICS+G4", alpha) if "+G" in model else "ICS"
            cmd = [
                iqtree_cmd,
                "--alisim", prefix,
                "--mdef", resolved_model_file,
                "-m", model_ics_str,
                "--length", str(length),
                "--seqtype", "AA",
                "-t", tree_random_spec,
                "--indel", f"{insert_rate},{delete_rate}",
                "--branch-scale", str(distance),
                "--redo"
            ]
            if seed is not None:
                cmd.extend(["-seed", str(seed)])
            res = subprocess.run(cmd, capture_output=True, text=True)
            if res.returncode != 0:
                raise RuntimeError(f"Error: AliSim execution failed with returncode {res.returncode}.\nCommand: {' '.join(cmd)}\nSTDOUT:\n{res.stdout}\nSTDERR:\n{res.stderr}")
        else:
            # 1. Simulate LG portion on random Birth-Death tree WITH indels
            prefix_lg = prefix + "_lg"
            model_lg_str = format_model_string(model, alpha)
            cmd_lg = [
                iqtree_cmd,
                "--alisim", prefix_lg,
                "-m", model_lg_str,
                "--length", str(k_lg),
                "-t", tree_random_spec,
                "--indel", f"{insert_rate},{delete_rate}",
                "--branch-scale", str(distance),
                "--redo"
            ]
            if seed is not None:
                cmd_lg.extend(["-seed", str(seed)])

            res_lg = subprocess.run(cmd_lg, capture_output=True, text=True)
            if res_lg.returncode != 0:
                raise RuntimeError(f"Error: AliSim LG simulation failed with returncode {res_lg.returncode}.\nCommand: {' '.join(cmd_lg)}\nSTDOUT:\n{res_lg.stdout}\nSTDERR:\n{res_lg.stderr}")

            alisim_treefile = prefix_lg + ".treefile"
            if not (os.path.exists(alisim_treefile) and os.path.getsize(alisim_treefile) > 0):
                alt_treefile = prefix_lg + ".tree"
                if os.path.exists(alt_treefile) and os.path.getsize(alt_treefile) > 0:
                    alisim_treefile = alt_treefile
                else:
                    raise RuntimeError(f"Error: AliSim did not produce expected output treefile: {alisim_treefile}")

            # 2. Simulate ICS portion on the EXACT SAME true tree WITHOUT indels (ICS sites are invariant to indels)
            prefix_ics = prefix + "_ics"
            model_ics_str = format_model_string("ICS+G4", alpha) if "+G" in model else "ICS"
            cmd_ics = [
                iqtree_cmd,
                "--alisim", prefix_ics,
                "--mdef", resolved_model_file,
                "-m", model_ics_str,
                "--length", str(k_ics),
                "-t", alisim_treefile,
                "--seqtype", "AA",
                "--redo"
            ]
            if seed is not None:
                cmd_ics.extend(["-seed", str(seed + 10007)])

            res_ics = subprocess.run(cmd_ics, capture_output=True, text=True)
            if res_ics.returncode != 0:
                raise RuntimeError(f"Error: AliSim ICS simulation failed with returncode {res_ics.returncode}.\nCommand: {' '.join(cmd_ics)}\nSTDOUT:\n{res_ics.stdout}\nSTDERR:\n{res_ics.stderr}")

            # 3. Write True Tree
            with open(alisim_treefile, "r") as f_in, open(out_tree, "w") as f_out:
                tree_content = f_in.read()
                f_out.write(tree_content)

            # 4. Helper to load alignment records from AliSim output (.phy or .fa)
            def load_alignment_records(pref):
                fa_f = pref + ".fa"
                phy_f = pref + ".phy"
                if os.path.exists(fa_f) and os.path.getsize(fa_f) > 0:
                    return {rec.id.strip("_"): str(rec.seq) for rec in SeqIO.parse(fa_f, "fasta")}
                elif os.path.exists(phy_f) and os.path.getsize(phy_f) > 0:
                    return {rec.id.strip("_"): str(rec.seq) for rec in SeqIO.parse(phy_f, "phylip")}
                else:
                    raise RuntimeError(f"Error: Could not find alignment output for prefix {pref}")

            lg_aln = load_alignment_records(prefix_lg)
            ics_aln = load_alignment_records(prefix_ics)

            taxa_keys = list(lg_aln.keys())
            if len(taxa_keys) < num_taxa:
                raise RuntimeError(f"Error: AliSim generated only {len(taxa_keys)} taxa, expected {num_taxa}.")

            lg_cols = len(lg_aln[taxa_keys[0]])

            # 5. Interleave ICS sites evenly across the sequence length
            # Divides the LG alignment into (k_ics + 1) chunks and places 1 indel-free ICS site between chunks
            chunk_step = lg_cols / (k_ics + 1)
            split_indices = [int(round(i * chunk_step)) for i in range(k_ics + 2)]

            combined_msa = {}
            for tid in taxa_keys:
                seq_lg = lg_aln[tid]
                seq_ics = ics_aln[tid]
                chunks = []
                for i in range(k_ics):
                    chunks.append(seq_lg[split_indices[i]:split_indices[i+1]])
                    chunks.append(seq_ics[i])
                chunks.append(seq_lg[split_indices[k_ics]:split_indices[k_ics+1]])
                combined_msa[tid] = "".join(chunks)

            # Write combined unaligned FASTA (gaps removed)
            with open(out_fasta, "w") as out_f:
                for tid in taxa_keys:
                    seq_clean = combined_msa[tid].replace("-", "")
                    if len(seq_clean) < 1:
                        raise RuntimeError(f"Error: Taxon {tid} has empty combined sequence.")
                    out_f.write(f">{tid}\n{seq_clean}\n")

            # Write combined True MSA if requested
            if out_true_msa:
                with open(out_true_msa, "w") as out_m:
                    for tid in taxa_keys:
                        out_m.write(f">{tid}\n{combined_msa[tid]}\n")

            # Write True Patristic Distance Matrix if requested
            if out_true_matrix:
                export_true_patristic_matrix(out_tree, out_true_matrix, ordered_labels=taxa_keys)

            # Clean up all temporary files from both runs
            for p in [prefix_lg, prefix_ics]:
                for ext in [".phy", ".unaligned.fa", ".fa", ".treefile", ".tree", ".log", ".iqtree"]:
                    fpath = p + ext
                    if os.path.exists(fpath):
                        try:
                            os.remove(fpath)
                        except OSError:
                            pass
            return
    else:
        model_str = format_model_string(model, alpha)
        cmd = [
            iqtree_cmd,
            "--alisim", prefix,
            "-m", model_str,
            "--length", str(length),
            "-t", tree_random_spec,
            "--indel", f"{insert_rate},{delete_rate}",
            "--branch-scale", str(distance),
            "--redo"
        ]
        if seed is not None:
            cmd.extend(["-seed", str(seed)])

        res = subprocess.run(cmd, capture_output=True, text=True)
        if res.returncode != 0:
            raise RuntimeError(f"Error: AliSim execution failed with returncode {res.returncode}.\nCommand: {' '.join(cmd)}\nSTDOUT:\n{res.stdout}\nSTDERR:\n{res.stderr}")

    # 1. Process True Tree file
    alisim_treefile = prefix + ".treefile"
    if not (os.path.exists(alisim_treefile) and os.path.getsize(alisim_treefile) > 0):
        alt_treefile = prefix + ".tree"
        if os.path.exists(alt_treefile) and os.path.getsize(alt_treefile) > 0:
            alisim_treefile = alt_treefile
        else:
            raise RuntimeError(f"Error: AliSim did not produce expected output treefile: {alisim_treefile}\nCommand: {' '.join(cmd)}\nSTDOUT:\n{res.stdout}\nSTDERR:\n{res.stderr}")

    with open(alisim_treefile, "r") as f_in, open(out_tree, "w") as f_out:
        tree_content = f_in.read()
        f_out.write(tree_content)

    # 2. Process Sequences and True Alignment files
    # AliSim outputs .fa or .phy for True MSA, and .unaligned.fa for unaligned sequences
    alisim_aligned_fa = prefix + ".fa"
    alisim_unaligned_fa = prefix + ".unaligned.fa"
    alisim_phy = prefix + ".phy"

    aligned_records = []
    if os.path.exists(alisim_aligned_fa) and os.path.getsize(alisim_aligned_fa) > 0:
        aligned_records = list(SeqIO.parse(alisim_aligned_fa, "fasta"))
    elif os.path.exists(alisim_phy) and os.path.getsize(alisim_phy) > 0:
        aligned_records = list(SeqIO.parse(alisim_phy, "phylip"))

    if os.path.exists(alisim_unaligned_fa) and os.path.getsize(alisim_unaligned_fa) > 0:
        unaligned_records = list(SeqIO.parse(alisim_unaligned_fa, "fasta"))
    elif aligned_records:
        unaligned_records = aligned_records
    else:
        raise RuntimeError(f"Error: AliSim did not produce expected sequence files.\nCommand: {' '.join(cmd)}\nSTDOUT:\n{res.stdout}\nSTDERR:\n{res.stderr}")

    # Check taxa count
    min_required_len = 1
    if len(unaligned_records) < num_taxa:
        raise RuntimeError(f"Error: AliSim generated only {len(unaligned_records)} taxa, but expected {num_taxa}.")

    # Write unaligned FASTA (gaps stripped)
    with open(out_fasta, "w") as out_f:
        for rec in unaligned_records:
            clean_id = rec.id.strip("_")
            clean_seq = str(rec.seq).replace("-", "")
            if len(clean_seq) < min_required_len:
                raise RuntimeError(f"Error: Taxon {rec.id} has empty sequence in AliSim output.")
            out_f.write(f">{clean_id}\n{clean_seq}\n")

    # Write True MSA FASTA if requested
    if out_true_msa:
        with open(out_true_msa, "w") as out_m:
            if aligned_records:
                for rec in aligned_records:
                    clean_id = rec.id.strip("_")
                    out_m.write(f">{clean_id}\n{str(rec.seq)}\n")
            else:
                # If aligned FASTA wasn't present, write unaligned sequences as fallback
                for rec in unaligned_records:
                    clean_id = rec.id.strip("_")
                    out_m.write(f">{clean_id}\n{str(rec.seq)}\n")

    # Write True Patristic Distance Matrix if requested
    if out_true_matrix:
        taxa_labels = [rec.id.strip("_") for rec in unaligned_records]
        export_true_patristic_matrix(out_tree, out_true_matrix, ordered_labels=taxa_labels)

    # Clean up temporary AliSim generated and NEXUS files
    for ext in [".phy", ".unaligned.fa", ".fa", ".treefile", ".tree", ".log", ".iqtree", ".partition.nex"]:
        fpath = prefix + ext
        if os.path.exists(fpath):
            try:
                os.remove(fpath)
            except OSError:
                pass

def main():
    parser = argparse.ArgumentParser(description="Simulate True Tree and Protein Sequences with Indels using AliSim")
    parser.add_argument("--distance", type=float, default=1.0, help="Branch length scale / Evolutionary distance D (default: 1.0)")
    parser.add_argument("--length", type=int, default=300, help="Sequence length L (aa) (default: 300)")
    parser.add_argument("--num_taxa", type=int, default=32, help="Taxon count N (default: 32)")
    parser.add_argument("--birth_rate", type=float, default=0.1, help="Birth rate for Birth-Death tree model (default: 0.1)")
    parser.add_argument("--death_rate", type=float, default=0.05, help="Death rate for Birth-Death tree model (default: 0.05)")
    parser.add_argument("--insert_rate", type=float, default=0.05, help="Insertion rate relative to substitutions (default: 0.05)")
    parser.add_argument("--delete_rate", type=float, default=0.10, help="Deletion rate relative to substitutions (default: 0.10)")
    parser.add_argument("--indel_rate", type=float, default=None, help="Symmetric indel rate for both insertion and deletion (overrides insert/delete if specified)")
    parser.add_argument("--alpha", type=float, default=1.0, help="Gamma shape parameter alpha for LG+G4 (default: 1.0)")
    parser.add_argument("--model", type=str, default="LG+G4", help="Substitution model prefix (default: LG+G4)")
    parser.add_argument("--ics_prop", type=float, default=0.0, help="Proportion of Invariant Category Sites (ICS) under Dayhoff 6 classes (default: 0.0)")
    parser.add_argument("--ics_model_file", help="Path to pre-generated Invariant Category Sites (ICS) NEXUS model file (default: models/ics_model.nex)")
    parser.add_argument("--outtree", required=True, help="Output true Newick tree file")
    parser.add_argument("--outfasta", required=True, help="Output unaligned FASTA file")
    parser.add_argument("--outtrue_msa", help="Output True MSA (aligned) FASTA file")
    parser.add_argument("--outtrue_matrix", help="Output True Patristic Distance Matrix (PHYLIP format)")
    parser.add_argument("--seed", type=int, default=None, help="Random seed for reproducibility")
    args = parser.parse_args()

    insert_rate = args.insert_rate if args.indel_rate is None else args.indel_rate
    delete_rate = args.delete_rate if args.indel_rate is None else args.indel_rate

    base_seed = args.seed if args.seed is not None else random.randint(1, 1000000)
    max_trials = 10
    last_error = None

    for trial in range(max_trials):
        current_seed = base_seed + trial * 10007
        random.seed(current_seed)

        try:
            simulate_tree_and_sequences_with_alisim(
                out_tree=args.outtree,
                out_fasta=args.outfasta,
                out_true_msa=args.outtrue_msa,
                out_true_matrix=args.outtrue_matrix,
                num_taxa=args.num_taxa,
                distance=args.distance,
                length=args.length,
                birth_rate=args.birth_rate,
                death_rate=args.death_rate,
                insert_rate=insert_rate,
                delete_rate=delete_rate,
                model=args.model,
                alpha=args.alpha,
                ics_prop=args.ics_prop,
                ics_model_file=args.ics_model_file,
                seed=current_seed
            )
            used_model_str = f"{args.model}(ics_prop={args.ics_prop})" if args.ics_prop > 0 else format_model_string(args.model, args.alpha)
            print(f"Simulated data generated (trial {trial+1}/{max_trials}): Tree -> {args.outtree}, FASTA -> {args.outfasta} (N={args.num_taxa}, D={args.distance}, L={args.length}, bd={{{args.birth_rate}/{args.death_rate}}}, model={used_model_str}, indel={{{insert_rate},{delete_rate}}})")
            return
        except Exception as e:
            last_error = e
            # Clean up partial output if any
            for p in [args.outfasta, args.outtree, args.outtrue_msa, args.outtrue_matrix]:
                if p and os.path.exists(p):
                    try:
                        os.remove(p)
                    except OSError:
                        pass

    # If all trials failed, raise the last exception
    raise RuntimeError(f"Failed to generate valid sequence simulation after {max_trials} attempts. Last error: {last_error}")

if __name__ == "__main__":
    main()

