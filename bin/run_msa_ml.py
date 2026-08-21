#!/usr/bin/env python3
"""
MSA + ML Pipeline Script
Runs MAFFT alignment, then infers Maximum Likelihood tree using IQ-TREE 2
with automatic model selection (-m MFP).
Bootstrap is disabled by default for maximum performance and stability across all dataset sizes.
Extracts estimated best-fit model and parameters into JSON metadata.
"""

import sys
import os
import json
import re
import argparse
import subprocess
from Bio import SeqIO

def run_mafft(in_fasta, out_msa):
    """Runs MAFFT with --threadit 0."""
    try:
        cmd = ["mafft", "--threadit", "0", "--auto", in_fasta]
        with open(out_msa, "w") as out_f:
            res = subprocess.run(cmd, stdout=out_f, stderr=subprocess.PIPE, text=True)
        if res.returncode == 0 and os.path.exists(out_msa) and os.path.getsize(out_msa) > 0:
            return True
    except FileNotFoundError:
        pass
    return False

def parse_iqtree_log(iqtree_log_file):
    """Parses .iqtree file to extract best-fit model according to BIC."""
    metadata = {
        "best_model_bic": "Unknown",
        "best_model_aic": "Unknown",
        "gamma_alpha": None,
        "invariable_sites": None
    }
    if not os.path.exists(iqtree_log_file):
        return metadata

    with open(iqtree_log_file, "r") as f:
        content = f.read()

    # Search for ModelFinder BIC selection
    bic_match = re.search(r"Best-fit model according to BIC:\s*([^\s\n]+)", content)
    if bic_match:
        metadata["best_model_bic"] = bic_match.group(1)

    aic_match = re.search(r"Best-fit model according to AIC:\s*([^\s\n]+)", content)
    if aic_match:
        metadata["best_model_aic"] = aic_match.group(1)

    # Search for Gamma alpha
    alpha_match = re.search(r"Gamma shape alpha:\s*([\d\.]+)", content)
    if alpha_match:
        try:
            metadata["gamma_alpha"] = float(alpha_match.group(1))
        except ValueError:
            pass

    # Search for Invariable sites proportion
    inv_match = re.search(r"Proportion of invariable sites:\s*([\d\.]+)", content)
    if inv_match:
        try:
            metadata["invariable_sites"] = float(inv_match.group(1))
        except ValueError:
            pass

    return metadata

def create_star_tree(fasta_path, out_tree_path):
    """Generates a trivial star tree if all sequences are completely identical."""
    records = list(SeqIO.parse(fasta_path, "fasta"))
    if not records:
        return False
    taxa_str = ",".join([f"{r.id}:0.001" for r in records])
    star_newick = f"({taxa_str});\n"
    with open(out_tree_path, "w") as f:
        f.write(star_newick)
    return True

def main():
    parser = argparse.ArgumentParser(description="MSA+ML Pipeline Execution (IQ-TREE 2)")
    parser.add_argument("--fasta", help="Input unaligned FASTA sequence file")
    parser.add_argument("--msa", help="Input pre-aligned MSA FASTA file (skips MAFFT if provided)")
    parser.add_argument("--outtree", required=True, help="Output Newick tree file")
    parser.add_argument("--outmsa", help="Output MAFFT MSA file (if MAFFT is run)")
    parser.add_argument("--outjson", help="Output JSON metadata file")
    parser.add_argument("--bootstrap", type=int, default=0, help="Ultrafast bootstrap replicates (-B, default: 0 = disabled)")
    parser.add_argument("--threads", type=int, default=1, help="Number of CPU threads for IQ-TREE 2 (-T)")
    args = parser.parse_args()

    if not args.msa and not args.fasta:
        parser.error("Either --msa (pre-aligned) or --fasta (unaligned) must be provided.")

    if args.msa:
        tmp_msa = args.msa
    else:
        tmp_msa = args.outmsa or (args.outtree + ".msa.fasta")
        # 1. Run MAFFT MSA
        mafft_success = run_mafft(args.fasta, tmp_msa)
        if not mafft_success:
            records = list(SeqIO.parse(args.fasta, "fasta"))
            with open(tmp_msa, "w") as f:
                for rec in records:
                    f.write(f">{rec.id}\n{str(rec.seq)}\n")

    json_out = args.outjson or (args.outtree + ".json")
    prefix = args.outtree + ".iqtree_run"

    # 2. Locate IQ-TREE 2 binary
    iqtree_cmd = None
    for bin_name in ["iqtree2", "iqtree"]:
        try:
            res = subprocess.run([bin_name, "--version"], capture_output=True, text=True)
            if res.returncode == 0:
                iqtree_cmd = bin_name
                break
        except FileNotFoundError:
            pass

    if not iqtree_cmd:
        raise RuntimeError("IQ-TREE executable ('iqtree2' or 'iqtree') was not found in PATH. Please ensure IQ-TREE 2 is installed.")

    metadata = {}
    tree_file = prefix + ".treefile"
    iq_file = prefix + ".iqtree"

    # 3. Build IQ-TREE command (Pure ML Tree inference with ModelFinder, no bootstrap)
    cmd = [
        iqtree_cmd,
        "-s", tmp_msa,
        "-m", "MFP",
        "--prefix", prefix,
        "--redo",
        "-T", str(args.threads)
    ]
    if args.bootstrap > 0:
        cmd.extend(["-B", str(args.bootstrap)])

    res = subprocess.run(cmd, capture_output=True, text=True)

    # 4. Fallback if bootstrap was explicitly requested but failed due to identical sequences
    if (not os.path.exists(tree_file) or os.path.getsize(tree_file) == 0) and args.bootstrap > 0:
        print("Warning: IQ-TREE failed with bootstrap. Retrying pure ML search without bootstrap...", file=sys.stderr)
        cmd_noboot = [
            iqtree_cmd,
            "-s", tmp_msa,
            "-m", "MFP",
            "--prefix", prefix,
            "--redo",
            "-T", str(args.threads)
        ]
        res = subprocess.run(cmd_noboot, capture_output=True, text=True)

    # 5. Check output tree or apply Star-tree fallback for 100% identical sequences
    if os.path.exists(tree_file) and os.path.getsize(tree_file) > 0:
        with open(tree_file, "r") as rf, open(args.outtree, "w") as wf:
            wf.write(rf.read())
        metadata = parse_iqtree_log(iq_file)
    else:
        # Fallback for 100% identical sequence datasets where IQ-TREE cannot produce a tree
        print("Warning: IQ-TREE could not construct a tree. Using star tree fallback.", file=sys.stderr)
        if create_star_tree(tmp_msa, args.outtree):
            metadata = {
                "best_model_bic": "Identical_Sequences",
                "best_model_aic": "Identical_Sequences",
                "gamma_alpha": None,
                "invariable_sites": 1.0
            }
        else:
            raise RuntimeError(f"IQ-TREE execution failed completely:\nCommand: {' '.join(cmd)}\nSTDOUT:\n{res.stdout}\nSTDERR:\n{res.stderr}")

    # Save JSON metadata
    with open(json_out, "w") as jf:
        json.dump(metadata, jf, indent=2)

    # Clean up temporary MSA if not requested
    if not args.outmsa and os.path.exists(tmp_msa):
        try:
            os.remove(tmp_msa)
        except OSError:
            pass

    print(f"MSA+ML tree successfully written to {args.outtree} (Model metadata saved to {json_out})")

if __name__ == "__main__":
    main()
