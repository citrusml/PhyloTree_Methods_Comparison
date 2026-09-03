#!/usr/bin/env python3
"""
Run MSA+BI (Bayesian Inference with MrBayes v3.2.6+)
Performs Bayesian phylogenetic tree reconstruction on aligned protein sequences (MSA)
using MrBayes, following the benchmark protocol of Matsui & Iwasaki (2020, Systematic Biology).

MrBayes settings:
  lset rates=gamma;
  prset aamodelpr=fixed(lg);
  mcmcp ngen=100000 samplefreq=1000 nruns=1 nchains=4;
  sumt burnin=20 contype=halfcompat;
"""

import os
import sys
import shutil
import argparse
import subprocess
import time
import json
from typing import List
from Bio import SeqIO
import dendropy

def convert_fasta_to_mrbayes_nexus(
    fasta_path: str,
    nexus_path: str,
    model: str = "lg",
    rates: str = "gamma",
    ngen: int = 100000,
    samplefreq: int = 1000,
    nruns: int = 1,
    nchains: int = 4,
    burnin: int = 20,
    contype: str = "halfcompat"
) -> None:
    """Convert an aligned FASTA file into a NEXUS file with an embedded MrBayes execution block."""
    records = list(SeqIO.parse(fasta_path, "fasta"))
    if not records:
        raise ValueError(f"No sequences found in FASTA: {fasta_path}")
    
    ntax = len(records)
    nchar = len(records[0].seq)
    
    with open(nexus_path, "w") as f:
        f.write("#NEXUS\n\n")
        f.write("begin data;\n")
        f.write(f"    dimensions ntax={ntax} nchar={nchar};\n")
        f.write("    format datatype=protein gap=- missing=?;\n")
        f.write("    matrix\n")
        for r in records:
            seq_str = str(r.seq).strip()
            f.write(f"    {r.id}    {seq_str}\n")
        f.write("    ;\nend;\n\n")
        
        f.write("begin mrbayes;\n")
        f.write("    set autoclose=yes nowarn=yes;\n")
        f.write(f"    lset rates={rates};\n")
        f.write(f"    prset aamodelpr=fixed({model.lower()});\n")
        f.write(f"    mcmcp ngen={ngen} samplefreq={samplefreq} nruns={nruns} nchains={nchains};\n")
        f.write("    mcmc;\n")
        f.write(f"    sumt burnin={burnin} contype={contype};\n")
        f.write("end;\n")

def find_mrbayes_binary() -> str:
    """Locate the MrBayes executable ('mb' or 'mrbayes')."""
    for candidate in ["mb", "mrbayes"]:
        path = shutil.which(candidate)
        if path:
            return path
    # Search common micromamba/conda environments
    conda_prefix = os.environ.get("CONDA_PREFIX")
    if conda_prefix:
        mb_path = os.path.join(conda_prefix, "bin", "mb")
        if os.path.exists(mb_path):
            return mb_path
    home = os.path.expanduser("~")
    for env_name in ["phylomethod_env", "base"]:
        test_path = os.path.join(home, ".micromamba", "envs", env_name, "bin", "mb")
        if os.path.exists(test_path):
            return test_path
    raise FileNotFoundError("MrBayes executable ('mb') not found in PATH or conda environments.")

def extract_consensus_tree(con_tre_path: str, outtree_path: str) -> None:
    """Read MrBayes consensus tree (.con.tre) and export it as a clean Newick tree."""
    if not os.path.exists(con_tre_path):
        raise FileNotFoundError(f"Consensus tree file not generated: {con_tre_path}")
    
    tns = dendropy.TaxonNamespace()
    trees = dendropy.TreeList.get(path=con_tre_path, schema="nexus", taxon_namespace=tns, preserve_underscores=True)
    if not trees:
        raise ValueError(f"No trees parsed from {con_tre_path}")
    
    # In MrBayes con.tre, the majority-rule consensus tree is the primary tree
    consensus_tree = trees[0]
    
    # Clean taxon labels and strip square bracket annotations if necessary
    newick_str = consensus_tree.as_string(schema="newick").strip()
    
    # Remove leading [&U] if present to ensure maximum compatibility
    if newick_str.startswith("[&U]"):
        newick_str = newick_str[4:].strip()
    
    with open(outtree_path, "w") as f:
        f.write(newick_str + "\n")

def main():
    parser = argparse.ArgumentParser(description="Run MSA + MrBayes Bayesian Inference (syz049 benchmark)")
    parser.add_argument("--msa", required=True, help="Input aligned FASTA file")
    parser.add_argument("--outtree", required=True, help="Output Newick tree file (.nwk)")
    parser.add_argument("--outjson", help="Output metadata JSON file")
    parser.add_argument("--model", default="lg", help="Amino acid substitution model (default: lg)")
    parser.add_argument("--rates", default="gamma", help="Rate variation across sites (default: gamma)")
    parser.add_argument("--ngen", type=int, default=100000, help="Number of MCMC generations (default: 100000)")
    parser.add_argument("--samplefreq", type=int, default=1000, help="Sampling frequency (default: 1000)")
    parser.add_argument("--burnin", type=int, default=20, help="Burn-in samples to discard (default: 20)")
    parser.add_argument("--nchains", type=int, default=4, help="Number of MCMC chains (default: 4)")
    parser.add_argument("--nruns", type=int, default=1, help="Number of independent runs (default: 1)")
    parser.add_argument("--contype", default="halfcompat", help="Consensus tree type (default: halfcompat)")
    parser.add_argument("--threads", type=int, default=2, help="Number of threads / cores (default: 2)")
    args = parser.parse_args()

    mb_bin = find_mrbayes_binary()
    
    # Work directory setup (use directory of outtree or local dir)
    outdir = os.path.dirname(os.path.abspath(args.outtree)) or "."
    base_name = os.path.splitext(os.path.basename(args.msa))[0]
    nexus_file = os.path.join(outdir, f"{base_name}_mb.nex")
    log_file = os.path.join(outdir, f"{base_name}_mb.log")
    con_tre_file = f"{nexus_file}.con.tre"

    t0 = time.time()
    
    # 1. Convert FASTA to NEXUS with embedded MrBayes block
    convert_fasta_to_mrbayes_nexus(
        fasta_path=args.msa,
        nexus_path=nexus_file,
        model=args.model,
        rates=args.rates,
        ngen=args.ngen,
        samplefreq=args.samplefreq,
        nruns=args.nruns,
        nchains=args.nchains,
        burnin=args.burnin,
        contype=args.contype
    )

    # 2. Run MrBayes
    cmd = [mb_bin, os.path.basename(nexus_file)]
    with open(log_file, "w") as lf:
        proc = subprocess.run(cmd, cwd=outdir, stdout=lf, stderr=subprocess.STDOUT)
    
    if proc.returncode != 0:
        with open(log_file, "r") as lf:
            tail = "".join(lf.readlines()[-30:])
        raise RuntimeError(f"MrBayes execution failed (exit code {proc.returncode}). Log tail:\n{tail}")

    # 3. Extract consensus tree to Newick
    extract_consensus_tree(con_tre_file, args.outtree)
    elapsed = time.time() - t0

    # 4. Save metadata JSON if requested
    if args.outjson:
        metadata = {
            "pipeline": "MSA+BI",
            "model": args.model,
            "rates": args.rates,
            "ngen": args.ngen,
            "samplefreq": args.samplefreq,
            "burnin": args.burnin,
            "nchains": args.nchains,
            "elapsed_sec": round(elapsed, 2)
        }
        with open(args.outjson, "w") as jf:
            json.dump(metadata, jf, indent=2)

    print(f"[MSA+BI] Successfully completed MrBayes inference in {elapsed:.1f}s -> {args.outtree}")

if __name__ == "__main__":
    main()
