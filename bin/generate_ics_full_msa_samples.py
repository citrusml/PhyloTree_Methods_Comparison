#!/usr/bin/env python3
"""
Generate True Reference and Estimated MSAs for Exp18 (ICS Full Benchmark).

Conditions:
- Taxa: 32
- Distances (D): [0.1, 0.2, 0.5, 1.0, 1.5, 2.0]
- Lengths (L): [300, 500, 1000, 1500]
- Tree model: paper_yule (Matsui & Iwasaki 2020)
- Substitution model: 100% Invariant Category Sites (ICS+G4{1.0})
- Indel: Zipfian POW{1.7/50} (insert=0.10, delete=0.10)
- Estimated MSA tool: MAFFT (--auto)
"""

from typing import List, Dict, Any
import os
import sys
import shutil
import itertools
from pathlib import Path
import pandas as pd
from Bio import SeqIO

# Ensure local imports work
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from src.ics_simulator import generate_paper_tree, simulate_ics_full_sequences, align_with_mafft
from bin.calculate_alignment_metrics import compute_sp_score


def main() -> None:
    project_root = Path(__file__).resolve().parent.parent
    results_dir = project_root / "results" / "results_ics_full"
    sample_dir = results_dir / "sample_alignments"
    rep_dir = results_dir / "replications"
    sample_dir.mkdir(parents=True, exist_ok=True)
    rep_dir.mkdir(parents=True, exist_ok=True)

    mafft_bin = "/Users/aibarakazuki/.micromamba/envs/phylomethod_env/bin/mafft"
    iqtree_bin = str(project_root / "bin" / "iqtree2")
    model_nex = str(project_root / "models" / "ics_model.nex")

    distances: List[float] = [0.1, 0.2, 0.5, 1.0, 1.5, 2.0]
    lengths: List[int] = [300, 500, 1000, 1500]
    taxa: int = 32

    records: List[Dict[str, Any]] = []
    total = len(distances) * len(lengths)
    print(f"=== Starting MSA Generation for Exp18 ICS Full (Total conditions: {total}) ===")

    for idx, (dist, length) in enumerate(itertools.product(distances, lengths), 1):
        cond_name = f"D{dist}_L{length}"
        cond_dir = sample_dir / cond_name
        cond_rep_dir = rep_dir / f"{cond_name}_rep0"
        cond_dir.mkdir(parents=True, exist_ok=True)
        cond_rep_dir.mkdir(parents=True, exist_ok=True)

        seed = 1000 + idx * 77

        # File targets
        tree_path = str(cond_dir / "true_tree.nwk")
        sim_prefix = str(cond_dir / "sim")
        true_msa_dest = cond_dir / "true_msa.fasta"
        est_msa_dest = cond_dir / "estimated_msa.fasta"
        unaligned_dest = cond_dir / "seqs.fasta"

        # 1. Tree generation
        generate_paper_tree(
            taxa=taxa,
            scale=dist,
            seed=seed,
            outtree_path=tree_path
        )

        # 2. AliSim simulation
        unaligned_fa, true_msa_fa = simulate_ics_full_sequences(
            tree_path=tree_path,
            length=length,
            out_prefix=sim_prefix,
            seed=seed,
            model_nex_path=model_nex,
            iqtree_bin=iqtree_bin
        )

        # Rename to clean filenames in cond_dir
        if os.path.exists(unaligned_dest):
            os.remove(unaligned_dest)
        if os.path.exists(true_msa_dest):
            os.remove(true_msa_dest)
        os.rename(unaligned_fa, unaligned_dest)
        os.rename(true_msa_fa, true_msa_dest)

        # Clean intermediate alisim log / tree files
        for extra in cond_dir.glob("sim*"):
            try:
                os.remove(extra)
            except OSError:
                pass

        # 3. MAFFT alignment
        align_with_mafft(
            unaligned_fasta=str(unaligned_dest),
            out_msa_path=str(est_msa_dest),
            mafft_bin=mafft_bin,
            threads=2
        )

        # 4. Copy to standard replications directory (Nextflow compatible: rep0/true_msa.fasta, rep0/msa.fasta)
        shutil.copy2(true_msa_dest, cond_rep_dir / "true_msa.fasta")
        shutil.copy2(est_msa_dest, cond_rep_dir / "msa.fasta")
        shutil.copy2(tree_path, cond_rep_dir / "true_tree.nwk")
        shutil.copy2(unaligned_dest, cond_rep_dir / "seqs.fasta")

        # Also create flat convenience copies in sample_alignments: D{D}_L{L}_true_msa.fasta, etc.
        shutil.copy2(true_msa_dest, sample_dir / f"{cond_name}_true_msa.fasta")
        shutil.copy2(est_msa_dest, sample_dir / f"{cond_name}_estimated_msa.fasta")

        # 5. Metrics calculation
        sp_score, modeler_score, correct_pairs, total_true, total_est = compute_sp_score(
            str(true_msa_dest), str(est_msa_dest)
        )

        t_records = list(SeqIO.parse(str(true_msa_dest), "fasta"))
        e_records = list(SeqIO.parse(str(est_msa_dest), "fasta"))
        t_aln_len = len(t_records[0].seq) if t_records else 0
        e_aln_len = len(e_records[0].seq) if e_records else 0

        record = {
            "distance": dist,
            "length": length,
            "condition": cond_name,
            "seed": seed,
            "sp_score": round(sp_score, 4) if sp_score is not None else None,
            "modeler_score": round(modeler_score, 4) if modeler_score is not None else None,
            "true_msa_length": t_aln_len,
            "est_msa_length": e_aln_len,
            "true_msa_path": str(true_msa_dest.relative_to(project_root)),
            "est_msa_path": str(est_msa_dest.relative_to(project_root)),
            "replications_dir": str(cond_rep_dir.relative_to(project_root))
        }
        records.append(record)
        print(f"[{idx:02d}/{total:02d}] {cond_name}: SP Score = {record['sp_score']:.4f}, TrueLen = {t_aln_len}, EstLen = {e_aln_len}")

    df_summary = pd.DataFrame(records)
    csv_out = sample_dir / "alignment_samples_summary.csv"
    df_summary.to_csv(csv_out, index=False)
    # Also save to results_ics_full root
    df_summary.to_csv(results_dir / "sample_alignments_summary.csv", index=False)

    print(f"\n✓ Successfully generated all {total} true and estimated MSAs!")
    print(f"✓ Summary table saved to: {csv_out}")


if __name__ == "__main__":
    main()
