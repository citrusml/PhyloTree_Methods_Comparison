#!/usr/bin/env python3
"""
analyze_dataset_similarity.py
=============================
Analyzes Sequence Similarity Score (SSS) and Pairwise Identity for HOGs
defined in the dataset/ directory (Amniota, Eukaryota, Bacteria, Actinobacteria).

Default outputs:
- dataset/results/hog_similarity_summary.csv
- dataset/results/pairwise_similarity_details.csv
- dataset/results/dataset_level_summary.csv
"""

import os
import sys
import argparse
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, List, Tuple, Any, Optional

import numpy as np
import pandas as pd
from Bio import SeqIO

# Ensure script dir and project root are in sys.path
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent.parent
DATASET_DIR = SCRIPT_DIR.parent

for p in [str(SCRIPT_DIR), str(PROJECT_ROOT)]:
    if p not in sys.path:
        sys.path.insert(0, p)

try:
    from similarity_calculator import compute_pairwise_sss_and_identity, get_pairwise_aligner
except ImportError:
    from src.similarity_calculator import compute_pairwise_sss_and_identity, get_pairwise_aligner

DEFAULT_DATASET_DIR = DATASET_DIR
DEFAULT_OUT_DIR = DATASET_DIR / "results"

# Known local search locations for sequences
LOCAL_PATHS = {
    "amn": [
        Path("/Users/aibarakazuki/amplify/input_fasta/amniota"),
        Path("/Users/aibarakazuki/amplify/results/amniota/amn3/msa"),
        Path("/Users/aibarakazuki/amplify/archive/input_fasta_new/amniota"),
        PROJECT_ROOT / "data" / "amniota"
    ],
    "euk": [
        Path("/Users/aibarakazuki/amplify/archive/input_fasta_new/eukaryota"),
        Path("/Users/aibarakazuki/amplify/results/eukaryota_new/euk3/msa"),
        Path("/Users/aibarakazuki/amplify/input_fasta/eukaryota"),
        PROJECT_ROOT / "data" / "eukaryota"
    ],
    "bac": [
        Path("/Users/aibarakazuki/amplify/input_fasta/bacteria"),
        PROJECT_ROOT / "data" / "bacteria",
        Path("/Users/aibarakazuki/amplify/results/bacteria/bac3/msa"),
        Path("/Users/aibarakazuki/input_fasta/bacteria")
    ],
    "act": [
        Path("/Users/aibarakazuki/amplify/input_fasta/actinomycete"),
        PROJECT_ROOT / "data" / "actinomycete",
        Path("/Users/aibarakazuki/amplify/results/actinomycete/act3/msa"),
        Path("/Users/aibarakazuki/input_fasta/actinomycete")
    ]
}


def load_hogs_list(hogs_file: Path) -> List[str]:
    """Loads HOG identifiers from a TSV/text file."""
    hogs = []
    with open(hogs_file, "r", encoding="utf-8") as f:
        for line in f:
            h = line.strip().split("\t")[0].strip()
            if h:
                norm_h = h.replace("HOG:", "").replace("HOG_", "")
                hogs.append(norm_h)
    return hogs


def load_taxa_list(taxa_file: Path) -> Dict[str, str]:
    """Loads taxon mnemonic codes and their lineage descriptions from *_taxa.tsv."""
    taxa = {}
    with open(taxa_file, "r", encoding="utf-8") as f:
        for line in f:
            parts = line.strip().split("\t")
            if parts and parts[0]:
                code = parts[0].strip()
                desc = parts[2].strip() if len(parts) > 2 else ""
                taxa[code] = desc
    return taxa


def find_hog_fasta_file(dataset_key: str, hog_id: str, custom_dirs: Optional[List[Path]] = None) -> Optional[Path]:
    """
    Finds existing FASTA or MSA file for a given HOG ID.
    Checks candidate names: '{hog_id}.msa.fa', 'HOG_{hog_id}.msa.fa', '{hog_id}.fa', etc.
    """
    search_dirs = list(custom_dirs or []) + LOCAL_PATHS.get(dataset_key, [])

    candidate_names = [
        f"{hog_id}.msa.fa",
        f"HOG_{hog_id}.msa.fa",
        f"{hog_id}.fa",
        f"HOG_{hog_id}.fa",
        f"{hog_id}.fasta",
        f"HOG_{hog_id}.fasta",
        f"{hog_id}.amp.fa",
        f"HOG_{hog_id}.amp.fa"
    ]

    for d in search_dirs:
        if not d.exists():
            continue
        for name in candidate_names:
            p = d / name
            if p.is_file() and p.stat().st_size > 0:
                return p
    return None


def process_single_hog(args_tuple: Tuple[str, str, Path, Dict[str, str]]) -> Optional[Dict[str, Any]]:
    """
    Worker task to process a single HOG FASTA file.
    Filters sequences to match valid species from taxa_dict, then computes SSS and Identity.
    """
    dataset_key, hog_id, fasta_path, taxa_dict = args_tuple

    try:
        records = list(SeqIO.parse(str(fasta_path), "fasta"))
        if not records:
            return None

        valid_seqs = []
        valid_ids = []
        seq_lens = []

        for rec in records:
            header_id = rec.id.split()[0]
            matched_code = None
            for code in taxa_dict.keys():
                if header_id.startswith(code):
                    matched_code = code
                    break

            if matched_code is None:
                for code in taxa_dict.keys():
                    if code in header_id:
                        matched_code = code
                        break

            if matched_code is not None or len(taxa_dict) == 0:
                clean_seq = str(rec.seq).replace("-", "").strip().upper()
                if len(clean_seq) > 0:
                    valid_seqs.append(clean_seq)
                    valid_ids.append(header_id)
                    seq_lens.append(len(clean_seq))

        if len(valid_seqs) < 2:
            return None

        calc_result = compute_pairwise_sss_and_identity(valid_seqs, valid_ids)

        return {
            "dataset": dataset_key,
            "hog_id": hog_id,
            "fasta_file": fasta_path.name,
            "num_taxa": len(valid_seqs),
            "num_pairs": calc_result["num_pairs"],
            "mean_seq_len": float(np.mean(seq_lens)),
            "min_seq_len": int(np.min(seq_lens)),
            "max_seq_len": int(np.max(seq_lens)),
            "mean_sss": calc_result["mean_sss"],
            "median_sss": calc_result["median_sss"],
            "min_sss": calc_result["min_sss"],
            "max_sss": calc_result["max_sss"],
            "std_sss": calc_result["std_sss"],
            "mean_identity": calc_result["mean_identity"],
            "median_identity": calc_result["median_identity"],
            "min_identity": calc_result["min_identity"],
            "max_identity": calc_result["max_identity"],
            "std_identity": calc_result["std_identity"],
            "pairs": calc_result["pairs"]
        }
    except Exception as e:
        print(f"Error processing {dataset_key} {hog_id} ({fasta_path}): {e}", file=sys.stderr)
        return None


def run_analysis(
    dataset_dir: Path = DEFAULT_DATASET_DIR,
    out_dir: Path = DEFAULT_OUT_DIR,
    custom_seq_dirs: Optional[Dict[str, List[Path]]] = None,
    num_workers: int = 4
) -> None:
    """Executes the SSS and identity benchmark across all 4 datasets."""
    out_dir.mkdir(parents=True, exist_ok=True)

    datasets = ["amn", "euk", "bac", "act"]
    dataset_names = {
        "amn": "Amniota",
        "euk": "Eukaryota",
        "bac": "Bacteria",
        "act": "Actinobacteria"
    }

    all_hog_summaries: List[Dict[str, Any]] = []
    all_pairs_details: List[Dict[str, Any]] = []

    print("================================================================================")
    print("  HOG Sequence Similarity Score (SSS) & Identity Benchmark")
    print(f"  Dataset Directory : {dataset_dir}")
    print(f"  Output Directory  : {out_dir}")
    print("================================================================================\n")

    for d_key in datasets:
        d_name = dataset_names[d_key]
        hogs_file = dataset_dir / f"{d_key}_hogs.tsv"
        taxa_file = dataset_dir / f"{d_key}_taxa.tsv"

        if not hogs_file.exists() or not taxa_file.exists():
            print(f"[{d_name}] Skipped: {hogs_file.name} or {taxa_file.name} not found.")
            continue

        hogs_list = load_hogs_list(hogs_file)
        taxa_dict = load_taxa_list(taxa_file)
        print(f"[{d_name}] Target HOGs: {len(hogs_list)}, Defined Taxa: {len(taxa_dict)}")

        custom_dirs = custom_seq_dirs.get(d_key, []) if custom_seq_dirs else []
        tasks = []
        missing_count = 0

        for h in hogs_list:
            fpath = find_hog_fasta_file(d_key, h, custom_dirs)
            if fpath is not None:
                tasks.append((d_key, h, fpath, taxa_dict))
            else:
                missing_count += 1

        print(f"  -> Found sequence files for {len(tasks)} / {len(hogs_list)} HOGs.")
        if missing_count > 0:
            print(f"  -> {missing_count} HOGs do not have sequence files locally yet.")

        if not tasks:
            print(f"  -> Skipping calculations for {d_name} (no local files).\n")
            continue

        print(f"  -> Computing SSS & Identity using {num_workers} parallel workers...")
        t0 = time.time()
        completed = 0

        with ProcessPoolExecutor(max_workers=num_workers) as executor:
            future_to_hog = {executor.submit(process_single_hog, t): t[1] for t in tasks}
            for future in as_completed(future_to_hog):
                res = future.result()
                if res is not None:
                    pairs = res.pop("pairs")
                    all_hog_summaries.append(res)
                    for p in pairs:
                        p["dataset"] = d_key
                        p["hog_id"] = res["hog_id"]
                        all_pairs_details.append(p)
                    completed += 1

        elapsed = time.time() - t0
        print(f"  -> Completed {completed} HOGs in {elapsed:.2f}s.\n")

    if not all_hog_summaries:
        print("No HOGs were successfully processed. Please check sequence paths.")
        return

    df_hogs = pd.DataFrame(all_hog_summaries)
    df_pairs = pd.DataFrame(all_pairs_details)

    hog_csv_path = out_dir / "hog_similarity_summary.csv"
    pairs_csv_path = out_dir / "pairwise_similarity_details.csv"

    df_hogs.to_csv(hog_csv_path, index=False)
    print(f"Saved HOG-level summary to: {hog_csv_path} ({len(df_hogs)} records)")

    df_pairs.to_csv(pairs_csv_path, index=False)
    print(f"Saved pairwise details to: {pairs_csv_path} ({len(df_pairs)} pairs)")

    dataset_summary = []
    for d_key, group in df_hogs.groupby("dataset"):
        dataset_summary.append({
            "dataset": d_key,
            "name": dataset_names.get(d_key, d_key),
            "hogs_evaluated": len(group),
            "mean_sss": group["mean_sss"].mean(),
            "median_sss": group["median_sss"].median(),
            "min_sss": group["min_sss"].min(),
            "max_sss": group["max_sss"].max(),
            "mean_identity_pct": group["mean_identity"].mean() * 100.0,
            "median_identity_pct": group["median_identity"].median() * 100.0,
            "min_identity_pct": group["min_identity"].min() * 100.0,
            "max_identity_pct": group["max_identity"].max() * 100.0,
            "hogs_sss_le_0_10": int((group["mean_sss"] <= 0.10).sum()),
            "pct_sss_le_0_10": float((group["mean_sss"] <= 0.10).mean() * 100.0)
        })

    df_ds = pd.DataFrame(dataset_summary)
    ds_csv_path = out_dir / "dataset_level_summary.csv"
    df_ds.to_csv(ds_csv_path, index=False)
    print(f"Saved dataset-level summary to: {ds_csv_path}\n")

    print("=== Dataset Summary Overview ===")
    print(df_ds.to_string(index=False))
    print("================================\n")


def main():
    parser = argparse.ArgumentParser(description="Analyze SSS and Identity for HOGs in dataset/")
    parser.add_argument("--dataset_dir", type=Path, default=DEFAULT_DATASET_DIR, help="Path to dataset directory")
    parser.add_argument("--out_dir", type=Path, default=DEFAULT_OUT_DIR, help="Output directory")
    parser.add_argument("--workers", type=int, default=6, help="Number of worker processes")
    args = parser.parse_args()

    run_analysis(dataset_dir=args.dataset_dir, out_dir=args.out_dir, num_workers=args.workers)


if __name__ == "__main__":
    main()
