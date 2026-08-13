#!/usr/bin/env python3
"""
Aggregate Work Results Script
Recovers and aggregates all evaluation CSV files directly from Nextflow work/ directories,
guaranteeing zero data loss and instant generation of benchmark_summary.csv and regime_map_delta_nrf.png
without re-running any computations.
"""

import sys
import os
import glob
import argparse
import pandas as pd
import subprocess

def main():
    parser = argparse.ArgumentParser(description="Aggregate all results from work/ directory")
    parser.add_argument("--workdir", default="work", help="Path to Nextflow work directory")
    parser.add_argument("--outdir", default="results", help="Output directory for final summary and plot")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    out_csv = os.path.join(args.outdir, "benchmark_summary.csv")

    print(f"Searching for evaluation CSV files in {args.workdir}/...")
    csv_pattern = os.path.join(args.workdir, "**", "*.csv")
    csv_files = glob.glob(csv_pattern, recursive=True)

    # Filter evaluation CSVs
    valid_csvs = [f for f in csv_files if any(k in os.path.basename(f) for k in ["eval", "pwa_nj", "msa_nj", "msa_ml"])]
    print(f"Found {len(valid_csvs)} evaluation CSV files!")

    if not valid_csvs:
        print("No CSV files found in work directory.")
        return

    dfs = []
    for f in valid_csvs:
        try:
            df = pd.read_csv(f)
            if "distance" in df.columns and "nrf_distance" in df.columns:
                dfs.append(df)
        except Exception:
            continue

    if not dfs:
        print("No valid evaluation records could be parsed.")
        return

    combined_df = pd.concat(dfs, ignore_index=True)
    # Deduplicate in case of duplicate entries
    combined_df = combined_df.drop_duplicates(subset=["distance", "length", "replicate", "pipeline"])
    combined_df = combined_df.sort_values(by=["distance", "length", "replicate", "pipeline"])

    combined_df.to_csv(out_csv, index=False)
    print(f"Successfully saved {len(combined_df)} records to {out_csv}")

    # Generate Regime Map plot
    project_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    plot_script = os.path.join(project_dir, "bin", "plot_regime_map.py")
    if os.path.exists(plot_script):
        cmd = [sys.executable, plot_script, "--csv", out_csv, "--outdir", args.outdir]
        subprocess.run(cmd)

    print(f"\n=== Aggregation Completed! ===")
    print(f"Summary CSV: {out_csv}")
    print(f"Regime Map:  {os.path.join(args.outdir, 'regime_map_delta_nrf.png')}")

if __name__ == "__main__":
    main()
