nextflow.enable.dsl=2

process COLLECT_AND_PLOT {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path csv_files

    output:
    path("benchmark_summary.csv")
    path("sss_summary.csv"), optional: true
    path("regime_map_delta_nrf.png"), optional: true
    path("nrf_boxplots.png"), optional: true
    path("nrf_vs_sss_curves.png"), optional: true
    path("sss_distribution_by_distance.png"), optional: true
    path("sss_breakdown_report.csv"), optional: true
    path("summary_statistics.csv"), optional: true
    path("method_comparisons"), optional: true
    path("length_analysis"), optional: true
    path("replications"), optional: true

    script:
    """
    python3 -c '
import os, glob, pandas as pd

method_files = [f for f in glob.glob("chunk_*.csv") if not f.startswith("chunk_sss_")]
if method_files:
    df_methods = pd.concat([pd.read_csv(f) for f in method_files], ignore_index=True)
else:
    df_methods = pd.DataFrame()

sss_files = glob.glob("chunk_sss_*.csv")
if sss_files:
    df_sss = pd.concat([pd.read_csv(f) for f in sss_files], ignore_index=True).drop_duplicates(subset=["distance", "length", "replicate"])
    df_sss.to_csv("sss_summary.csv", index=False)
    if not df_methods.empty:
        df_merged = pd.merge(df_methods, df_sss, on=["distance", "length", "replicate"], how="left")
    else:
        df_merged = df_sss
else:
    df_merged = df_methods

df_merged.to_csv("benchmark_summary.csv", index=False)
print(f"[COLLECT_AND_PLOT] Aggregated {len(df_methods)} method rows and {len(sss_files)} SSS chunks into benchmark_summary.csv")
'
    python3 ${moduleDir}/../bin/plot/generate_all_reports.py --csv benchmark_summary.csv --outdir .
    """
}
