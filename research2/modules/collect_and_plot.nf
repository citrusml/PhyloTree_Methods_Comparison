nextflow.enable.dsl=2

process COLLECT_AND_PLOT_RESEARCH2 {
    tag "summary"

    publishDir "${moduleDir}/../results", mode: 'copy'

    input:
    path(csv_files)

    output:
    path("all_properties_research2.csv"), emit: all_csv
    path("summary_table_research2.csv"), emit: summary_csv
    path("summary_report.md"), emit: report
    path("plots/*.png"), emit: plots

    script:
    """
    mkdir -p plots

    # 1. Combine all CSVs
    python3 -c '
import glob, os, sys
import pandas as pd

files = sorted(glob.glob("props_*.csv"))
if not files:
    print("No props CSV files found!", file=sys.stderr)
    sys.exit(1)

dfs = [pd.read_csv(f) for f in files]
all_df = pd.concat(dfs, ignore_index=True)
all_df.to_csv("all_properties_research2.csv", index=False)
print(f"Combined {len(files)} CSVs into all_properties_research2.csv ({len(all_df)} rows)")
'

    # 2. Run plotting and table generation
    python3 ${moduleDir}/../scripts/plot_research2_comparison.py \\
        --csv all_properties_research2.csv \\
        --outdir plots

    # Copy summaries to root for publishing
    cp plots/summary_table_research2.csv .
    cp plots/summary_report.md .
    """
}
