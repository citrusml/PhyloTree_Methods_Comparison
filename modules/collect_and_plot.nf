nextflow.enable.dsl=2

process COLLECT_AND_PLOT {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path csv_files

    output:
    path("benchmark_summary.csv")
    path("regime_map_delta_nrf.png"), optional: true
    path("nrf_boxplots.png"), optional: true
    path("summary_statistics.csv"), optional: true
    path("method_comparisons"), optional: true
    path("length_analysis"), optional: true
    path("replications"), optional: true

    script:
    """
    head -n 1 \$(ls *.csv | head -1) > benchmark_summary.csv
    for f in *.csv; do tail -n +2 "\$f"; done >> benchmark_summary.csv
    python3 ${moduleDir}/../bin/plot/generate_all_reports.py --csv benchmark_summary.csv --outdir .
    """
}
