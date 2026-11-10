nextflow.enable.dsl=2

process COLLECT_AND_PLOT {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path csv_files

    output:
    path("benchmark_summary.csv")
    path("alignment_summary.csv"), optional: true
    path("sss_summary.csv"), optional: true
    path("regime_map_delta_nrf.png"), optional: true
    path("nrf_boxplots.png"), optional: true
    path("nrf_vs_sss_curves.png"), optional: true
    path("sp_score_vs_distance.png"), optional: true
    path("sp_score_vs_mae.png"), optional: true
    path("sp_score_vs_tree_accuracy.png"), optional: true
    path("error_propagation_analysis.png"), optional: true
    path("sss_distribution_by_distance.png"), optional: true
    path("sss_breakdown_report.csv"), optional: true
    path("summary_statistics.csv"), optional: true
    path("method_comparisons"), optional: true
    path("length_analysis"), optional: true
    path("replications"), optional: true

    script:
    """
    export MPLBACKEND=Agg

    # Resolve Python binary (prioritize phylomethod_env)
    PY_BIN=""
    for cand in \
        "${CONDA_PREFIX}/bin/python3" \
        "${CONDA_PREFIX}/bin/python" \
        "$HOME/.micromamba/envs/phylomethod_env/bin/python3" \
        "$HOME/micromamba/envs/phylomethod_env/bin/python3" \
        "${MAMBA_ROOT_PREFIX:-$HOME/micromamba}/envs/phylomethod_env/bin/python3" \
        "\$(which python3 2>/dev/null)" \
        "/usr/bin/python3"; do
        if [ -n "\$cand" ] && [ -x "\$cand" ]; then
            if "\$cand" -c "import seaborn" >/dev/null 2>&1; then
                PY_BIN="\$cand"
                break
            elif [ -z "\$PY_BIN" ]; then
                PY_BIN="\$cand"
            fi
        fi
    done

    if [ -z "\$PY_BIN" ]; then
        PY_BIN="python3"
    fi

    # Ensure seaborn is available
    if ! "\$PY_BIN" -c "import seaborn" >/dev/null 2>&1; then
        echo "[COLLECT_AND_PLOT] Warning: seaborn missing in \$PY_BIN, installing via pip..."
        "\$PY_BIN" -m pip install --quiet seaborn || true
    fi

    echo "[COLLECT_AND_PLOT] Using Python binary: \$PY_BIN"

    "\$PY_BIN" -c '
import os, glob, sys, traceback
import pandas as pd

def safe_read(f):
    try:
        if os.path.exists(f) and os.path.getsize(f) > 0:
            return pd.read_csv(f)
    except Exception as e:
        print(f"[Warning] Skipping corrupt file {f}: {e}", file=sys.stderr)
    return None

try:
    method_files = [f for f in glob.glob("chunk_*.csv") if not f.startswith("chunk_sss_") and not f.startswith("chunk_alignment_")]
    dfs = [safe_read(f) for f in method_files]
    dfs = [d for d in dfs if d is not None and not d.empty]
    if dfs:
        df_methods = pd.concat(dfs, ignore_index=True)
    else:
        df_methods = pd.DataFrame()

    aln_files = glob.glob("chunk_alignment_*.csv")
    aln_dfs = [safe_read(f) for f in aln_files]
    aln_dfs = [d for d in aln_dfs if d is not None and not d.empty]
    if aln_dfs:
        df_aln = pd.concat(aln_dfs, ignore_index=True).drop_duplicates(subset=["distance", "length", "replicate"])
        df_aln.to_csv("alignment_summary.csv", index=False)
        if not df_methods.empty:
            df_methods["distance"] = df_methods["distance"].round(4)
            df_aln["distance"] = df_aln["distance"].round(4)
            df_methods = pd.merge(df_methods, df_aln, on=["distance", "length", "replicate"], how="left")
        else:
            df_methods = df_aln

    sss_files = glob.glob("chunk_sss_*.csv")
    sss_dfs = [safe_read(f) for f in sss_files]
    sss_dfs = [d for d in sss_dfs if d is not None and not d.empty]
    if sss_dfs:
        df_sss = pd.concat(sss_dfs, ignore_index=True).drop_duplicates(subset=["distance", "length", "replicate"])
        df_sss.to_csv("sss_summary.csv", index=False)
        if not df_methods.empty:
            df_methods["distance"] = df_methods["distance"].round(4)
            df_sss["distance"] = df_sss["distance"].round(4)
            drop_cols = [c for c in ["num_taxa"] if c in df_sss.columns and c in df_methods.columns]
            if drop_cols:
                df_sss = df_sss.drop(columns=drop_cols)
            df_merged = pd.merge(df_methods, df_sss, on=["distance", "length", "replicate"], how="left")
        else:
            df_merged = df_sss
    else:
        df_merged = df_methods

    df_merged.to_csv("benchmark_summary.csv", index=False)
    print(f"[COLLECT_AND_PLOT] Aggregated {len(df_methods)} method rows, {len(aln_files)} alignment chunks, {len(sss_files)} SSS chunks into benchmark_summary.csv")
except Exception as e:
    traceback.print_exc()
    sys.exit(1)
'
    "\$PY_BIN" ${moduleDir}/../bin/plot/generate_all_reports.py --csv benchmark_summary.csv --outdir .
    """
}
