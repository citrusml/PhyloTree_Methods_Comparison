#!/bin/bash
#BSUB -J gb_nelsi_v
#BSUB -q mafft
#BSUB -n 4
#BSUB -M 16GB
#BSUB -W 72:00
#BSUB -o logs/graph_benchmark/nelsi_v_bsub.out
#BSUB -e logs/graph_benchmark/nelsi_v_bsub.err

set -euo pipefail

export PATH="$HOME/bin:$PATH"
eval "$(micromamba shell hook --shell bash)"
micromamba activate phylomethod_env
export PATH="$CONDA_PREFIX/bin:$HOME/bin:$PATH"
export NELSI_RSCRIPT="${MAMBA_ROOT_PREFIX:-$HOME/.micromamba}/envs/nelsi_env/bin/Rscript"

if [ -n "${LS_SUBCWD:-}" ] && [ -f "${LS_SUBCWD}/main.nf" ]; then
    REPO_ROOT="${LS_SUBCWD}"
elif [ -f "/lustre10/home/citrusml2004/PhyloTree_Methods_Comparison/main.nf" ]; then
    REPO_ROOT="/lustre10/home/citrusml2004/PhyloTree_Methods_Comparison"
elif [ -f "${HOME}/PhyloTree_Methods_Comparison/main.nf" ]; then
    REPO_ROOT="${HOME}/PhyloTree_Methods_Comparison"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
fi
cd "${REPO_ROOT}"

mkdir -p logs/graph_benchmark results/graph_benchmark
export NXF_OPTS="-Xms2g -Xmx12g"

run_nelsi_v() {
    local v_val="$1"
    local v_tag="$2"
    local outdir="$3"
    local workdir="work_gb_nelsi_${v_tag}"

    mkdir -p "${outdir}"
    echo "=== Running NELSI v Sweep: v=${v_val} (N=32, L=500, LBA=1.0) ==="
    nextflow run main.nf \
        -c next_configs/graph_benchmark/base_small_medium.config \
        -work-dir "${workdir}" \
        -profile supercomputer \
        --nelsi_v "${v_val}" \
        --outdir "${outdir}" \
        -resume \
        > "logs/graph_benchmark/nelsi_${v_tag}.log" 2>&1
}

run_nelsi_v 0.0 "v00" "results/graph_benchmark/nelsi_N32_L500_v00_lba1"
run_nelsi_v 0.3 "v03" "results/graph_benchmark/nelsi_N32_L500_v03_lba1"
# v=0.6 (デフォルト点) は run_taxa_16_32_64.sh の default_N32_L500_v06_lba1 で実行済みのため省略
run_nelsi_v 1.0 "v10" "results/graph_benchmark/nelsi_N32_L500_v10_lba1"

echo "=== NELSI v Sweep (0.0, 0.3, 1.0; v=0.6 via default) completed ==="
