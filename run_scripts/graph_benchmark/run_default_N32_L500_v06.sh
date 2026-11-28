#!/bin/bash
#BSUB -J gb_default_N32
#BSUB -q mafft
#BSUB -n 4
#BSUB -M 16GB
#BSUB -W 72:00
#BSUB -o logs/graph_benchmark/default_N32_bsub.out
#BSUB -e logs/graph_benchmark/default_N32_bsub.err

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

mkdir -p logs/graph_benchmark results/graph_benchmark/default_N32_L500_v06_lba1
export NXF_OPTS="-Xms2g -Xmx12g"

echo "=== Running Graph Benchmark Default (N=32, L=500, v=0.6, LBA=1.0) ==="
nextflow run main.nf \
    -c next_configs/graph_benchmark/base_small_medium.config \
    -work-dir work_gb_default_N32 \
    -profile supercomputer \
    -resume \
    > logs/graph_benchmark/default_N32.log 2>&1
