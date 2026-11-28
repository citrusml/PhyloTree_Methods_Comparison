#!/bin/bash
#BSUB -J gb_taxa_128_256
#BSUB -q mafft
#BSUB -n 4
#BSUB -M 16GB
#BSUB -W 72:00
#BSUB -o logs/graph_benchmark/taxa_128_256_bsub.out
#BSUB -e logs/graph_benchmark/taxa_128_256_bsub.err

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

run_large_taxa() {
    local n_taxa="$1"
    local chunk="$2"
    local outdir="results/graph_benchmark/taxa_N${n_taxa}_L500_v06_lba1"
    local workdir="work_gb_taxa_N${n_taxa}"

    mkdir -p "${outdir}"
    echo "=== Running Large Taxa Sweep: N=${n_taxa} (L=500, v=0.6, LBA=1.0) ==="
    nextflow run main.nf \
        -c next_configs/graph_benchmark/large_taxa_128_256.config \
        -work-dir "${workdir}" \
        -profile supercomputer \
        --taxa "${n_taxa}" \
        --chunk_size "${chunk}" \
        --outdir "${outdir}" \
        -resume \
        > "logs/graph_benchmark/taxa_N${n_taxa}.log" 2>&1
}

run_large_taxa 128 2
run_large_taxa 256 1

echo "=== Taxa Sweep (128, 256) completed ==="
