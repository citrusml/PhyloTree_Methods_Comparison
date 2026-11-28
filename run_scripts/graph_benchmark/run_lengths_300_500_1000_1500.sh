#!/bin/bash
#BSUB -J gb_lengths
#BSUB -q mafft
#BSUB -n 4
#BSUB -M 16GB
#BSUB -W 72:00
#BSUB -o logs/graph_benchmark/lengths_bsub.out
#BSUB -e logs/graph_benchmark/lengths_bsub.err

set -euo pipefail

export PATH="$HOME/bin:$PATH"
eval "$(micromamba shell hook --shell bash)"
micromamba activate phylomethod_env
export PATH="$CONDA_PREFIX/bin:$HOME/bin:$PATH"
export NELSI_RSCRIPT="${MAMBA_ROOT_PREFIX:-$HOME/.micromamba}/envs/nelsi_env/bin/Rscript"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
cd "${REPO_ROOT}"

mkdir -p logs/graph_benchmark results/graph_benchmark
export NXF_OPTS="-Xms2g -Xmx12g"

run_length() {
    local seq_len="$1"
    local outdir="$2"
    local workdir="work_gb_length_L${seq_len}"

    mkdir -p "${outdir}"
    echo "=== Running Length Sweep: L=${seq_len} (N=32, v=0.6, LBA=1.0) ==="
    nextflow run main.nf \
        -c next_configs/graph_benchmark/base_small_medium.config \
        -work-dir "${workdir}" \
        -profile supercomputer \
        --lengths "${seq_len}" \
        --outdir "${outdir}" \
        -resume \
        > "logs/graph_benchmark/length_L${seq_len}.log" 2>&1
}

run_length 300  "results/graph_benchmark/length_N32_L300_v06_lba1"
# L=500 (デフォルト点) は run_taxa_16_32_64.sh の default_N32_L500_v06_lba1 で実行済みのため省略
run_length 1000 "results/graph_benchmark/length_N32_L1000_v06_lba1"
run_length 1500 "results/graph_benchmark/length_N32_L1500_v06_lba1"

echo "=== Length Sweep (300, 1000, 1500; L=500 via default) completed ==="
