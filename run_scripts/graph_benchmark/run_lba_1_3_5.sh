#!/bin/bash
#BSUB -J gb_lba
#BSUB -q mafft
#BSUB -n 4
#BSUB -M 16GB
#BSUB -W 72:00
#BSUB -o logs/graph_benchmark/lba_bsub.out
#BSUB -e logs/graph_benchmark/lba_bsub.err

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

run_lba() {
    local lba_val="$1"
    local lba_tag="$2"
    local outdir="$3"
    local workdir="work_gb_lba_${lba_tag}"

    mkdir -p "${outdir}"
    echo "=== Running LBA Sweep: lba_ratio=${lba_val} (N=32, L=500, v=0.6) ==="
    nextflow run main.nf \
        -c next_configs/graph_benchmark/base_small_medium.config \
        -work-dir "${workdir}" \
        -profile supercomputer \
        --lba_ratio "${lba_val}" \
        --outdir "${outdir}" \
        -resume \
        > "logs/graph_benchmark/lba_${lba_tag}.log" 2>&1
}

# LBA=1.0 (デフォルト点) は run_taxa_16_32_64.sh の default_N32_L500_v06_lba1 で実行済みのため省略
run_lba 3.0 "lba3" "results/graph_benchmark/lba_N32_L500_v06_lba3"
run_lba 5.0 "lba5" "results/graph_benchmark/lba_N32_L500_v06_lba5"

echo "=== LBA Sweep (3x, 5x; 1x via default) completed ==="
