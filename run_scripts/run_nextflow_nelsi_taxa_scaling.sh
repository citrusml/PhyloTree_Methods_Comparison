#!/bin/bash
#BSUB -J nextflow_nelsi_taxa_scaling
#BSUB -q mafft
#BSUB -n 1
#BSUB -M 16GB
#BSUB -W 120:00
#BSUB -o logs/nextflow_nelsi_taxa_scaling_bsub.out
#BSUB -e logs/nextflow_nelsi_taxa_scaling_bsub.err

export PATH="$HOME/bin:$PATH"
eval "$(micromamba shell hook --shell bash)"
micromamba activate phylomethod_env
export PATH="$CONDA_PREFIX/bin:$HOME/bin:$PATH"

# NELSI 用 Rscript のパスを指定（実際の envs パスに合わせて変更）
export NELSI_RSCRIPT="${MAMBA_ROOT_PREFIX:-$HOME/.micromamba}/envs/nelsi_env/bin/Rscript"

cd /lustre10/home/citrusml2004/PhyloTree_Methods_Comparison
mkdir -p logs results/results_nelsi_taxa32_L500 results/results_nelsi_taxa64_L500 results/results_nelsi_taxa128_L500

bash bin/install_indelible.sh
export NXF_OPTS="-Xms2g -Xmx12g"

echo "=== Running NELSI Taxa=32, L=500 ==="
nextflow run main.nf \
    -c next_configs/nextflow_nelsi_taxa32.config \
    -work-dir work_nelsi_taxa32 \
    -profile supercomputer \
    -resume \
    > logs/nextflow_nelsi_taxa32.log 2>&1

echo "=== Running NELSI Taxa=64, L=500 ==="
nextflow run main.nf \
    -c next_configs/nextflow_nelsi_taxa64.config \
    -work-dir work_nelsi_taxa64 \
    -profile supercomputer \
    -resume \
    > logs/nextflow_nelsi_taxa64.log 2>&1

echo "=== Running NELSI Taxa=128, L=500 ==="
nextflow run main.nf \
    -c next_configs/nextflow_nelsi_taxa128.config \
    -work-dir work_nelsi_taxa128 \
    -profile supercomputer \
    -resume \
    > logs/nextflow_nelsi_taxa128.log 2>&1

