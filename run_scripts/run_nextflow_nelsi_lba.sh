#!/bin/bash
#BSUB -J nextflow_nelsi_lba
#BSUB -q mafft
#BSUB -n 1
#BSUB -M 16GB
#BSUB -W 120:00
#BSUB -o logs/nextflow_nelsi_lba_bsub.out
#BSUB -e logs/nextflow_nelsi_lba_bsub.err

export PATH="$HOME/bin:$PATH"
eval "$(micromamba shell hook --shell bash)"
micromamba activate phylomethod_env
export PATH="$CONDA_PREFIX/bin:$HOME/bin:$PATH"

# NELSI 用 Rscript のパスを指定（実際の envs パスに合わせて変更）
export NELSI_RSCRIPT="${MAMBA_ROOT_PREFIX:-$HOME/.micromamba}/envs/nelsi_env/bin/Rscript"

cd /lustre10/home/citrusml2004/PhyloTree_Methods_Comparison
mkdir -p logs results/results_nelsi_lba3 results/results_nelsi_lba5

bash bin/install_indelible.sh
export NXF_OPTS="-Xms2g -Xmx12g"

echo "=== Running NELSI Yule LBA (taxa=20, lba_ratio=3.0) ==="
nextflow run main.nf \
    -c next_configs/nextflow_nelsi_lba3.config \
    -work-dir work_nelsi_lba3 \
    -profile supercomputer \
    -resume \
    > logs/nextflow_nelsi_lba3.log 2>&1

echo "=== Running NELSI Yule LBA (taxa=20, lba_ratio=5.0) ==="
nextflow run main.nf \
    -c next_configs/nextflow_nelsi_lba5.config \
    -work-dir work_nelsi_lba5 \
    -profile supercomputer \
    -resume \
    > logs/nextflow_nelsi_lba5.log 2>&1
