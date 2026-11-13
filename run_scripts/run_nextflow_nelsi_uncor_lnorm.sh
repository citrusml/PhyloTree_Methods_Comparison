#!/bin/bash
#BSUB -J nextflow_nelsi_uncor_lnorm
#BSUB -q mafft
#BSUB -n 1
#BSUB -M 16GB
#BSUB -W 120:00
#BSUB -o logs/nextflow_nelsi_uncor_lnorm_bsub.out
#BSUB -e logs/nextflow_nelsi_uncor_lnorm_bsub.err

export PATH="$HOME/bin:$PATH"
eval "$(micromamba shell hook --shell bash)"
micromamba activate phylomethod_env
export PATH="$CONDA_PREFIX/bin:$HOME/bin:$PATH"

# NELSI 用 Rscript のパスを指定（実際の envs パスに合わせて変更）
export NELSI_RSCRIPT="${MAMBA_ROOT_PREFIX:-$HOME/.micromamba}/envs/nelsi_env/bin/Rscript"

cd /lustre10/home/citrusml2004/PhyloTree_Methods_Comparison
mkdir -p logs results/results_nelsi_uncor_lnorm

bash bin/install_indelible.sh
export NXF_OPTS="-Xms2g -Xmx12g"

nextflow run main.nf \
    -c next_configs/nextflow_nelsi_uncor_lnorm.config \
    -work-dir work_nelsi_uncor_lnorm \
    -profile supercomputer \
    -resume \
    > logs/nextflow_nelsi_uncor_lnorm.log 2>&1

