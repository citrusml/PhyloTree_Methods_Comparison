#!/bin/bash
#BSUB -J nextflow_nelsi_ics_full
#BSUB -q mafft
#BSUB -n 1
#BSUB -M 16GB
#BSUB -W 72:00
#BSUB -o logs/nextflow_nelsi_ics_full_bsub.out
#BSUB -e logs/nextflow_nelsi_ics_full_bsub.err

export PATH="$HOME/bin:$PATH"
eval "$(micromamba shell hook --shell bash)"
micromamba activate phylomethod_env
export PATH="$CONDA_PREFIX/bin:$HOME/bin:$PATH"

# NELSI 用 Rscript のパスを指定（実際の envs パスに合わせて変更）
export NELSI_RSCRIPT="${MAMBA_ROOT_PREFIX:-$HOME/.micromamba}/envs/nelsi_env/bin/Rscript"

cd /lustre10/home/citrusml2004/PhyloTree_Methods_Comparison
mkdir -p logs results/results_nelsi_ics_full_v03 results/results_nelsi_ics_full_v06

export NXF_OPTS="-Xms2g -Xmx12g"

echo "=== Running NELSI (v=0.3) x ICS Full (taxa=20, L=500, rep=20) ==="
nextflow run main.nf \
    -c next_configs/nextflow_nelsi_ics_full_v03.config \
    -work-dir work_nelsi_ics_full_v03 \
    -profile supercomputer \
    -resume \
    > logs/nextflow_nelsi_ics_full_v03.log 2>&1

echo "=== Running NELSI (v=0.6) x ICS Full (taxa=20, L=500, rep=20) ==="
nextflow run main.nf \
    -c next_configs/nextflow_nelsi_ics_full_v06.config \
    -work-dir work_nelsi_ics_full_v06 \
    -profile supercomputer \
    -resume \
    > logs/nextflow_nelsi_ics_full_v06.log 2>&1
