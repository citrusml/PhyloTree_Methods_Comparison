#!/bin/bash
#BSUB -J nextflow_ics_full
#BSUB -q mafft
#BSUB -n 1
#BSUB -M 16GB
#BSUB -W 72:00
#BSUB -o logs/nextflow_ics_full_bsub.out
#BSUB -e logs/nextflow_ics_full_bsub.err
export PATH="$HOME/bin:$PATH"
eval "$(micromamba shell hook --shell bash)"
micromamba activate phylomethod_env

cd /lustre10/home/citrusml2004/PhyloTree_Methods_Comparison
mkdir -p logs results/results_ics_full

# JVMヒープを明示的に指定
export NXF_OPTS="-Xms2g -Xmx12g"

# リアルタイムに進捗ログを出力して実行
nextflow run next_main/main_ics_full.nf -c next_configs/nextflow_ics_full.config -profile supercomputer -resume > logs/nextflow_ics_full.log 2>&1
