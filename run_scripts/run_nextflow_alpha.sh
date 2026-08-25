#!/bin/bash
#BSUB -J nextflow_alpha
#BSUB -q mafft
#BSUB -n 1
#BSUB -M 8GB
#BSUB -W 48:00
#BSUB -o /dev/null
#BSUB -e /dev/null
export PATH="$HOME/bin:$PATH"
eval "$(micromamba shell hook --shell bash)"
micromamba activate phylomethod_env

cd /lustre10/home/citrusml2004/PhyloTree_Methods_Comparison
mkdir -p logs results/results_alpha

# リアルタイムに進捗ログを出力して実行 (いつでも tail -f logs/nextflow_alpha.log で確認可能)
nextflow run next_main/main_alpha.nf -c next_configs/nextflow_alpha.config -profile supercomputer -resume > logs/nextflow_alpha.log 2>&1
