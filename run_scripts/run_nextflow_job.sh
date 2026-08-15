#!/bin/bash
#BSUB -J nextflow_main
#BSUB -q mafft
#BSUB -n 1
#BSUB -M 8GB
#BSUB -W 24:00
#BSUB -o /dev/null
#BSUB -e /dev/null
export PATH="$HOME/bin:$PATH"
eval "$(micromamba shell hook --shell bash)"
micromamba activate phylomethod_env

cd /lustre10/home/citrusml2004/PhyloTree_Methods_Comparison
mkdir -p logs

# リアルタイムに進捗ログを出力して実行 (いつでも tail -f logs/nextflow_live.log で確認可能)
nextflow run /next_main/main.nf -c next_configs/nextflow.config -profile supercomputer -resume > logs/nextflow.log 2>&1