#!/bin/bash
#BSUB -J nextflow_taxon
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
mkdir -p logs results/results_taxon

# リアルタイムに進捗ログを出力して実行 (いつでも tail -f logs/nextflow_taxon.log で確認可能)
nextflow run next_main/main_taxon.nf -c next_configs/nextflow_taxon.config -profile supercomputer -resume > logs/nextflow_taxon.log 2>&1
