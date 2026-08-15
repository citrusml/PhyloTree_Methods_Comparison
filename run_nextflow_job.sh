#!/bin/bash
#BSUB -J nextflow_main
#BSUB -q mafft
#BSUB -n 1
#BSUB -M 8GB
#BSUB -W 72:00
#BSUB -o /dev/null
#BSUB -e /dev/null
export PATH="$HOME/bin:$PATH"
eval "$(micromamba shell hook --shell bash)"
micromamba activate phylomethod_env

cd /lustre10/home/citrusml2004/PhyloTree_Methods_Comparison

# ロックファイルクリーンアップ
rm -f .nextflow/cache/*/db/LOCK

# リアルタイムに進捗ログを出力して実行
nextflow run main.nf -c nextflow.config -profile supercomputer -resume