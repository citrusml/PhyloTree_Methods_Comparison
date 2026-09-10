#!/bin/bash
#BSUB -J research2_benchmark
#BSUB -q mafft
#BSUB -n 1
#BSUB -M 16GB
#BSUB -W 120:00
#BSUB -o logs/research2_bsub.out
#BSUB -e logs/research2_bsub.err

# 環境設定
export PATH="$HOME/bin:$PATH"
eval "$(micromamba shell hook --shell bash)"
micromamba activate phylomethod_env

cd /lustre10/home/citrusml2004/PhyloTree_Methods_Comparison
mkdir -p logs research2/results research2/plots

# INDELible バイナリの配備・確認
bash bin/install_indelible.sh

# JVMヒープを明示的に指定
export NXF_OPTS="-Xms2g -Xmx12g"

echo "[$(date)] research2 (AliSim vs INDELible 比較実験) 開始"
nextflow run research2/main_research2.nf \
    -c research2/configs/research2.config \
    -profile supercomputer \
    -resume \
    > logs/research2_nextflow.log 2>&1

echo "[$(date)] research2 完了 (exit code: $?)"
