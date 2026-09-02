#!/bin/bash
#BSUB -J nextflow_zipfian_indel
#BSUB -q mafft
#BSUB -n 1
#BSUB -M 16GB
#BSUB -W 120:00
#BSUB -o logs/nextflow_zipfian_indel_bsub.out
#BSUB -e logs/nextflow_zipfian_indel_bsub.err

# 環境設定
export PATH="$HOME/bin:$PATH"
eval "$(micromamba shell hook --shell bash)"
micromamba activate phylomethod_env

cd /lustre10/home/citrusml2004/PhyloTree_Methods_Comparison
mkdir -p logs results/results_zipfian_indel

# JVMヒープを明示的に指定
export NXF_OPTS="-Xms2g -Xmx12g"

# 実験14: Zipfian べき乗則 Indel + 高ギャップペナルティ
# 論文 (Matsui & Iwasaki 2020) の条件を再現:
#   - Indel 長さ分布: Zipfian (POW{1.7/50}, べき指数 a=1.7, 最大 50 残基)
#   - 挿入率 / 欠失率: 0.05 / 0.05 (置換率に対する相対値)
#   - ギャップペナルティ: Gap Open=20, Extension=2 (PhyPA準拠)
echo "[$(date)] 実験14 (Zipfian Indel Benchmark) 開始"
nextflow run next_main/main.nf \
    -c next_configs/nextflow_zipfian_indel.config \
    -profile supercomputer \
    -resume \
    > logs/nextflow_zipfian_indel.log 2>&1

echo "[$(date)] 実験14 完了 (exit code: $?)"
