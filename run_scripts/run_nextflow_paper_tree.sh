#!/bin/bash
#BSUB -J nextflow_paper_tree
#BSUB -q mafft
#BSUB -n 1
#BSUB -M 16GB
#BSUB -W 120:00
#BSUB -o logs/nextflow_paper_tree_bsub.out
#BSUB -e logs/nextflow_paper_tree_bsub.err

# 環境設定
export PATH="$HOME/bin:$PATH"
eval "$(micromamba shell hook --shell bash)"
micromamba activate phylomethod_env

cd /lustre10/home/citrusml2004/PhyloTree_Methods_Comparison
mkdir -p logs results/results_paper_tree

# JVMヒープを明示的に指定
export NXF_OPTS="-Xms2g -Xmx12g"

# 実験15: 論文準拠 系統樹生成モデル (Yule過程 + 対数抽出枝長) ベンチマーク
# 論文 (Matsui & Iwasaki 2020) の条件を再現:
#   - 系統樹生成モデル: 後退 Yule 過程 (Backward Yule) + 対数分布枝長 (1 - ln(u*(e-1)+1))
#   - Indel 長さ分布: Zipfian (POW{1.7/50}, べき指数 a=1.7, 最大 50 残基)
#   - 挿入率 / 欠失率: 0.10 / 0.10 (置換率に対する相対値: 10% 引き上げ、論文 Fig. 3a 条件)
#   - 比較手法: PWA+NJ, MSA+NJ, MSA+ML (IQ-TREE 2), MSA+RAXML (RAxML -f d), GS
#   - ギャップペナルティ: Gap Open=20, Extension=2 (PhyPA準拠)
echo "[$(date)] 実験15 (Paper Tree Benchmark) 開始"
nextflow run main.nf \
    -c next_configs/nextflow_paper_tree.config \
    -profile supercomputer \
    -resume \
    > logs/nextflow_paper_tree.log 2>&1

echo "[$(date)] 実験15 完了 (exit code: $?)"
