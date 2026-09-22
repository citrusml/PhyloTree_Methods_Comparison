#!/bin/bash
#BSUB -J nextflow_paper_tree3
#BSUB -q mafft
#BSUB -n 1
#BSUB -M 16GB
#BSUB -W 120:00
#BSUB -o logs/nextflow_paper_tree3_bsub.out
#BSUB -e logs/nextflow_paper_tree3_bsub.err

# 環境設定
export PATH="$HOME/bin:$PATH"
eval "$(micromamba shell hook --shell bash)"
micromamba activate phylomethod_env

cd /lustre10/home/citrusml2004/PhyloTree_Methods_Comparison
mkdir -p logs results/results_paper_tree3

# INDELible バイナリの配備・確認
bash bin/install_indelible.sh

# JVMヒープを明示的に指定
export NXF_OPTS="-Xms2g -Xmx12g"

# 実験17: 論文完全再現 + 誤差伝播・中間生データ保存ベンチマーク (Paper Tree 3)
# 論文 (Matsui & Iwasaki 2020) の条件を厳密再現:
#   - 系統樹生成モデル: 後退 Yule 過程 (Backward Yule) + 対数分布枝長 (1 - ln(u*(e-1)+1))
#   - 配列シミュレータ: INDELible v1.03 (連続時間マルコフ過程)
#   - Indel 長さ分布: Zipfian (POW 1.7 50, べき指数 a=1.7, 最大 50 残基)
#   - 挿入率 / 欠失率: 0.10 / 0.10 (論文 Fig. 3a 最強条件)
#   - アミノ酸置換モデル: WAG + Gamma4 (alpha=1.0)
#   - 比較手法: PWA+NJ, MSA+NJ, MSA+ML (IQ-TREE 2), MSA+RAXML (RAxML -f d), GS
#   - 新機能: SP score (f_D) 計算, レプリケーション生データ永続保存, パトリスティック距離誤差伝播解析
echo "[$(date)] 実験17 (Paper Tree Benchmark 3 - 原著完全再現) 開始"
nextflow run main.nf \
    -c next_configs/nextflow_paper_tree3.config \
    -profile supercomputer \
    -resume \
    > logs/nextflow_paper_tree3.log 2>&1

echo "[$(date)] 実験17 完了 (exit code: $?)"
