#!/bin/bash
#BSUB -J nextflow_paper_tree4
#BSUB -q mafft
#BSUB -n 1
#BSUB -M 16GB
#BSUB -W 120:00
#BSUB -o logs/nextflow_paper_tree4_bsub.out
#BSUB -e logs/nextflow_paper_tree4_bsub.err

# 環境設定
export PATH="$HOME/bin:$PATH"
eval "$(micromamba shell hook --shell bash)"
micromamba activate phylomethod_env

cd /lustre10/home/citrusml2004/PhyloTree_Methods_Comparison
mkdir -p logs results/results_paper_tree4

# INDELible バイナリの配備・確認
bash bin/install_indelible.sh

# BioPerl の配備・確認 (未インストールの場合は自動インストール)
if ! perl -MBio::Tree::RandomFactory -e 1 2>/dev/null; then
    echo "[$(date)] Installing perl-bioperl into phylomethod_env..."
    micromamba install -y -c bioconda -c conda-forge perl-bioperl
fi

# JVMヒープを明示的に指定
export NXF_OPTS="-Xms2g -Xmx12g"

# 実験18: 論文距離行列・PSA+NJ 完全一致ベンチマーク (Paper Tree 4)
# 論文 (Matsui & Iwasaki 2020, Systematic Biology, Fig. 3a) の距離計算・実験条件を厳密再現:
#   - 系統樹生成モデル: 後退 Yule 過程 (Backward Yule) + 対数分布枝長 (1 - ln(u*(e-1)+1))
#   - 配列シミュレータ: INDELible v1.03 (連続時間マルコフ過程, WAG + Gamma4, alpha=1.0)
#   - Indel 長さ分布: Zipfian (POW 1.7 50, べき指数 a=1.7, 最大 50 残基), 挿入率/欠失率 0.10/0.10
#   - PSA + NJ (PWA+NJ / MMSEQS+NJ):
#       MMseqs2 (-s 7.5 -e 10) による対称化相対配列類似度スコア w_ij を算出し、
#       d_ij = 1 - w_ij (PWA+NJ) および d_ij = -ln(w_ij) (MMSEQS+NJ) で距離行列を構築して NJ 推定
#   - MSA + NJ (MSA+NJ / MSA_LINSI+NJ / TRUE_MSA+NJ):
#       R phangorn::dist.ml に準拠した WAG + Gamma4 (alpha=1.0) 最尤距離により距離行列を構築して NJ 推定
#   - アライメント: MAFFT --auto および MAFFT L-INS-i (mafft-linsi) の双方を比較
echo "[$(date)] 実験18 (Paper Tree Benchmark 4 - 論文距離行列・PSA+NJ 完全一致版) 開始"
nextflow run main.nf \
    -c next_configs/nextflow_paper_tree4.config \
    -profile supercomputer \
    -resume \
    > logs/nextflow_paper_tree4.log 2>&1

echo "[$(date)] 実験18 完了 (exit code: $?)"

