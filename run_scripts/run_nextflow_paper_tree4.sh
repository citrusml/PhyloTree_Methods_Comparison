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
export PATH="$CONDA_PREFIX/bin:$HOME/bin:$PATH"

cd /lustre10/home/citrusml2004/PhyloTree_Methods_Comparison
mkdir -p logs results/results_paper_tree4

# Nextflow の確認・自動復旧 (前回の perl-bioperl 導入で競合削除された場合の自動再インストール)
if ! command -v nextflow >/dev/null 2>&1; then
    echo "[$(date)] Restoring nextflow into phylomethod_env..."
    micromamba install -y -n phylomethod_env -c bioconda -c conda-forge nextflow
    hash -r
fi

# BioPerl は Nextflow (openjdk/openssl) との依存競合を防ぐため独立環境 bioperl_env に配備
BIOPERL_Root="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}/envs/bioperl_env"
if [ ! -x "${BIOPERL_Root}/bin/perl" ] && [ ! -x "$HOME/.micromamba/envs/bioperl_env/bin/perl" ]; then
    echo "[$(date)] Creating isolated bioperl_env for Bio::Tree::RandomFactory..."
    micromamba create -y -n bioperl_env -c bioconda -c conda-forge perl-bioperl || true
fi

# INDELible バイナリの配備・確認
bash bin/install_indelible.sh

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

