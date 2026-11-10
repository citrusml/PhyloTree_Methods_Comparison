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
mkdir -p "$HOME/bin"
export PATH="$CONDA_PREFIX/bin:$HOME/bin:$PATH"

cd /lustre10/home/citrusml2004/PhyloTree_Methods_Comparison
mkdir -p logs results/results_paper_tree4

# 1. phylomethod_env 内に混入した perl-bioperl の除去 (openjdk/mmseqs2 との依存競合の原因を排除)
if micromamba list -n phylomethod_env 2>/dev/null | grep -q "perl-bioperl"; then
    echo "[$(date)] Removing conflicting perl-bioperl from phylomethod_env..."
    micromamba remove -y -n phylomethod_env perl-bioperl || true
fi

# 2. phylomethod_env 必須ツールの一括確認・自動復旧 (nextflow, mmseqs2, mafft, rapidnj, fastme, iqtree, raxml)
MISSING_PKGS=""
command -v nextflow >/dev/null 2>&1 || MISSING_PKGS="$MISSING_PKGS nextflow"
command -v mmseqs   >/dev/null 2>&1 || MISSING_PKGS="$MISSING_PKGS mmseqs2"
command -v mafft    >/dev/null 2>&1 || MISSING_PKGS="$MISSING_PKGS mafft"
command -v rapidnj  >/dev/null 2>&1 || MISSING_PKGS="$MISSING_PKGS rapidnj"
command -v fastme   >/dev/null 2>&1 || MISSING_PKGS="$MISSING_PKGS fastme"
(command -v iqtree2 >/dev/null 2>&1 || command -v iqtree >/dev/null 2>&1) || MISSING_PKGS="$MISSING_PKGS iqtree"
(command -v raxmlHPC-PTHREADS >/dev/null 2>&1 || command -v raxmlHPC >/dev/null 2>&1) || MISSING_PKGS="$MISSING_PKGS raxml"

if [ -n "$MISSING_PKGS" ]; then
    echo "[$(date)] Restoring missing packages into phylomethod_env:$MISSING_PKGS"
    micromamba install -y -n phylomethod_env -c conda-forge -c bioconda \
        nextflow mmseqs2 mafft rapidnj fastme iqtree raxml || true
    hash -r
fi

# 3. 万一 phylomethod_env 内で mmseqs2 が競合した場合の独立環境 mmseqs_env 配備 + シンボリックリンク作成
if ! command -v mmseqs >/dev/null 2>&1; then
    echo "[$(date)] Creating isolated mmseqs_env for MMseqs2..."
    micromamba create -y -n mmseqs_env -c conda-forge -c bioconda mmseqs2 || true
    for cand in "${MAMBA_ROOT_PREFIX:-$HOME/micromamba}/envs/mmseqs_env/bin/mmseqs" "$HOME/.micromamba/envs/mmseqs_env/bin/mmseqs"; do
        if [ -x "$cand" ]; then
            ln -sf "$cand" "$HOME/bin/mmseqs"
            [ -n "$CONDA_PREFIX" ] && ln -sf "$cand" "$CONDA_PREFIX/bin/mmseqs"
            break
        fi
    done
    hash -r
fi

# 4. BioPerl は独立環境 bioperl_env に配備
BIOPERL_Root="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}/envs/bioperl_env"
if [ ! -x "${BIOPERL_Root}/bin/perl" ] && [ ! -x "$HOME/.micromamba/envs/bioperl_env/bin/perl" ]; then
    echo "[$(date)] Creating isolated bioperl_env for Bio::Tree::RandomFactory..."
    micromamba create -y -n bioperl_env -c bioconda -c conda-forge perl-bioperl || true
fi

# 5. INDELible バイナリの配備・確認
bash bin/install_indelible.sh

echo "[$(date)] Tool verification:"
echo "  nextflow: $(command -v nextflow)"
echo "  mmseqs:   $(command -v mmseqs)"
echo "  mafft:    $(command -v mafft)"
echo "  rapidnj:  $(command -v rapidnj)"

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

