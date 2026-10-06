#!/bin/bash
#BSUB -J nextflow_exp19
#BSUB -q mafft
#BSUB -n 1
#BSUB -M 16GB
#BSUB -W 120:00
#BSUB -o logs/nextflow_exp19_bsub.out
#BSUB -e logs/nextflow_exp19_bsub.err

# 環境設定
export PATH="$HOME/bin:$PATH"
eval "$(micromamba shell hook --shell bash)"
micromamba activate phylomethod_env

if [ -d "$CONDA_PREFIX/lib/jvm" ] && [ -z "$JAVA_HOME" ]; then
    export JAVA_HOME="$CONDA_PREFIX/lib/jvm"
fi

cd /lustre10/home/citrusml2004/PhyloTree_Methods_Comparison
mkdir -p logs results/results_exp19 results/results_exp19_N64 results/results_exp19_N128

# INDELible バイナリの配備・確認
bash bin/install_indelible.sh

# BioPerl の配備・確認
if ! perl -MBio::Tree::RandomFactory -e 1 2>/dev/null; then
    echo "[$(date)] Installing perl-bioperl into phylomethod_env..."
    micromamba install -y -c bioconda -c conda-forge perl-bioperl
fi

# JVMヒープ
export NXF_OPTS="-Xms2g -Xmx12g"

# 実験19 (Exp19): 論文完全再現 ＋ MMseqs2由来PSA+NJ ＋ MAFFT-LINSI vs --auto ＋ WAG+G4固定モデル ＋ 多Taxonスケーリング
# Stage 1: N = 20 (D = 0.1 .. 3.0, L = 1000, 50 reps)
echo "[$(date)] 実験19 (Stage 1: N=20 Main Benchmark) 開始"
nextflow run main.nf \
    -c next_configs/nextflow_paper_tree_exp19.config \
    -profile supercomputer \
    -resume \
    > logs/nextflow_exp19_N20.log 2>&1

echo "[$(date)] 実験19 Stage 1 (N=20) 完了 (exit code: $?)"

# Stage 2: N = 64 多Taxonスケーリング (D = 0.5, 1.0, 2.0, 30 reps)
if [ -f "next_configs/nextflow_paper_tree_exp19_N64.config" ]; then
    echo "[$(date)] 実験19 (Stage 2: N=64 Taxon Scaling) 開始"
    nextflow run main.nf \
        -c next_configs/nextflow_paper_tree_exp19_N64.config \
        -profile supercomputer \
        -resume \
        > logs/nextflow_exp19_N64.log 2>&1
    echo "[$(date)] 実験19 Stage 2 (N=64) 完了 (exit code: $?)"
fi

# Stage 3: N = 128 多Taxonスケーリング (D = 0.5, 1.0, 2.0, 20 reps)
if [ -f "next_configs/nextflow_paper_tree_exp19_N128.config" ]; then
    echo "[$(date)] 実験19 (Stage 3: N=128 Taxon Scaling) 開始"
    nextflow run main.nf \
        -c next_configs/nextflow_paper_tree_exp19_N128.config \
        -profile supercomputer \
        -resume \
        > logs/nextflow_exp19_N128.log 2>&1
    echo "[$(date)] 実験19 Stage 3 (N=128) 完了 (exit code: $?)"
fi

echo "[$(date)] 実験19 (全ステージ) 完了"
