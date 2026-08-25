#!/bin/bash
#BSUB -J nextflow_alpha
#BSUB -q mafft
#BSUB -n 1
#BSUB -M 16GB                          # 16GB: 実験3は条件数が多くJVMヒープが増大するため増量
#BSUB -W 72:00
#BSUB -o logs/nextflow_alpha_bsub.out  # BSUBジョブ自体のstdoutを記録（OOMなど診断に必須）
#BSUB -e logs/nextflow_alpha_bsub.err  # BSUBジョブ自体のstderrを記録
export PATH="$HOME/bin:$PATH"
eval "$(micromamba shell hook --shell bash)"
micromamba activate phylomethod_env

cd /lustre10/home/citrusml2004/PhyloTree_Methods_Comparison
mkdir -p logs results/results_alpha

# JVMヒープを明示的に指定（デフォルトの自動検出に任せない）
export NXF_OPTS="-Xms2g -Xmx12g"

# リアルタイムに進捗ログを出力して実行 (いつでも tail -f logs/nextflow_alpha.log で確認可能)
nextflow run next_main/main_alpha.nf -c next_configs/nextflow_alpha.config -profile supercomputer -resume > logs/nextflow_alpha.log 2>&1
