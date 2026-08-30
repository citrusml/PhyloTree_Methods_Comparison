#!/bin/bash
#BSUB -J nextflow_fastme_options
#BSUB -q mafft
#BSUB -n 1
#BSUB -M 16GB
#BSUB -W 72:00
#BSUB -o logs/nextflow_fastme_options_bsub.out
#BSUB -e logs/nextflow_fastme_options_bsub.err
export PATH="$HOME/bin:$PATH"
eval "$(micromamba shell hook --shell bash)"
micromamba activate phylomethod_env

cd /lustre10/home/citrusml2004/PhyloTree_Methods_Comparison
mkdir -p logs results/results_fastme_options

# JVMヒープを明示的に指定
export NXF_OPTS="-Xms2g -Xmx12g"

# リアルタイムに進捗ログを出力して実行 (いつでも tail -f logs/nextflow_fastme_options.log で確認可能)
nextflow run next_main/main_fastme_options.nf -c next_configs/nextflow_fastme_options.config -profile supercomputer -resume > logs/nextflow_fastme_options.log 2>&1
