#!/bin/bash
#BSUB -J eval_nelsi_lba
#BSUB -q mafft
#BSUB -n 8
#BSUB -M 16GB
#BSUB -W 04:00
#BSUB -o logs/eval_nelsi_lba_bsub.out
#BSUB -e logs/eval_nelsi_lba_bsub.err

export PATH="$HOME/bin:$PATH"
eval "$(micromamba shell hook --shell bash)"
micromamba activate phylomethod_env
export PATH="$CONDA_PREFIX/bin:$HOME/bin:$PATH"

cd /lustre10/home/citrusml2004/PhyloTree_Methods_Comparison
mkdir -p logs results/nelsi_lba_specific_plots

DIRS=""
for d in results/results_nelsi_yule results/results_nelsi_lba3 results/results_nelsi_lbs3 results/results_nelsi_lba5 results/results_nelsi_lbs5; do
    if [ -d "$d" ]; then
        DIRS="$DIRS $d"
    fi
done

python3 analysis/evaluate_lba_metrics.py \
    --results_dirs $DIRS \
    --lba_pair T1 T20 \
    --workers 8 \
    --outcsv results/nelsi_lba_specific_metrics.csv \
    --summary_csv results/nelsi_lba_specific_summary.csv \
    --plot_dir results/nelsi_lba_specific_plots \
    > logs/eval_nelsi_lba.log 2>&1
