#!/bin/bash
#BSUB -J nelsi_wag_all
#BSUB -q mafft
#BSUB -n 4
#BSUB -M 16GB
#BSUB -W 72:00
#BSUB -o logs/nextflow_nelsi_wag_all_bsub.out
#BSUB -e logs/nextflow_nelsi_wag_all_bsub.err

set -euo pipefail

export PATH="$HOME/bin:$PATH"
eval "$(micromamba shell hook --shell bash)"
micromamba activate phylomethod_env
export PATH="$CONDA_PREFIX/bin:$HOME/bin:$PATH"

export NELSI_RSCRIPT="${MAMBA_ROOT_PREFIX:-$HOME/.micromamba}/envs/nelsi_env/bin/Rscript"

cd /lustre10/home/citrusml2004/PhyloTree_Methods_Comparison
mkdir -p logs \
    results/results_nelsi_wag_v00 \
    results/results_nelsi_wag_v03 \
    results/results_nelsi_wag_v06 \
    results/results_nelsi_wag_lba3 \
    results/results_nelsi_wag_lba5

export NXF_OPTS="-Xms2g -Xmx12g"

run_cond() {
    local suffix="$1"
    local cond_label="$2"
    local lba_flag="$3"

    echo "=== Running WAG+G4 (No ICS): ${cond_label} ==="
    nextflow run main.nf \
        -c "next_configs/nextflow_nelsi_wag_${suffix}.config" \
        -work-dir "work_nelsi_wag_${suffix}" \
        -profile supercomputer \
        -resume \
        > "logs/nextflow_nelsi_wag_${suffix}.log" 2>&1

    if [ -n "${lba_flag}" ]; then
        python3 analysis/eval_additivity_on_replications.py \
            --results_dir "results/results_nelsi_wag_${suffix}" \
            --condition "${cond_label}" \
            --lba_pair T1 T20 \
            --workers 4 \
            > "logs/eval_additivity_wag_${suffix}.log" 2>&1
    else
        python3 analysis/eval_additivity_on_replications.py \
            --results_dir "results/results_nelsi_wag_${suffix}" \
            --condition "${cond_label}" \
            --workers 4 \
            > "logs/eval_additivity_wag_${suffix}.log" 2>&1
    fi
}

run_cond "v00"  "wag_v0.0_lba1" ""
run_cond "v03"  "wag_v0.3_lba1" ""
run_cond "v06"  "wag_v0.6_lba1" ""
run_cond "lba3" "wag_v0.3_lba3" "lba"
run_cond "lba5" "wag_v0.3_lba5" "lba"

echo "=== All WAG+G4 (No ICS) x NELSI x LBA evaluations completed ==="
