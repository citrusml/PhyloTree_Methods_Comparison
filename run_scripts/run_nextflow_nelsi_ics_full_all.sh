#!/bin/bash
#BSUB -J nelsi_ics_full_all
#BSUB -q mafft
#BSUB -n 4
#BSUB -M 16GB
#BSUB -W 72:00
#BSUB -o logs/nextflow_nelsi_ics_full_all_bsub.out
#BSUB -e logs/nextflow_nelsi_ics_full_all_bsub.err

set -euo pipefail

export PATH="$HOME/bin:$PATH"
eval "$(micromamba shell hook --shell bash)"
micromamba activate phylomethod_env
export PATH="$CONDA_PREFIX/bin:$HOME/bin:$PATH"

export NELSI_RSCRIPT="${MAMBA_ROOT_PREFIX:-$HOME/.micromamba}/envs/nelsi_env/bin/Rscript"

cd /lustre10/home/citrusml2004/PhyloTree_Methods_Comparison
mkdir -p logs \
    results/results_nelsi_ics_full_v00 \
    results/results_nelsi_ics_full_lba3 \
    results/results_nelsi_ics_full_lba5

export NXF_OPTS="-Xms2g -Xmx12g"

echo "=== Step 1: Reuse existing v=0.3 and v=0.6 (evaluate 7 additivity/FastME methods on existing PHYLIP matrices) ==="
python3 analysis/eval_additivity_on_replications.py \
    --results_dir results/results_nelsi_ics_full_v03 \
    --condition v0.3_lba1 \
    --workers 4 \
    > logs/eval_additivity_v03.log 2>&1

python3 analysis/eval_additivity_on_replications.py \
    --results_dir results/results_nelsi_ics_full_v06 \
    --condition v0.6_lba1 \
    --workers 4 \
    > logs/eval_additivity_v06.log 2>&1

echo "=== Step 2: Run new Nextflow benchmarks (v=0.0, LBA 3x, LBA 5x) + Additivity evaluation ==="

echo "--- [2a] NELSI v=0.0, LBA=1x ---"
nextflow run main.nf \
    -c next_configs/nextflow_nelsi_ics_full_v00.config \
    -work-dir work_nelsi_ics_full_v00 \
    -profile supercomputer \
    -resume \
    > logs/nextflow_nelsi_ics_full_v00.log 2>&1

python3 analysis/eval_additivity_on_replications.py \
    --results_dir results/results_nelsi_ics_full_v00 \
    --condition v0.0_lba1 \
    --workers 4 \
    > logs/eval_additivity_v00.log 2>&1

echo "--- [2b] NELSI v=0.3, LBA=3x ---"
nextflow run main.nf \
    -c next_configs/nextflow_nelsi_ics_full_lba3.config \
    -work-dir work_nelsi_ics_full_lba3 \
    -profile supercomputer \
    -resume \
    > logs/nextflow_nelsi_ics_full_lba3.log 2>&1

python3 analysis/eval_additivity_on_replications.py \
    --results_dir results/results_nelsi_ics_full_lba3 \
    --condition v0.3_lba3 \
    --lba_pair T1 T20 \
    --workers 4 \
    > logs/eval_additivity_lba3.log 2>&1

echo "--- [2c] NELSI v=0.3, LBA=5x ---"
nextflow run main.nf \
    -c next_configs/nextflow_nelsi_ics_full_lba5.config \
    -work-dir work_nelsi_ics_full_lba5 \
    -profile supercomputer \
    -resume \
    > logs/nextflow_nelsi_ics_full_lba5.log 2>&1

python3 analysis/eval_additivity_on_replications.py \
    --results_dir results/results_nelsi_ics_full_lba5 \
    --condition v0.3_lba5 \
    --lba_pair T1 T20 \
    --workers 4 \
    > logs/eval_additivity_lba5.log 2>&1

echo "=== All ICS Full x NELSI x LBA evaluations completed ==="
