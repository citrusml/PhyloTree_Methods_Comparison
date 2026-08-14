#!/bin/bash
#BSUB -J nextflow_main
#BSUB -q mafft
#BSUB -n 1
#BSUB -M 8GB
#BSUB -W 24:00
#BSUB -o /lustre10/home/citrusml2004/PhyloTree_Methods_Comparison/logs/nextflow_lsf_%J.log
#BSUB -e /lustre10/home/citrusml2004/PhyloTree_Methods_Comparison/logs/nextflow_lsf_%J.err
source ~/.bashrc
micromamba activate phylomethod_env
cd /lustre10/home/citrusml2004/PhyloTree_Methods_Comparison
# ロックファイルクリーンアップ
rm -f .nextflow/cache/*/db/LOCK
nextflow run main.nf -profile supercomputer -resume -ansi-log false