#!/bin/bash
#BSUB -J nextflow_ics_full
#BSUB -q mafft
#BSUB -n 1
#BSUB -M 16GB
#BSUB -W 120:00
#BSUB -o logs/nextflow_ics_full_bsub.out
#BSUB -e logs/nextflow_ics_full_bsub.err

# 環境設定
export PATH="$HOME/bin:$PATH"
eval "$(micromamba shell hook --shell bash)"
micromamba activate phylomethod_env

cd /lustre10/home/citrusml2004/PhyloTree_Methods_Comparison
mkdir -p logs results/results_ics_full

# JVMヒープを明示的に指定
export NXF_OPTS="-Xms2g -Xmx12g"

# 実験17: 全長 ICS モデル (Dayhoff 6クラス不変) + Indel ベンチマーク
#   - 置換モデル: ICS+G4 (alpha=1.0)
#   - シミュレータ: AliSim (--mdef models/ics_model.nex)
#   - 系統樹生成モデル: 後退 Yule 過程 (paper_yule)
#   - Indel: Zipfian べき乗則 (insert=0.10, delete=0.10, POW{1.7/50})
#   - 比較手法: PWA+NJ, MSA+NJ, MSA+ML, MSA+RAXML, GS, TRUE対照群, SSS
echo "[$(date)] 実験17 (ICS Full Benchmark) 開始"
nextflow run main.nf \
    -c next_configs/nextflow_ics_full.config \
    -profile supercomputer \
    -resume \
    > logs/nextflow_ics_full.log 2>&1

echo "[$(date)] 実験17 完了 (exit code: $?)"
