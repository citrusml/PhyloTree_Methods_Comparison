#!/bin/bash
# ==============================================================================
# ローカル動作検証用テスト実行スクリプト
# ==============================================================================
# 実験15の設定（next_configs/nextflow_paper_tree.config）の構文および
# パイプライン全モジュール（SIMULATE, SSS, PWA, MAFFT, NJ, ML, RAxML, GS, 評価）
# が正常に動作するかを小規模条件（8 taxa, 2 distances, 2 reps）で高速に検証します。

set -euo pipefail

eval "$(micromamba shell hook --shell bash)"
micromamba activate phylomethod_env

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
cd "$PROJECT_ROOT"

echo "[$(date)] ローカル結合テスト開始 (-profile test)"
nextflow run main.nf \
    -c next_configs/nextflow_paper_tree.config \
    -profile test \
    --outdir results/test_local_run

echo "[$(date)] ローカル結合テスト完了"
