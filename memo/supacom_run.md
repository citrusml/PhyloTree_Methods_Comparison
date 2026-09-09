# スパコン（mafftsrv）実行ガイド

## 前提条件

- スパコン: `mafftsrv1` 〜 `mafftsrv5`（計 160 CPUコア / LSF スケジューラ）
- リポジトリ: `https://github.com/citrusml/PhyloTree_Methods_Comparison`
- パッケージマネージャ: `micromamba`

---

## Step 0: 初回セットアップ（初回のみ）

### 0-1. リポジトリのクローン

```bash
cd ~
git clone https://github.com/citrusml/PhyloTree_Methods_Comparison.git
cd PhyloTree_Methods_Comparison
```

### 0-2. micromamba のインストール（未インストールの場合）

```bash
curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj bin/micromamba
mkdir -p ~/bin
mv bin/micromamba ~/bin/micromamba
export PATH="$HOME/bin:$PATH"

# ~/.bashrc に追加して永続化
echo 'export PATH="$HOME/bin:$PATH"' >> ~/.bashrc
source ~/.bashrc

# バージョン確認
micromamba --version
```

### 0-3. 解析環境の構築

```bash
cd ~/PhyloTree_Methods_Comparison
micromamba create -f environment.yml -y
```

> ✅ `phylomethod_env` 環境が作成されます（mafft, iqtree, rapidnj, fastme, nextflow, python 等を含む）

### 0-4. Nextflow の動作確認

```bash
micromamba run -n phylomethod_env nextflow -version
```

---

## Step 1: 動作確認テスト（ローカル小規模実行）

スパコンにログイン後、本番実行前に小規模テストを行います。

```bash
cd ~/PhyloTree_Methods_Comparison
micromamba activate phylomethod_env

# ローカルプロファイルで小規模テスト（D=2条件 × L=2条件 × 5回反復）
nextflow run main.nf -profile local
```

> 成功すると `results/benchmark_summary.csv` と `results/regime_map_delta_nrf.png` が生成されます。

---

## Step 2: 本番実行（Phase 1 全パラメータ）

### 2-1. 実行コマンド（LSF バッチジョブ - 推奨）

```bash
cd ~/PhyloTree_Methods_Comparison
mkdir -p logs

# LSF バッチジョブとして投入（SSH切断の影響を受けず安定実行）
bsub < run_nextflow_job.sh
```

またはターミナルで直接バックグラウンド実行する場合：
```bash
micromamba activate phylomethod_env
mkdir -p logs

# 進捗ログを logs/run_nextflow.log に出力しながら実行
nohup nextflow run main.nf -profile supercomputer -resume > logs/run_nextflow.log 2>&1 < /dev/null &
```

### 2-2. 実行内容

| パラメータ | 値 |
|---|---|
| 進化距離 D | 0.1, 0.5, 1.0, 2.0, 3.0（5条件） |
| 配列長 L | 100, 300, 500, 1000, 1500 aa（5条件） |
| 試行回数 | 100 replicates |
| タキソン数 | N = 16 |
| パイプライン | PWA+NJ, MSA+NJ, MSA+ML（各 2,500 ジョブ） |
| **総 LSF ジョブ数** | **約 10,000 ジョブ** |
| **推定壁時計時間** | **約 6〜8 時間** |

### 2-3. LSF リソース割り当て（自動）

| プロセス | CPU | メモリ | 時間上限 |
|---|---|---|---|
| SIMULATE_DATA | 1 | 2 GB | 30 min |
| RUN_PWA_NJ | 1 | 4 GB | 1 h |
| RUN_MSA_NJ | 2 | 4 GB | 1 h |
| RUN_MSA_ML | 4 | 8 GB | **4 h** |
| COLLECT_AND_PLOT | 4 | 16 GB | 1 h |

---

## Step 3: 実行状況・進捗ログの確認

### Nextflow 進捗ログのリアルタイム確認

Nextflow の標準プログレス表示（各プロセスの完了率や実行状況）をリアルタイムで確認できます：

```bash
# bsub で投入した場合の進捗ログ確認
tail -f logs/nextflow_lsf_*.log

# nohup で実行した場合の進捗ログ確認
tail -f logs/run_nextflow.log
```

**出力例**:
```text
executor >  lsf (160)
[a1/b2c3d4] process > SIMULATE_DATA (D=0.1_L=100_rep=1) [100%] 2500 of 2500 ✔
[e5/f6g7h8] process > RUN_PWA_NJ (D=0.1_L=100_rep=1)    [ 40%] 1000 of 2500
[i9/j0k1l2] process > RUN_MSA_NJ (D=0.1_L=100_rep=1)    [ 40%] 1000 of 2500
[m3/n4o5p6] process > RUN_MSA_ML (D=0.1_L=100_rep=1)    [ 20%]  500 of 2500
```

### LSF ジョブ一覧

```bash
bjobs
```

### 終了後の出力確認

```bash
ls results/
# benchmark_summary.csv  - 全実験結果の集計表
# regime_map_delta_nrf.png - 相図ヒートマップ
# replications/          - 各レプリケーションごとの木、アライメント、距離行列
```

---

## Step 4: 結果の取得

スパコンからローカルへ結果を転送します（ローカル端末で実行）：

```bash
scp -r citrusml2004@mafftweb:~/PhyloTree_Methods_Comparison/results/ ./results_phase1/
```

---

## トラブルシューティング

### 途中でジョブが失敗した場合

Nextflow は自動でリジュームが可能です：

```bash
# -resume オプションで完了済みジョブをスキップして再実行
nextflow run main.nf -profile supercomputer -resume
```

### LSF キュー名（`mafft`）について

`nextflow.config` にスパコンの標準キュー **`queue = 'mafft'`** を設定済みです。通常は何も意識せずそのまま実行できます。必要に応じて `bqueues` で確認の上、別のキューを指定することも可能です：

```bash
# 利用可能なキュー一覧の確認
bqueues

# 別のキューを指定して実行したい場合
nextflow run main.nf -profile supercomputer --process.queue admin -resume
```

### メモリ不足エラーが出た場合

`nextflow.config` の `withName: RUN_MSA_ML { memory = '16 GB' }` に増量して再実行してください。

---

## ファイル構成（参考）

```
PhyloTree_Methods_Comparison/
├── main.nf               # Nextflow メインワークフロー
├── nextflow.config       # LSF クラスター設定
├── environment.yml       # micromamba 環境定義
├── bin/
│   ├── simulate_data.py  # 樹形・配列シミュレーション（AliSim / LG+G）
│   ├── run_pwa_nj.py     # PWA+NJ パイプライン（LG スコア行列 + Poisson 距離 + RapidNJ）
│   ├── run_msa_nj.py     # MSA+NJ 対照群（MAFFT + Poisson 距離 + RapidNJ）
│   ├── run_msa_ml.py     # MSA+ML パイプライン（MAFFT + IQ-TREE 2 ModelFinder）
│   ├── evaluate_trees.py # RF 距離・nRF 評価（DendroPy）
│   └── plot_regime_map.py# 相図ヒートマップ自動生成（Seaborn）
└── results/              # 出力ディレクトリ（実行後に生成）
    ├── benchmark_summary.csv
    └── regime_map_delta_nrf.png
```
