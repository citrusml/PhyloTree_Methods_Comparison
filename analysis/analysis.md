
## 方針

解析は基本的にipynbファイルでlocalで行う。
まずは今までの結果*_summary.csvをまとめる。これまでの実験では変数がそれぞれ異なったので*_summary.csvのカラムが異なっているが、それを統一し、一つのcsvファイルにまとめる。

そしてそれぞれの条件でのsimulationを1つだけ実行し、どのようなmsa, treeができているかを確認する。このmsaやtreeはanalysisのexampleフォルダに条件ごとに入れること。

### 実験1（基本パラメータ空間: 距離 $D$ $\times$ 配列長 $L$）
- Taxon 数 $N = 32$
- 進化距離 $D \in [0.1, 0.5, 1.0, 2.0, 3.0]$
- 配列長 $L \in [100, 300, 500, 1000, 1500]$
- Replicates: 100
- 置換モデル: `LG+G4`
- $\alpha = 1.0$ (固定)
- 実行パイプライン: `PWA+NJ`, `MSA+NJ`, `MSA+ML`, `TRUE_DIST`, `TRUE_MSA+NJ`, `TRUE_MSA+ML`

### 実験2（Taxon数 スケーリング実験）
- Taxon 数 $N \in [8, 16, 64, 128]$
- 進化距離 $D \in [0.1, 0.5, 1.0, 2.0, 3.0]$
- 配列長 $L \in [100, 500, 1000]$
- Replicates: 100
- 置換モデル: `LG+G4`
- $\alpha = 1.0$ (固定)
- 実行パイプライン: `PWA+NJ`, `MSA+NJ`, `MSA+ML`

### 実験3（ガンマ形状母数 $\alpha$ / サイト間速度不均一性実験）
- Taxon 数 $N = 32$
- 進化距離 $D \in [0.1, 0.5, 1.0, 2.0, 3.0]$
- 配列長 $L \in [100, 500, 1000]$
- ガンマ形状母数 $\alpha \in [0.25, 0.5, 1.0, 2.0]$
- Replicates: 100
- 置換モデル: `LG+G4`
- 実行パイプライン: `PWA+NJ`, `MSA+NJ`, `MSA+ML`

### 実験4（真のペアワイズアライメント TRUE_PWA+NJ 実験）
- Taxon 数 $N = 32$
- 進化距離 $D \in [0.1, 0.5, 1.0, 2.0, 3.0]$
- 配列長 $L \in [100, 300, 500, 1000, 1500]$
- Replicates: 100
- $\alpha = 1.0$ (固定)
- 置換モデル: `LG+G4`
- 実行パイプライン: `TRUE_PWA+NJ`

### 実験5（ガンマ補正距離 Gamma Distance Benchmark）
- Taxon 数 $N = 32$
- 進化距離 $D \in [0.1, 0.5, 1.0, 2.0, 3.0]$
- 配列長 $L \in [100, 500, 1000]$
- ガンマ形状母数 $\alpha \in [0.25, 0.5, 1.0, 2.0]$
- 距離計算モデル: `gamma_poisson`（$\alpha_{\text{dist}} = 1.0$ 固定）
- Replicates: 100
- 置換モデル: `LG+G4`
- 実行パイプライン: `PWA+NJ`, `MSA+NJ`, `MSA+ML`

#### 4. 実験6（ICS条件ベンチマーク）の実験設定案
- Taxon 数: $N = 32$
- 進化距離: $D \in [0.1, 0.5, 1.0, 2.0, 3.0]$
- 配列長: $L \in [100, 300, 500, 1000]$
- ICS 比率: $\mathrm{ics\_prop} \in [0.0, 0.05, 0.1, 0.2]$
- ガンマ形状母数: $\alpha = 1.0$ (固定)
- Replicates: 100
- 実行パイプライン: `PWA+NJ`, `MSA+NJ`, `MSA+ML`

### 実験7（FastME 距離法ベンチマーク: 距離 $D$ $\times$ 配列長 $L$）
- Taxon 数: $N = 32$
- 進化距離: $D \in [0.1, 0.5, 1.0, 2.0, 3.0]$
- 配列長: $L \in [100, 300, 500, 1000, 1500]$
- ガンマ形状母数: $\alpha = 1.0$ (固定)
- 置換モデル: `LG+G4`
- 距離モデル: Poisson 距離
- 系統樹構築ツール: `fastme`
- Replicates: 100
- 実行パイプライン: `PWA+FastME`, `MSA+FastME`  <- これの名前を`FastME`を`FastME_NoOption`にしておいて

### 実験8（高進化距離領域ベンチマーク: 距離 $D \in [4.0, 5.0, 6.0]$ $\times$ 配列長 $L$）
- Taxon 数: $N = 32$
- 進化距離: $D \in [4.0, 5.0, 6.0]$
- 配列長: $L \in [300, 500, 1000, 1500]$
- ガンマ形状母数: $\alpha = 1.0$ (固定)
- 置換モデル: `LG+G4`
- Replicates: 100
- 実行パイプライン: `PWA+NJ`, `MSA+NJ`, `MSA+ML`

### 実験9（LGmodelを用いる）
- Taxon 数: $N = 32$
- 進化距離: $D \in [0.1, 0.5, 1.0, 2.0, 3.0]$
- 配列長: $L \in [100, 300, 500, 1000, 1500]$
- ガンマ形状母数: $\alpha = 1.0$ (固定)
- 置換モデル: `LG`
- Replicates: 100
- 実行パイプライン: `PWA+NJ`, `MSA+NJ`, `MSA+ML`

### 実験10（FastME オプション付き: LG+Gamma & SPR 最適化）
- Taxon 数: $N = 32$
- 進化距離: $D \in [0.1, 0.5, 1.0, 2.0, 3.0]$
- 配列長: $L \in [100, 300, 500, 1000, 1500]$
- ガンマ形状母数: $\alpha = 1.0$ (固定)
- 置換モデル: `LG+G4`
- Replicates: 100
- 実行パイプライン:
  - `MSA+FastME_LG_G`: MAFFT MSA $\to$ `fastme -i msa.phy -pL -g1.0 -s -q` (FastME 内蔵 LG+G 距離 + SPR トポロジー改善 + 三角不等式補正)
  - `PWA+FastME_SPR`: PWA Poisson 距離行列 $\to$ `fastme -i pwa_matrix.phy -s -q` (PWA 距離行列に対する SPR トポロジー改善 + 三角不等式補正)

### 実験11（全長 ICS モデル＋Indel: Full-length Invariant Category Sites with Indels）
- Taxon 数: $N = 32$
- 進化距離: $D \in [0.1, 0.5, 1.0, 2.0, 3.0]$
- 配列長: $L \in [100, 300, 500, 1000, 1500]$
- ICS 比率: $\mathrm{ics\_prop} = 1.0$ (全長 100% Dayhoff 6分類不変サイトモデル)
- 挿入率: $\mathrm{insert\_rate} = 0.05$, 欠失率: $\mathrm{delete\_rate} = 0.10$ (Indel 許容)
- ガンマ形状母数: $\alpha = 1.0$ (固定)
- 置換モデル: `ICS+G4`
- Replicates: 100
- 実行パイプライン: `PWA+NJ` (PSA+NJ), `MSA+NJ`, `MSA+ML`
- 実行スクリプト: `run_scripts/run_nextflow_ics_full.sh` (出力先: `results/results_ics_full`)



