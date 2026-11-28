# グラフ作成用データ所在・欠損まとめ

## 共通デフォルト条件（基準条件）

- **進化距離 ($D$)**: `0.1, 0.2, 0.5, 1.0, 1.5, 2.0, 3.0`
- **Taxa 数 ($N$)**: `32`
- **配列長 ($L$)**: `500`
- **Indel 率 / 時計モデル / LBA**: `indel = 0.1`, `nelsi v = 0.6`, `lba_ratio = 1.0` (`1*`), `100 rep`
- **対象手法（7項目）**:
  1. `MSA+ML`（IQ-TREE 2、モデル推定 / `MFP`）
  2. `MSA+ML`（IQ-TREE 2、モデル指定: 正解モデル `WAG+G4`）※ FastTree ではなく **IQ-TREE 2 が必須**
  3. `MSA+NJ`（`poisson`）
  4. `PSA+NJ`（`poisson`、`BLOSUM62, open = 10, extend = 1`）
  5. `PSA+NJ`（`sss linear`）
  6. `GS`（`poisson`）
  7. `NJ + FastME (NNI+SPR)`（上記 `MSA+NJ` / `PSA+NJ` の距離行列に対する FastME の NNI+SPR 探索）

> **方針**: Taxa 数・配列長・NELSI $v$・LBA の各比較では、上記のデフォルト条件から**対象パラメータのみを変更したデータ**を使用する。

---

## 軸別 データ所在・欠損状況

### 1. デフォルト条件（$N=32, L=500, v=0.6, \mathrm{LBA}=1\times$）
- **存在するデータ（4手法）**:
  - レプリケート生データ: `/Users/aibarakazuki/JoiningBasedOnAdditivity/benchmark_results.csv`（`n_taxa == 32, nelsi_v == 0.6, indel_rate == 0.1` で 9,800 行）
  - 集計サマリー: `/Users/aibarakazuki/JoiningBasedOnAdditivity/results/nelsi_v06_v10/summary_by_condition.csv`
  - 揃っている手法:
    - `MSA+NJ (poisson)` → `MSA_NJ`
    - `PSA+NJ (poisson, BLOSUM62 open=10 extend=1)` → `Poisson_NJ`
    - `PSA+NJ (sss linear)` → `SSS_Linear_NJ`
    - `GS (poisson)` → `GS_Poisson`
- **存在しないデータ（3手法）**:
  - `MSA+ML`（IQ-TREE 2: `MFP`）および `MSA+ML`（IQ-TREE 2: `WAG+G4`）※ `MSA_ML` は FastTree `-wag -gamma` のみ
  - `NJ + FastME (NNI+SPR)`（全距離行列）※ `SSS_Linear_FastME_NNI`（NNI のみ）しか存在しない

---

### 2. Taxa 数スイープ（$N \in \{16, 32, 64, 128, 256\}$, $L=500, v=0.6, \mathrm{LBA}=1\times$）

#### 存在するデータ（$N = 32, 64$ の 4 手法）
- **所在**:
  - `/Users/aibarakazuki/JoiningBasedOnAdditivity/benchmark_results.csv`（`nelsi_v == 0.6, indel_rate == 0.1`, $N \in \{32, 64\}$ で計 19,600 行）
  - `/Users/aibarakazuki/JoiningBasedOnAdditivity/results/nelsi_v06_v10/summary_by_condition.csv`
- **揃っている手法**: `MSA_NJ`, `Poisson_NJ`, `SSS_Linear_NJ`, `GS_Poisson`

#### 存在しないデータ
| 対象範囲 | 欠損しているデータ | 備考 |
| :--- | :--- | :--- |
| **$N = 16, 128, 256$** | **全 7 手法・全距離** | データ自体が未存在 |
| **$N = 32, 64$** | `MSA+ML`（IQ-TREE 2: `MFP` / `WAG+G4`） | FastTree `-wag -gamma` のみ |
| **$N = 32, 64$** | `NJ + FastME (NNI+SPR)`（全距離行列） | `SSS_Linear_FastME_NNI`（NNI のみ）しか存在しない |

---

### 3. 配列長スイープ（$L \in \{300, 500, 1000, 1500\}$, $N=32, v=0.6, \mathrm{LBA}=1\times$）
- **存在するデータ**:
  - **$L = 500$**: デフォルト条件と同一（4 手法が `JoiningBasedOnAdditivity` に存在、IQ-TREE 2 `MSA+ML` と `FastME (NNI+SPR)` が欠損）
- **存在しないデータ**:
  - **$L = 300, 1000, 1500$**: **全 7 手法・全距離**（$v = 0.6$ は $L = 500$ のみ存在）

---

### 4. NELSI $v$ スイープ（$v \in \{0.0, 0.3, 0.6, 1.0\}$, $N=32, L=500, \mathrm{LBA}=1\times$）
- **存在するデータ（$v = 0.0, 0.3, 0.6, 1.0$ 全 4 値・4 手法）**:
  - **所在**: `/Users/aibarakazuki/JoiningBasedOnAdditivity/benchmark_results.csv`（`n_taxa == 32, indel_rate == 0.1`, $v \in \{0.0, 0.3, 0.6, 1.0\}$ で各 `v` 9,800 行、計 39,200 行）および `summary_by_condition.csv`
  - **揃っている手法**: `MSA_NJ`, `Poisson_NJ`, `SSS_Linear_NJ`, `GS_Poisson`
- **存在しないデータ（$v = 0.0, 0.3, 0.6, 1.0$ 共通）**:
  - `MSA+ML`（IQ-TREE 2: `MFP` / `WAG+G4`）※ FastTree `-wag -gamma` のみ
  - `NJ + FastME (NNI+SPR)`（全距離行列）※ `SSS_Linear_FastME_NNI`（NNI のみ）しか存在しない
  - （参考）本リポジトリ `results/results_nelsi_taxa32_L500` には $v = 0.3$ の IQ-TREE 2 `WAG+G4` 固定版があるが、$D = 3.0$ なし・`GS` なし・`gap 20/2`・`MSA+NJ = wag_gamma`。

---

### 5. LBA スイープ（$\mathrm{LBA} \in \{1\times, 3\times, 5\times\}$, $N=32, L=500, v=0.6$）
- **存在するデータ**:
  - **$\mathrm{LBA} = 1\times$**: デフォルト条件と同一（4 手法が `JoiningBasedOnAdditivity` に存在、IQ-TREE 2 `MSA+ML` と `FastME (NNI+SPR)` が欠損）
- **存在しないデータ**:
  - **$\mathrm{LBA} = 3\times, 5\times$**: **全 7 手法・全距離**（本リポジトリの `results_nelsi_lbs3/5` は $N=20, v=0.3$ のため非該当）

---

## 新規データ収集用ファイル・スクリプト一覧

### 1. コピー済み既存データ（`graph/`）
- `graph/existing_joining_benchmark_results.csv`（`JoiningBasedOnAdditivity/benchmark_results.csv` のコピー）
- `graph/existing_joining_summary_by_condition.csv`（`JoiningBasedOnAdditivity/results/nelsi_v06_v10/summary_by_condition.csv` のコピー）

### 2. Nextflow 設定ファイル（`next_configs/graph_benchmark/`）
- `next_configs/graph_benchmark/base_small_medium.config`: $N \le 64$ 用デフォルト設定（$D \in [0.1..3.0]$, $N=32$, $L=500$, `indel=0.1`, `nelsi_v=0.6`, `lba_ratio=1.0`, `100 rep`, 全 7 手法有効、`BLOSUM62 open=10 extend=1`, `nj_spr=true`）
- `next_configs/graph_benchmark/large_taxa_128_256.config`: $N=128, 256$ 用設定（`chunk_size=2`、メモリ・実行時間枠を拡張）

### 3. 実行スクリプト（`run_scripts/graph_benchmark/`）
※ デフォルト点（$N=32, L=500, v=0.6, \mathrm{LBA}=1\times$）の出力先を `results/graph_benchmark/default_N32_L500_v06_lba1` に統一しているため、`run_default_N32_L500_v06.sh` を先に実行すると各スイープの `-resume` で重複計算をスキップできます。

| スクリプト名 | 対象条件 | 主な出力先 (`results/graph_benchmark/...`) |
| :--- | :--- | :--- |
| `run_scripts/graph_benchmark/run_default_N32_L500_v06.sh` | デフォルト（$N=32, L=500, v=0.6, \mathrm{LBA}=1\times$） | `default_N32_L500_v06_lba1` |
| `run_scripts/graph_benchmark/run_taxa_16_32_64.sh` | Taxa スイープ前半（$N \in \{16, 32, 64\}$） | `taxa_N16_L500_v06_lba1`, `default_...`, `taxa_N64_L500_v06_lba1` |
| `run_scripts/graph_benchmark/run_taxa_128_256.sh` | Taxa スイープ後半（$N \in \{128, 256\}$, `chunk_size=2, 1`） | `taxa_N128_L500_v06_lba1`, `taxa_N256_L500_v06_lba1` |
| `run_scripts/graph_benchmark/run_lengths_300_500_1000_1500.sh` | 配列長スイープ（$L \in \{300, 1000, 1500\}$、$L=500$ は `run_taxa_16_32_64.sh` で実行） | `length_N32_L{300,1000,1500}_v06_lba1` |
| `run_scripts/graph_benchmark/run_nelsi_v_00_03_06_10.sh` | NELSI $v$ スイープ（$v \in \{0.0, 0.3, 1.0\}$、$v=0.6$ は `run_taxa_16_32_64.sh` で実行） | `nelsi_N32_L500_v{00,03,10}_lba1` |
| `run_scripts/graph_benchmark/run_lba_1_3_5.sh` | LBA スイープ（$\mathrm{LBA} \in \{3\times, 5\times\}$、$1\times$ は `run_taxa_16_32_64.sh` で実行） | `lba_N32_L500_v06_lba{3,5}` |

### 4. パイプライン側の拡張（FastME NNI+SPR 対応）
- `bin/run_nj_nni.py`（`--spr` フラグ追加で `fastme -n B -s` を実行）
- `modules/run_nj_nni.nf`（`params.nj_spr = true` 時に `--spr` を付与）
- `main.nf`（`params.nj_spr = true` 時に出力パイプライン名を `PWA+NJ+NNI+SPR`, `MMSEQS+NJ+NNI+SPR`, `MSA+NJ+NNI+SPR` に設定）

