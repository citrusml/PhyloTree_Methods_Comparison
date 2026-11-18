# Phylogenetic Benchmark Results: Simulator Evolution & SSS-Based Evaluation

本ドキュメントは、本リポジトリの [`results/`](./results) 配下に蓄積された系統樹推定ベンチマーク実験（[`analysis/all_experiments_summary.csv`](./analysis/all_experiments_summary.csv) 統合済み全 23 実験・総計 **266,740 系統樹推定**、および外部比較として [`JoiningBasedOnAdditivity`](../JoiningBasedOnAdditivity/results_hpc_indel0.1/main/hpc_benchmark_summary.csv) の BioPerl Taxa スケーリング結果）について、**Sequence Similarity Score ($SSS$, $w$)**・**進化距離 ($D$)**・**配列長 ($L$)**・**Taxa数 ($N$)** を軸としたトポロジー誤差（$nRF$）の解析結果、および **シミュレーターの変遷による影響** をまとめたレポートです。

---

## 1. 全実験の概括比較テーブル

| # | 実験ディレクトリ | 系統樹モデル | 配列・Indel シミュレータ | 主な検証パラメータ | 比較パイプライン | SSS 範囲 $[w_{\min} .. w_{\max}]$ | 主な発見・挙動 |
|---|---|---|---|---|---|---:|---|
| 1 | [`results_zipfian_indel`](./results/results_zipfian_indel) | Birth-Death | AliSim (`LG+G4`, `POW{1.7/50}`) | $L \in [300..1500]$, $D \in [0.1..6.0]$ | MSA+ML, MSA+NJ, PSA+NJ, TRUE群 | $[0.0078 .. 0.9318]$ | 長大インデルにより $w < 0.06$ で MSA+ML が急激に崩壊し、PSA+NJ が逆転 |
| 2 | [`results_paper_tree_pre_fix`](./results/results_paper_tree_pre_fix) | 後退 Yule (BioPerl) | AliSim (`WAG+G4`, `POW{1.7/50}`) | $N=32$, $L \in [500..1500]$, $D \in [0.1..2.0]$ | GS, MSA+ML, MSA+RAXML, PSA+NJ, TRUE群 | $[0.0054 .. 0.6745]$ | 論文系統樹モデル（初期版）。超低類似度域での GS 法の安定性を確認 |
| 3 | [`results_alpha`](./results/results_alpha) | Birth-Death | AliSim (`LG+G4`, 幾何 Indel) | $\alpha \in [0.25..2.0]$, $L \in [100..1000]$ | MSA+ML, MSA+NJ, PSA+NJ | $[0.0063 .. 0.9537]$ | $\alpha=0.25$（極度の速度不均一）で SSS 減衰が最も激しく、ML の $+G4$ 適合が有効 |
| 4 | [`results_gamma`](./results/results_gamma) | Birth-Death | AliSim (`LG+G4`, 幾何 Indel) | $\alpha \in [0.25..2.0]$, `gamma_poisson` | MSA+ML, MSA+NJ, PSA+NJ | $[0.0063 .. 0.9537]$ | ガンマ補正距離により高距離域での距離過小評価が緩和 |
| 5 | [`results_taxon`](./results/results_taxon) | Birth-Death | AliSim (`LG+G4`, 幾何 Indel) | $N \in [8, 16, 32, 64, 128]$, $L \in [100..1000]$ | MSA+ML, MSA+NJ, PSA+NJ | $[0.0000 .. 0.9861]$ | $N$ 増加に伴い探索空間が拡大し、低 SSS 側での $nRF$ 悪化が加速 |
| 6 | [`results_ics`](./results/results_ics) | Birth-Death | AliSim (`LG+G4_ICS`, 幾何 Indel) | $ICS \in [0.0..0.20]$, $L \in [100..1000]$ | MSA+ML, MSA+NJ, PSA+NJ | $[0.0072 .. 0.9525]$ | 保存ブロック（ICS）が局所整列のアンカーとなり、遠縁領域の SSS を底上げ |
| 7 | [`results_ics_full`](./results/results_ics_full) | Birth-Death / 後退 Yule | AliSim (`ICS+G4`, $ICS=1.0$) | $L \in [100..1500]$, $D \in [0.1..3.0]$ | MSA+ML, MSA+NJ, MSA+RAXML, PSA+NJ | $[0.0145 .. 0.9451]$ | 全長保存ドメイン下では変異速度低下により全距離帯で極めて高い精度を維持 |
| 8 | [`results_high_dist`](./results/results_high_dist) | Birth-Death | AliSim (`LG+G4`, 幾何 Indel) | 極限距離 $D \in [4.0, 5.0, 6.0]$ | MSA+ML, MSA+NJ, PSA+NJ | $[0.0078 .. 0.0910]$ | 全反復が Twilight Zone ($w < 0.10$) に集中し、情報限界に到達 |
| 9 | [`results_high_gap`](./results/results_high_gap) | Birth-Death | AliSim (`LG+G4`, Gap 20/2) | $D \in [0.1..6.0]$, FastME 比較 | MSA+ML, MSA+NJ, PSA+NJ, FastME群 | $[0.0078 .. 0.9318]$ | ギャップペナルティ引き上げ後も、幾何 Indel 下では ML 優位を維持 |
| 10 | [`results_simple`](./results/results_simple) | Birth-Death | AliSim (`LG`, 均一速度) | $L \in [100..1500]$, $D \in [0.1..3.0]$ | MSA+ML, MSA+NJ, PSA+NJ | $[0.0035 .. 0.9537]$ | サイト間速度変動がない場合、高類似度域の推定精度が全体的に向上 |
| 11 | [`results_new`](./results/results_new) | Birth-Death | AliSim (`LG+G4`, 幾何 Indel) | $D \in [0.1..3.0]$, 真の距離行列対照 | MSA+ML, MSA+NJ, PSA+NJ, `TRUE_DIST+NJ` | $[0.0064 .. 0.9494]$ | `TRUE_DIST+NJ` によりアライメント誤差と距離推定誤差の寄与を分離評価 |
| 12 | [`results_fastme_NoOption`](./results/results_fastme_NoOption) | Birth-Death | AliSim (`LG+G4`, 幾何 Indel) | FastME デフォルト vs RapidNJ | MSA+FastME, PSA+FastME | $[0.0064 .. 0.9494]$ | 長大配列で FastME が NJ より良好なトポロジーを構築 |
| 13 | [`results_fastme_options`](./results/results_fastme_options) | Birth-Death | AliSim (`LG+G4`, 幾何 Indel) | FastME (SPR 探索 + LG+G 補正) | MSA+FastME_LG_G, PSA+FastME_SPR | $[0.0078 .. 0.9318]$ | SPR 局所探索の導入により低 SSS 域でのトポロジー探索能が向上 |
| 14 | [`results_true_pwa`](./results/results_true_pwa) | Birth-Death | AliSim (`LG+G4`, 幾何 Indel) | 真のペアワイズ整列対照 | `TRUE_PSA+NJ` | $[0.0064 .. 0.9494]$ | PSA におけるアラインメント誤差を除いた理論的上限値を画定 |
| 15 | [`results_paper_tree2`](./results/results_paper_tree2) | 後退 Yule (BioPerl) | **INDELible v1.03** (`WAG+G4`, `POW 1.7 50`) | $N=20, L=1000, D \in [0.1..2.0]$, 60 reps | GS, MSA+ML, MSA+RAXML, MSA+NJ, PSA+NJ, TRUE群 | $[0.0097 .. 0.6516]$ | **論文 Fig. 3a 完全再現**。$w \le 0.03$ で GS が ML/RAxML と同等以上の耐性を実証 |
| 16 | [`results_paper_tree3`](./results/results_paper_tree3) | 後退 Yule (BioPerl) | **INDELible v1.03** (`WAG+G4`, `POW 1.7 50`) | $N=20, L \in [500, 1000, 1500]$, 100 reps | GS, MSA+ML, MSA+RAXML, MSA+NJ, PSA+NJ (NW), TRUE群 | $[0.0086 .. 0.6516]$ | 16,200 系統樹推定。SP score を介した 4 段階誤差伝播（$D \to SP \to \text{MAE} \to nRF$）を解明 |
| 17 | [`results_paper_tree4`](./results/results_paper_tree4) | 後退 Yule (BioPerl) | **INDELible v1.03** (`WAG+G4`, `POW 1.7 50`) | $N=20, L \in [500, 1000, 1500]$, 100 reps | `PSA+NJ` (MMseqs2 Linear), `MMSEQS+NJ` (Log), `MSA+NJ` (WAG+G4), TRUE群 | $[0.0086 .. 0.6516]$ | **有界線形距離 $d=1-w_{ij}$** が NW 飽和と対数発散を克服し、$D=2.0$ で最良精度 |
| 18 | [`results_nelsi_yule`](./results/results_nelsi_yule) | **Yule + NELSI 自己相関** | **INDELible v1.03** (`WAG+G4`, `POW 1.7 50`) | $N=20, L \in [500..1500]$, $v=0.3$, root-to-tip=$D$ | GS, MMSEQS+NJ (log), MSA+ML/NJ/RAXML, PSA+NJ (sss linear), TRUE群 | $[0.0119 .. 0.9006]$ | 高 SSS 域で MSA+ML/RAXML と NJ (sss linear) が高精度。系統間速度差により `TRUE_MSA+NJ` は悪化 |
| 19 | [`results_nelsi_taxa32_L500`](./results/results_nelsi_taxa32_L500) | **Yule + NELSI 自己相関** | **INDELible v1.03** (`WAG+G4`, `POW 1.7 50`) | $N=32, L=500$, $v=0.3$, root-to-tip=$D$ | MMSEQS+NJ (sss linear), PSA+NJ (poisson), MSA+NJ/ML/RAXML, 各+NNI | $[0.0125 .. 0.8923]$ | `NJ (sss linear)` (`0.1325`) が `NJ (poisson)` (`0.1837`)・`MSA+NJ` (`0.1800`) を大きく上回る |
| 20 | [`results_nelsi_taxa64_L500`](./results/results_nelsi_taxa64_L500) | **Yule + NELSI 自己相関** | **INDELible v1.03** (`WAG+G4`, `POW 1.7 50`) | $N=64, L=500$, $v=0.3$, root-to-tip=$D$ | 同上（5手法 + NJ+NNI） | $[0.0111 .. 0.8608]$ | `MSA+ML` (`0.1214`)・`NJ (sss linear)` (`0.1458`) とも $N=64$ でも劣化極めて小 |
| 21 | [`results_nelsi_taxa128_L500`](./results/results_nelsi_taxa128_L500) | **Yule + NELSI 自己相関** | **INDELible v1.03** (`WAG+G4`, `POW 1.7 50`) | $N=128, L=500$, $v=0.3$, root-to-tip=$D$ | 同上（5手法 + NJ+NNI） | $[0.0088 .. 0.8299]$ | BioPerl と異なり `MSA+ML` (`0.1177`) が悪化せず維持。NNI は NJ と全 1,800 例で完全一致 |

---

## 2. カテゴリ別 詳細結果とプロット

### カテゴリ A: 論文再現・現実的 Indel モデル & 系統樹シミュレーター比較

#### A-1. `results_zipfian_indel` (AliSim + べき乗則 Indel 長分布)
- **条件**: べき指数 $a=1.7$・最大 50 残基の Zipfian インデル、Indel 率 $(0.05, 0.10)$、$L \in [300, 500, 1000, 1500]$、$D \in [0.1 .. 6.0]$
- **グラフ**: [レジームマップ (`regime_map_delta_nrf.png`)](./results/results_zipfian_indel/regime_map_delta_nrf.png) / [箱ひげ図 (`nrf_boxplots.png`)](./results/results_zipfian_indel/nrf_boxplots.png)

![results_zipfian_indel レジームマップ](./results/results_zipfian_indel/regime_map_delta_nrf.png)

- **要点**: 長大インデルにより多重配列アラインメント（MAFFT）のギャップ誤配置が急増し、$w < 0.06$ で `MSA+ML` が崩壊して `PSA+NJ` との逆転（Crossover）が発生します。

---

#### A-2. `results_paper_tree2` (INDELible v1.03 + 後退 Yule 系統樹による論文 Fig. 3a 再現)
- **条件**: INDELible v1.03（WAG+G4, $\alpha=1.0$, `POW 1.7 50`, 挿入/欠失率 0.10）、BioPerl 後退 Yule 系統樹、$N=20, L=1000, D \in [0.1..2.0]$, 60 反復（3,240 推定）
- **グラフ**: [SSS 感度曲線 (`nrf_vs_sss_curves.png`)](./results/results_paper_tree2/nrf_vs_sss_curves.png) / [レジームマップ (`regime_map_delta_nrf.png`)](./results/results_paper_tree2/regime_map_delta_nrf.png)

![results_paper_tree2 SSS 感度曲線](./results/results_paper_tree2/nrf_vs_sss_curves.png)

- **要点**:
  - Extreme 領域（$SSS \le 0.03$, 96 レプリケート）において、`GS`（平均 $nRF = 0.3768$）が `MSA+ML`（$0.3768$）および `MSA+RAXML`（$0.3811$）と同等以上の耐性を示し、論文 Fig. 3a の傾向を再現しました。
  - IQ-TREE 2 (`MSA+ML`) と RAxML (`MSA+RAXML`) は全距離帯でほぼ一致した感度曲線を描きます。

---

#### A-3. `results_paper_tree3` (100 反復大規模検証 & 4 段階誤差伝播カスケード解析)
- **条件**: `results_paper_tree2` と同一の INDELible + Zipfian + 後退 Yule 条件を $L \in [500, 1000, 1500]$・100 反復へ拡張（計 16,200 系統樹推定、SSS 範囲 $[0.0086 .. 0.6516]$、中央値 $0.0573$）。
- **グラフ**: [誤差伝播解析 (`error_propagation_analysis.png`)](./results/results_paper_tree3/error_propagation_analysis.png) / [SSS 感度曲線 (`nrf_vs_sss_curves.png`)](./results/results_paper_tree3/nrf_vs_sss_curves.png)

![results_paper_tree3 誤差伝播解析](./results/results_paper_tree3/error_propagation_analysis.png)

- **要点**:
  - 配列長別平均 $nRF$（$L=500 / 1000 / 1500$）は、`MSA+ML` が `0.1580 / 0.1393 / 0.1319`、`GS` が `0.1984 / 0.1768 / 0.1794`、`PSA+NJ (NW)` が `0.2004 / 0.1820 / 0.1757` となりました。
  - $D$ 増大 → SP スコア低下 → 距離行列 MAE 悪化 → $nRF$ 悪化の連鎖を定量化し、`GS` が SP 低下に対して最も緩やかな回帰傾きを持つことを示しました。

---

#### A-4. `results_paper_tree4` (MMseqs2 局所類似度 + 有界線形距離 & WAG+G4 ML 距離)
- **条件**: 原著論文（Matsui & Iwasaki 2020, p.268）の距離計算定義に完全統一。MMseqs2 (`-s 7.5 -e 10`) 相対類似度 $w_{ij}$ からの線形距離（`PWA+NJ`: $1-w_{ij}$）・対数距離（`MMSEQS+NJ`: $-\ln w_{ij}$）および `WAG+G4` 最尤距離（`MSA+NJ`）を評価（計 9,000 系統樹推定）。
- **グラフ**: [論文 Fig. 3a 統合プロット (`fig2_reproduction_fig3a_synthesis.png`)](./results/results_paper_tree4/plots/fig2_reproduction_fig3a_synthesis.png) / [線形 vs 対数直接比較 (`fig3_psa_linear_vs_log_comparison.png`)](./results/results_paper_tree4/plots/fig3_psa_linear_vs_log_comparison.png)

![results_paper_tree4 統合プロット](./results/results_paper_tree4/plots/fig2_reproduction_fig3a_synthesis.png)

- **要点**:
  - 配列長別平均 $nRF$（$L=500 / 1000 / 1500$）は、`PWA+NJ` (MMseqs2 Linear) が `0.1598 / 0.1386 / 0.1309`、`MSA+NJ` (WAG+G4) が `0.1694 / 0.1426 / 0.1372`、`MMSEQS+NJ` (Log) が `0.2124 / 0.1839 / 0.1688`。
  - NW 大域整列のギャップペナルティ飽和および対数変換 $-\ln w_{ij}$ の $w_{ij} \to 0$ 発散を、有界な $[0, 1]$ 線形距離が完全に抑え込み、遠縁領域（$D=1.5 \sim 2.0$）で `MSA+NJ` / `MSA+ML` を逆転しました。

---

#### A-5. `results_nelsi_yule` (Yule トポロジー + NELSI 自己相関緩和時計モデル)
- **条件**: トポロジーを純粋出生 Yule 過程（`ape::rphylo`）、枝長を **NELSI** の自己相関モデル（`simulate.autocor.kishino`, $v = 0.3$, `scale_mode = "root_to_tip"`）で生成し、INDELible v1.03 で配列進化（$N=20, D \in [0.1..2.0], L \in [500, 1000, 1500]$, 100 反復 $\times$ 10 パイプライン = 18,000 系統樹推定、SSS 範囲 $[0.0119 .. 0.9006]$、平均 `0.3178`、中央値 `0.2137`）。
- **グラフ**:
  - [SSS 感度曲線 (`nrf_vs_sss_by_condition.png`)](./results/results_nelsi_yule/nrf_vs_sss_by_condition.png)
  - [NELSI vs BioPerl 配列長別比較（エラーバー & バイオリンプロット）](./results/results_nelsi_yule/nelsi_vs_bioperl_by_length.png)

![results_nelsi_yule SSS 感度曲線](./results/results_nelsi_yule/nrf_vs_sss_by_condition.png)

![NELSI vs BioPerl 配列長別比較](./results/results_nelsi_yule/nelsi_vs_bioperl_by_length.png)

- **配列長別平均 $nRF$（$L = 500 / 1000 / 1500$）の対比**:
  - `MSA+ML`: NELSI `0.1097 / 0.0839 / 0.0697` vs BioPerl (PT3) `0.1580 / 0.1393 / 0.1319`
  - `PSA+NJ` (MMseqs2 Linear): NELSI `0.1232 / 0.0908 / 0.0775` vs BioPerl (PT4) `0.1598 / 0.1386 / 0.1309`
  - `MSA+NJ` (WAG+G4): NELSI `0.1651 / 0.1221 / 0.1061` vs BioPerl (PT4) `0.1694 / 0.1426 / 0.1372`
  - `GS`: NELSI `0.2058 / 0.1853 / 0.1675` vs BioPerl (PT3) `0.1984 / 0.1768 / 0.1794`
  - `TRUE_MSA+NJ`: NELSI **`0.1342 / 0.0944 / 0.0735`** vs BioPerl (PT4) **`0.0775 / 0.0509 / 0.0435`**
  - `TRUE_MSA+ML`: NELSI `0.0725 / 0.0501 / 0.0354` vs BioPerl (PT3) `0.0905 / 0.0492 / 0.0272`

---

#### A-6. NELSI Taxa数スケーリング（$L=500, N \in \{20, 32, 64, 128\}$）& BioPerl との比較解析
- **対象データ**: `results_nelsi_yule`（$N=20, L=500$ 抽出）、`results_nelsi_taxa32_L500`（$N=32$）、`results_nelsi_taxa64_L500`（$N=64$）、`results_nelsi_taxa128_L500`（$N=128$）。※注記: 距離モデルを厳密に揃えると、$N=20$ の `PSA+NJ` と $N=32..128$ の `MMSEQS+NJ` がともに **SSS 線形距離（`NJ (sss linear)`）**、$N=32..128$ の `PSA+NJ` が **ポアソン距離（`NJ (poisson)`）**、$N=20$ の `MMSEQS+NJ` が **`NJ (sss log)`** に対応します。
- **全距離（$D \in [0.1..2.0]$）統合平均 $nRF \pm SD$（$L=500$）**:
  - `MSA+ML`: `0.1097±0.1091` ($N=20$) → `0.1208±0.0963` ($N=32$) → `0.1214±0.0851` ($N=64$) → **`0.1177±0.0717`** ($N=128$)
  - `MSA+RAXML`: `0.1111±0.1093` ($N=20$) → `0.1203±0.0950` ($N=32$) → `0.1225±0.0847` ($N=64$) → `0.1191±0.0720` ($N=128$)
  - **`NJ (sss linear)`**: `0.1232±0.1025` ($N=20$) → `0.1325±0.0857` ($N=32$) → `0.1458±0.0659` ($N=64$) → **`0.1475±0.0487`** ($N=128$)
  - `MSA+NJ (wag_gamma)`: `0.1651±0.1208` ($N=20$) → `0.1800±0.1062` ($N=32$) → `0.1960±0.0968` ($N=64$) → `0.2155±0.0909` ($N=128$)
  - `NJ (poisson)`: `N/A` ($N=20$) → `0.1837±0.0891` ($N=32$) → `0.1977±0.0697` ($N=64$) → `0.2002±0.0503` ($N=128$)
  - **`NJ+NNI`（FastME BalME）**: $N=32, 64, 128$ の全条件（各 600 反復 × 3 手法 = 1,800 件）で元の `NJ` と **0/600 件の差（完全一致）**。

![Taxa数別 log(SSS) vs nRF](./results/results_nelsi_taxa128_L500/nelsi_taxa_scaling_nrf_vs_log_sss.png)

![Taxa数別 4手法の距離別 Error Bar & Box Plot](./results/results_nelsi_taxa128_L500/nelsi_taxa_scaling_by_distance_4methods_err_box.png)

- **$N=20, L=500$ における NELSI vs BioPerl 詳細比較・勝敗分析**:
  - $D=0.1 \sim 0.5$ の近〜中分岐では BioPerl の方が低 $nRF$ となる一方、$D=1.0 \sim 2.0$ では BioPerl が Twilight Zone（平均 SSS: $D=1.0$ で `0.0405`, $D=2.0$ で `0.0205`）へ落ち込み全手法の誤差が `0.43〜0.50` まで跳ね上がります。対照的に NELSI は $D=2.0$ でも平均 SSS `0.0538` を保ち、`NJ (sss linear)`（`0.216±0.123`）は `MSA+ML`（`0.231±0.143`）を上回ります。
  - $N=20$（$L \in [500, 1000, 1500]$ 全 1,800 ペア）での直接勝敗では、`NJ (sss linear)` は NELSI で `461/1800 (25.6%)` 勝（引き分け `37.5%`）、BioPerl で `344/1800 (19.1%)` 勝（引き分け `58.5%`）を記録し、高分岐 $D=1.5 \sim 2.0$ の重心では対角線 $y=x$ の下側（NJ 勝利側）へ入ります。

![Taxa=20 L=500 NELSI vs BioPerl 距離別比較](./results/results_nelsi_yule/nelsi_vs_bioperl_taxa20_L500_by_distance_err_box.png)

![Taxa=20 L=500 log(SSS) vs nRF 比較](./results/results_nelsi_yule/nelsi_vs_bioperl_taxa20_L500_nrf_vs_log_sss.png)

![Taxa=20 NJ (sss linear) vs MSA+ML 勝敗散布図](./results/results_nelsi_yule/nelsi_vs_bioperl_taxa20_nj_vs_ml_scatter.png)

- **配列数増加（$N=20 \to 32 \to 64 \to 128$）による `MSA+ML` 劣化の対比（BioPerl vs NELSI）**:
  - 同一条件（$D \in \{0.5, 1.0, 2.0\}, L=500, \text{indel}=0.1$）で [`JoiningBasedOnAdditivity`](../JoiningBasedOnAdditivity/results_hpc_indel0.1/main/hpc_benchmark_summary.csv) の BioPerl 結果と比較すると、**BioPerl では $N=20 \to 128$ で `MSA+ML` が `0.204 → 0.393`（`+0.190` 悪化）と急激に崩壊**する一方、`NJ (sss linear)` の悪化は `0.226 → 0.303`（`+0.077`）に留まり、**$N \ge 64$ で `NJ (sss linear)` が `MSA+ML` を逆転**します。
  - 一方、**NELSI では `MSA+ML`（`0.131 → 0.124`, `-0.007`）・`NJ (sss linear)`（`0.138 → 0.150`, `+0.012`）ともに配列数が増えてもほとんど悪化しません**。その物理的理由は平均 SSS の推移にあり、BioPerl では $N$ 増大に伴い平均 SSS が `0.047 → 0.008` へとトワイライトゾーン深部に落ち込むため MAFFT が崩壊するのに対し、NELSI（`root_to_tip` スケーリング）では $N=128$ でも平均 SSS `0.108` が保たれアライメントが崩壊しないためです。

![BioPerl vs NELSI Taxa数増加による悪化比較](./results/results_nelsi_taxa128_L500/nelsi_vs_bioperl_taxa_scaling_degradation.png)

---

### カテゴリ B: 配列の生物学的異質性・保存構造（AliSim 条件群）

- **`results_alpha` / `results_gamma`** (サイト間速度変動 $\alpha \in [0.25..2.0]$):
  - [alpha_scaling_curves.png](./results/results_alpha/alpha_scaling_curves.png) / [gamma_distance_alpha_scaling.png](./results/results_gamma/gamma_distance_alpha_scaling.png)
  - 速度不均一が強い（$\alpha = 0.25$）ほど高速置換サイトと保存サイトの二極化が進み、$+G4$ モデルを扱う `MSA+ML` の優位性が顕著になります。
- **`results_ics` / `results_ics_full`** (保存ドメイン比率 $ICS \in [0.0..1.0]$):
  - [ics_scaling_benchmark.png](./results/results_ics/ics_scaling_benchmark.png) / [ics_full_benchmark.png](./results/results_ics_full/ics_full_benchmark.png)
  - 部分保存領域（ICS）が存在すると MMseqs2 局所整列のアンカーが維持され、極限距離でも SSS のゼロ落ち込みが緩和されます。
- **`results_simple`** (均一速度 `LG` モデル):
  - [scaling_curves_high_dist.png](./results/results_simple/scaling_curves_high_dist.png)

---

### カテゴリ C: スケーラビリティ・極限距離・距離法アルゴリズム比較

- **`results_taxon`** ($N \in [8, 16, 32, 64, 128]$): [taxon_scaling_curves.png](./results/results_taxon/taxon_scaling_curves.png)
  - $N$ 増大に伴い探索空間が指数的に広がり、低 SSS 領域での精度劣化が加速します。
- **`results_high_dist`** ($D \in [4.0, 5.0, 6.0]$): [scaling_curves_high_dist.png](./results/results_high_dist/scaling_curves_high_dist.png)
  - SSS 最大値 `0.0910` 以下の極限領域では相同シグナルがほぼ失われ、全手法が $nRF > 0.70$ に収束します。
- **`results_high_gap` / `results_new` / `results_true_pwa`**:
  - [results_high_gap/regime_map_delta_nrf.png](./results/results_high_gap/regime_map_delta_nrf.png) / [results_new/regime_map_delta_nrf.png](./results/results_new/regime_map_delta_nrf.png) / [results_true_pwa/true_pwa_scaling_curves.png](./results/results_true_pwa/true_pwa_scaling_curves.png)
  - `TRUE_DIST+NJ`（理論下限）との対比により、低 SSS 域の誤差増大が主にアラインメント崩壊に起因することを分離確認しました。
- **`results_fastme_NoOption` / `results_fastme_options`** (FastME vs RapidNJ):
  - [scaling_curves_fastme.png](./results/results_fastme_NoOption/scaling_curves_fastme.png) / [scaling_curves_fastme_options.png](./results/results_fastme_options/scaling_curves_fastme_options.png)
  - SPR 局所探索＋LG+G 補正距離を用いた FastME は、単純な NJ より低類似度域で良好な探索性能を示します。

---

## 3. 総合的な結論とシミュレーターの変化に関する考察

本プロジェクトでは、配列進化シミュレーター・インデル分布・距離推定器・系統樹シミュレーターを段階的に発展させてきました。各段階での「シミュレーターの変化」が手法間の優劣と推定メカニズムに与えた影響を、以下の 4 点に総括します。

### 1. Indel モデルの変化（幾何分布 → Zipfian べき乗則）が引き起こす「アライメント相転移」
- **幾何分布 Indel（`results_alpha`, `results_new` 等）**: インデル長が短く配列全体に均等に生じるため、MAFFT は配列全長にわたって安定して相同アンカーを配置できます。この条件下では、置換モデル（`LG+G4`）を統計的に活用できる **`MSA+ML` が全距離帯で一貫して最強** となります。
- **Zipfian べき乗則 Indel（`POW 1.7 50`, `results_zipfian_indel` 以降）**: 最大 50 残基に及ぶ重尾（Heavy-tailed）の挿入・欠失が入ると、$SSS < 0.06$ の遠縁領域で多重整列がギャップの強制整列（過剰整列バイアス）を起こします。その結果、**アライメント誤差が尤度推定の前提を破壊し、`MSA+ML` の精度が急落して局所整列・グラフ分割に基づく `PSA+NJ` や `GS` に逆転される Crossover 現象** が普遍的に出現します。

### 2. 配列シミュレーターと距離計量の変化（AliSim → INDELible v1.03、NW 大域整列 → MMseqs2 局所線形距離）
- **AliSim から INDELible v1.03 への移行（`results_paper_tree2` 以降）**: 連続時間マルコフジャンプ過程で置換とインデルを厳密に共進化させる INDELible v1.03（WAG+G4, `POW 1.7 50`）の採用により、Matsui & Iwasaki (2020) の物理条件が忠実に再現され、極限領域（$SSS \le 0.03$）における `GS` の対 MSA 耐性（`0.3768`）が実証されました。
- **NW 大域整列から MMseqs2 局所線形距離への移行（`results_paper_tree3` → `results_paper_tree4`）**:
  - 大域 Needleman-Wunsch（`PT3`）は長大インデルの累積ギャップペナルティにより高分岐域で距離が飽和しました。
  - 一方、MMseqs2 局所整列相対類似度 $w_{ij}$ からの **有界線形非類似度 $d = 1 - w_{ij}$（`PT4`）** は、対数変換 $d = -\ln w_{ij}$ が $w_{ij} \to 0$ 近傍で引き起こす距離発散と三角形不等式の破綻を回避し、保存ドメインのシグナルのみを安定して抽出するため、$D = 1.5 \sim 2.0$ の遠縁領域で `MSA+NJ`（WAG+G4）および `MSA+ML` を上回る精度を達成しました。

### 3. 系統樹シミュレーターの変化（Birth-Death → BioPerl 後退 Yule → NELSI 自己相関緩和時計）
同じ配列シミュレーター（INDELible v1.03）と同一のパラメータグリッド（$N=20, D \in [0.1..2.0], L \in [500, 1000, 1500]$）の下でも、**系統樹シミュレーターの変更は「配列類似度の分布」と「手法ごとのボトルネック」を根本から変化**させました。

1. **BioPerl 後退 Yule（`PT3`, `PT4`）— 「アライメント崩壊」支配型**:
   - 後退 Yule 過程は各枝長を $1 - \ln(u(e-1)+1)$ から抽出し直接 $D$ を掛けるため、同一 $D$ でもルート〜葉の累積進化距離が大きくなり、SSS は平均 `0.1554`（中央値 `0.0573`）と急速にトワイライトゾーンへ突入します。
   - そのため、正解アライメント下では NJ と ML の差が小さい一方（`TRUE_MSA+NJ`: `0.0775 / 0.0509 / 0.0435` vs `TRUE_MSA+ML`: `0.0905 / 0.0492 / 0.0272`）、推定アライメントでは **MAFFT の誤整列がボトルネック** となり、多重整列を回避する `GS` や MMseqs2 線形 `PSA+NJ` が強く優位に立ちます。

2. **NELSI 自己相関緩和時計（Yule + Kishino モデル, $v=0.3$, `results_nelsi_yule`）— 「系統間速度差」支配型**:
   - 純粋出生 Yule 時間木に対し、祖先から子孫へ進化速度が対数正規に自己相関しながら変化する緩和時計を適用し、平均 root-to-tip 距離を $D$ にスケーリングしました。
   - 第一に、同一 $D$ あたりの平均配列類似度が **SSS 平均 `0.3178`（中央値 `0.2137`）** と BioPerl 時の約 2 倍以上に高く保たれるため、アライメント崩壊の割合が減り、`MSA+ML`（`0.1097 / 0.0839 / 0.0697`）および `PSA+NJ`（`0.1232 / 0.0908 / 0.0775`）の全体誤差が大幅に低下しました。
   - 第二に、注目すべき逆転現象として、正解アライメント同士の比較において **`TRUE_MSA+NJ` の誤差が BioPerl 時（`0.0775 / 0.0509 / 0.0435`）から NELSI 時（`0.1342 / 0.0944 / 0.0735`）へと約 1.7〜1.9 倍悪化** し、`TRUE_MSA+ML`（`0.0725 / 0.0501 / 0.0354`）との間に明確な乖離が生じました。これは、自己相関緩和時計によって生じた **系統間の進化速度の不均一性（特定のクレードが速く・別のクレードが遅く進化する偏り）** により、ペアワイズ距離行列の推定分散や長枝・短枝混在の影響が NJ 法のトポロジー復元を直接劣化させる一方、全サイト・全枝の尤度を同時最適化する最尤法（ML）が速度不均一に対して高い頑健性を保つことを示しています。

### 4. 配列数スケーリング（$N=20 \to 128$）と「モデル依存のアライメント崩壊」の本質
今回の Taxa スケーリング比較（Section A-6）により、**「配列数が増えると ML が悪化するのか、それともシミュレーターの木スケール定義によって SSS が崩壊しているのか」** が明確に切り分けられました。
- **BioPerl で `MSA+ML` が配列数に対して脆弱な理由**: BioPerl（後退 Yule）および従来のスケーリングでは $N$ の増加とともに木全体の総距離が膨らみ、平均 SSS が `0.047 (N=20) → 0.008 (N=128)` まで低下します。このため $N \ge 64$ では MAFFT の多重整列がほぼ全壊し、最尤法がその誤整列を尤度に組み込むことで **`+0.190` もの精度劣化**を起こします。アライメントフリーな **`NJ (sss linear)` は MSA 崩壊の影響を受けないため、BioPerl 大規模木（$N=64, 128$）では ML を逆転**しました。
- **NELSI で `MSA+ML` と `NJ (sss linear)` がともに頑健な理由**: NELSI（`root_to_tip = D`）では $N=128$ でも平均 SSS `0.108` が維持され、さらに種数が増えることで長枝を分断するタクソンサンプリング効果（LBA 緩和）が働くため、`MSA+ML`（`-0.007`）・`NJ (sss linear)`（`+0.012`）ともに悪化しません。また NELSI 環境下の全距離法の中では、`MSA+NJ (wag_gamma)` や `NJ (poisson)` よりも **`NJ (sss linear)` が一貫して最良（$N=128$ で `0.1475` vs `0.2155 / 0.2002`）** であり、最深部 $D=2.0$ では `MSA+ML` をも上回りました。

### 5. 短い進化距離（高 SSS 領域, $D = 0.1 \sim 0.2$）で NELSI が BioPerl より悪化する理由
$D = 1.0 \sim 2.0$ の遠縁領域では NELSI の方が高精度であるにもかかわらず、**短い進化距離（$D = 0.1 \sim 0.2$、※配列類似度 $SSS$ としては大きい $SSS \approx 0.60 \sim 0.77$ の高類似度領域）では逆に NELSI の $nRF$ が BioPerl より悪く（高く）なります**（$N=20, L=500, D=0.1$ での `MSA+ML`: NELSI `0.0659` vs BioPerl `0.0012`、`NJ (sss linear)`: `0.0988` vs `0.0012`、`GS`: `0.2406` vs `0.1071`）。

実測の枝長分布および正解アライメント対照（`TRUE_MSA+*`）から、その原因は以下の **3 点（情報不足・ソフトポリトミー支配）** として説明できます：

1. **同一 $D$ あたりの実質進化距離が約 1/3 であり、内部枝に十分な置換が入らない（確率的情報限界）**:
   - $D = 0.1, N = 20$（20 反復実測）の木構造を比較すると、**BioPerl (`paper_yule`)** は正規化を行わない後退 Yule 過程であるため、名目 $D=0.1$ でも **root-to-tip 平均距離が `0.2927`（総枝長 `1.5868`）、内部枝長の中央値が `0.0587`（下位10%点 `0.0231`）** と十分に長く、$L = 500$ において期待置換数が 1 回未満となる内部枝（$b \cdot 500 < 1$）は **`0.0%`** です。
   - これに対し **NELSI (`nelsi_kishino`)** は平均 root-to-tip 距離を厳密に `0.1000` にスケーリングするため、**総枝長平均が `0.7181`、内部枝長の中央値が `0.0137`（下位10%点 `0.0029`）** と BioPerl の約 1/4 に縮みます。その結果、$L = 500$ では **全内部枝の `7.8%` で期待置換数が 1 回未満（`2.8%` で 0.5 回未満）** となり、分岐後にアミノ酸置換が 1 度も起きず配列情報から原理的に復元不能な「ソフト・ポリトミー（zero-information branch）」が構造的に生じます。

2. **アライメント誤差ではなく、`TRUE_MSA` の段階ですでに誤差が生じている**:
   - $D = 0.1, L = 500$ における NELSI の **`TRUE_MSA+ML` の平均 $nRF$ は `0.0635`** であり、推定アライメントを用いた `MSA+ML`（`0.0659`）との差はわずか `0.0024` に過ぎません（BioPerl では `TRUE_MSA+ML = 0.0012`）。
   - さらに NELSI では $D = 0.1 \to 0.2$ へ進化距離を**少し伸ばした方が全手法で $nRF$ が明確に低下**します（`MSA+ML`: `0.0659 → 0.0565`、`NJ (sss linear)`: `0.0988 → 0.0729`、`GS`: `0.2406 → 0.2147`）。これは **$D = 0.1$ がアライメント崩壊（遠縁エラー）ではなく「置換イベント不足による系統シグナル欠乏領域（Under-divergence regime）」** であることを直接裏付けています。

3. **Yule 時間木＋緩和時計による極短内部枝の発生と、距離法・`GS` での識別限界**:
   - 純粋出生 Yule 時間木では先端付近や連続種分化ノードに短い内部枝が集中しやすく、自己相関緩和時計によってそこに低置換速度が重なると、同程度の $SSS$（たとえば $SSS \approx 0.6 \sim 0.8$）であっても局所的に未解決分岐が残ります。
   - このとき 1〜数残基の違いしか生じない微細分岐に対して、ペアワイズ距離行列や局所類似度重みを使う `TRUE_MSA+NJ`（`0.0976`）・`GS`（`0.2406`）は微小な距離雑音・MMseqs2 スコアの離散化の影響を受けやすく、サイト単位の尤度を用いる `TRUE_MSA+ML`（`0.0635`）よりも精度低下が顕著になります（`log(SSS) vs nRF` 曲線で見られる $SSS \to 0.9$ 端の右上がり「J字フック」の主要因）。

### 6. 研究全体としての設計指針
- **極端な長大インデルかつ超高分岐（$SSS \le 0.06$）が支配的なデータ**（BioPerl 後退 Yule 相当や大規模高発散群）では、多重整列の崩壊を防ぐ **MMseqs2 局所類似度＋有界線形距離 (`NJ (sss linear)`)** または **Graph Splitting (`GS`)** が最も頑健です。
- **現実的な時間木トポロジー＋系統間速度変動（自己相関緩和時計）を伴う中〜高分岐データ**（NELSI Yule 相当）では、アライメントが破綻しない限り **最尤法（`MSA+ML` / `MSA+RAXML`）** が系統間の速度差を吸収して最良の精度を発揮します。ただし超高分岐端（$D=2.0$ 付近）や大規模データにおける高速代替法としては、**`NJ (sss linear)` がポアソン距離 NJ・WAG+G4 MSA+NJ を大きく引き離す最良の距離法** となります。

---

## 4. 今後確かめるべき課題・ネクストステップ

これまでの比較解析から、以下の 5 点を次に検証することで研究全体の主張をより強固にできます。

1. **無相関緩和時計モデル（`nelsi_uncor_lnorm`）や速度変動幅 $v$ を変えた検証**:
   - 自己相関モデル（`nelsi_kishino`）では近縁クレード間で速度が似るため局所的な速度差が穏やかでした。すでに作成済みの [`nextflow_nelsi_uncor_lnorm.config`](./next_configs/nextflow_nelsi_uncor_lnorm.config) を実行し、枝ごとに独立に速度が跳ぶ無相関対数正規モデル（UCLN）で **`NJ (sss linear)` と `MSA+ML` の優劣・系統間速度差耐性がどう変化するか** を確かめる必要があります。
2. **NELSI 環境下での長枝牽引（LBA: `lba_ratio = 3.0, 5.0`）の影響検証**:
   - Yule＋自己相関による速度変動に意図的な LBA を重ねたとき、距離法（`NJ (sss linear)`, `NJ (poisson)`）と最尤法（`MSA+ML`）のどちらが先に LBA に引きずられるかを、準備済みの [`run_nextflow_nelsi_lba.sh`](./run_scripts/run_nextflow_nelsi_lba.sh) および [`analyze_nelsi_lba.py`](./analysis/analyze_nelsi_lba.py) で検証する必要があります。
3. **NELSI における極限低 SSS 領域（$D = 3.0 \sim 6.0$ や $L < 500$）での Crossover 確認**:
   - 今回 $D=2.0$ 時点で `NJ (sss linear)` が平均 $nRF$ で `MSA+ML` を上回りましたが（$N=20, L=500$: `0.216` vs `0.231`）、NELSI では $D=2.0$ でも平均 SSS が `0.0538` に留まっており、BioPerl のような $SSS \le 0.03$ の完全トワイライトゾーンにはほとんど入っていません。$D=3.0 \sim 5.0$ まで伸ばした際、NELSI でもアライメント崩壊による明確な逆転が起きるかを確認すべきです。
4. **`N=20` での NELSI + Poisson 距離の補完と超大規模 Taxa（$N=256, 512$）の検証**:
   - $N=20$ 条件で未実施となっている NELSI の `NJ (poisson)` データを補完するとともに、$N=256 \sim 512$ でも `MSA+ML` の優位性と `NJ (sss linear)` のスケーラビリティ・実行時間優位がどこまで保たれるかを確かめる価値があります。
5. **BalME NNI が NJ と一致した理由の掘り下げと、SPR 探索・距離加法補正の効果検証**:
   - $N=32..128$ の全 1,800 反復で FastME BalME NNI が NJ と `0/600` 件の差（完全一致）となったことから、単純な近傍交換（NNI）だけでは距離行列由来の誤トポロジーを脱出できないことが示唆されます。**SPR 探索（`fastme -s`）** や **`NJ (sss linear)` に対する加法補正・非線形距離変換** が、$D=2.0$ 付近の精度をさらに引き上げられるかを検証する必要があります。
