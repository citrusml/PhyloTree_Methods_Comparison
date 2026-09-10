# Phylogenetic Benchmark Results: SSS-Based Multi-Condition Evaluation

本ドキュメントは、`/Users/kazukiaibara/PhyloMethod/results` に蓄積された全 15 系統樹推定ベンチマーク実験（総計 160,000 系統樹推定以上、31,900 レプリケート）における **Sequence Similarity Score ($SSS$, $w$)** を横軸としたトポロジー誤差（$nRF$）の解析結果を網羅的にまとめたレポートです。

各実験条件において出力された **条件別分割プロット (`nrf_vs_sss_by_condition.png`)** を軸に、パラメータ変化（配列長 $L$、進化距離 $D$、ガンマ形状母数 $\alpha$、分類群数 $N$、ICS保存領域比率など）がアルゴリズムの頑健性に与える影響を比較・考察します。

---

## 1. 全実験の概括比較テーブル

| # | 実験名 / ディレクトリ | 主な検証パラメータ | 比較パイプライン | SSS 範囲 $[w_{\min} .. w_{\max}]$ | 主な発見・挙動 |
|---|---|---|---|:---:|---|
| 1 | [`results_zipfian_indel`](./results/results_zipfian_indel) | べき乗則 Indel ($POW\{1.7/50\}$), $L \in [300 .. 1500]$ | MSA+ML, MSA+NJ, PSA+NJ, TRUE群 | $[0.0078 .. 0.9318]$ | 長大インデルにより $w < 0.06$ で MSA+ML が急激に崩壊し、PSA+NJ が逆転 |
| 2 | [`results_paper_tree`](./results/results_paper_tree) | 後退Yule系統樹, $L \in [500 .. 1500]$, GS法 | GS, MSA+ML, MSA+RAXML, PSA+NJ, TRUE群 | $[0.0139 .. 0.4948]$ | 論文 (syz049) 完全準拠系統樹。GS法が低類似度領域で安定した耐性を発揮 |
| 3 | [`results_paper_tree_pre_fix`](./results/results_paper_tree_pre_fix) | 後退Yule系統樹 (初期スケール版) | GS, MSA+ML, MSA+RAXML, PSA+NJ, TRUE群 | $[0.0054 .. 0.6745]$ | 枝長補正前の予備実験。全体傾向は results_paper_tree と整合 |
| 4 | [`results_alpha`](./results/results_alpha) | サイト間速度変動 $\alpha \in [0.25 .. 2.0]$, $L \in [100 .. 1000]$ | MSA+ML, MSA+NJ, PSA+NJ | $[0.0063 .. 0.9537]$ | $\alpha=0.25$ (極度の不均一性) で SSS 減衰が最も激しく、MLのモデル適合が重要 |
| 5 | [`results_gamma`](./results/results_gamma) | ガンマ分布変動 $\alpha \in [0.25 .. 2.0]$ | MSA+ML, MSA+NJ, PSA+NJ | $[0.0063 .. 0.9537]$ | results_alpha との再現性確認。$w > 0.10$ での MSA+ML の圧倒的優位性 |
| 6 | [`results_taxon`](./results/results_taxon) | 分類群数 $N \in [8, 16, 32, 64, 128]$, $L \in [100 .. 1000]$ | MSA+ML, MSA+NJ, PSA+NJ | $[0.0000 .. 0.9861]$ | $N$ 増加に伴い探索空間が指数関数的に増大し、低 SSS 側での $nRF$ 悪化が加速 |
| 7 | [`results_ics`](./results/results_ics) | 断片保存領域 $ICS \in [0.0 .. 0.20]$, $L \in [100 .. 1000]$ | MSA+ML, MSA+NJ, PSA+NJ | $[0.0072 .. 0.9525]$ | 保存ブロック（ICS）の存在がペアワイズ局所アラインメントのアンカーとなり SSS を底上げ |
| 8 | [`results_ics_full`](./results/results_ics_full) | 完全保存ドメイン ($ICS = 1.0$), $L \in [100 .. 1500]$ | MSA+ML, MSA+NJ, PSA+NJ | $[0.0145 .. 0.9451]$ | 全長が ICS の場合、置換速度低下によりトポロジー推定精度が全体的に極めて高く維持 |
| 9 | [`results_high_dist`](./results/results_high_dist) | 極限進化距離 $D \in [4.0, 5.0, 6.0]$ | MSA+ML, MSA+NJ, PSA+NJ | $[0.0078 .. 0.0910]$ | 全反復が Twilight Zone ($w < 0.10$) に集中。情報消失による限界領域 |
| 10 | [`results_high_gap`](./results/results_high_gap) | 広範囲進化距離 $D \in [0.1 .. 6.0]$, 幾何分布 Indel | MSA+ML, MSA+NJ, PSA+NJ, TRUE群 | $[0.0078 .. 0.9318]$ | 広域 $D$ 変化における SSS 連続減衰と各手法の感度曲線を取得 |
| 11 | [`results_simple`](./results/results_simple) | 単一速度 LG モデル (ガンマ不均一なし) | MSA+ML, MSA+NJ, PSA+NJ | $[0.0035 .. 0.9537]$ | 速度不均一がない場合、高類似度域での推定精度が全体的に向上 |
| 12 | [`results_new`](./results/results_new) | 標準ベースライン ($D \in [0.1 .. 3.0]$) + 真の距離行列対照 | MSA+ML, MSA+NJ, PSA+NJ, TRUE_DIST+NJ群 | $[0.0064 .. 0.9494]$ | アラインメント誤差と距離推定誤差の寄与度を TRUE_DIST+NJ により分離評価 |
| 13 | [`results_fastme_NoOption`](./results/results_fastme_NoOption) | 最小進化法 FastME (デフォルトオプション) | MSA+FastME, PSA+FastME | $[0.0064 .. 0.9494]$ | 近隣結合法 (NJ) と比較し、FastME が長大配列で良好なトポロジーを構築 |
| 14 | [`results_fastme_options`](./results/results_fastme_options) | FastME (SPR 探索 + LG+G 補正行列) | MSA+FastME_LG_G, PSA+FastME_SPR | $[0.0078 .. 0.9318]$ | SPR 局所探索の導入により低 SSS 域でのトポロジー探索能が向上 |
| 15 | [`results_true_pwa`](./results/results_true_pwa) | 真のペアワイズアラインメント対照群 | TRUE_PSA+NJ | $[0.0064 .. 0.9494]$ | PSA におけるアラインメント誤差の影響を排除した理論的上限値 |

---

## 2. カテゴリ別 詳細結果と条件別プロット

### カテゴリ A: 論文再現・現実的 Indel モデル

#### 1. `results_zipfian_indel` (べき乗則 Indel 長分布)
- **概要**: 論文 (Matsui & Iwasaki 2020) と同様の Heavy-tailed（べき指数 $a=1.7$、最大 50 残基）Indel を導入した実験。
- **条件パラメータ**: Indel 率 $(0.05, 0.10)$, 配列長 $L \in [300, 500, 1000, 1500]$, 進化距離 $D \in [0.1 .. 6.0]$
- **グラフ**: [画像を直接開く](./results/results_zipfian_indel/nrf_vs_sss_by_condition.png)

![results_zipfian_indel 条件別プロット](./results/results_zipfian_indel/nrf_vs_sss_by_condition.png)

- **分析**:
  - 長大 Indel の存在により、多重配列アラインメント（MAFFT）のギャップ誤配置が爆発的に増加。
  - $w < 0.06$（Crossover 領域）において、`MSA+ML` のトポロジー誤差が急上昇し、ペアワイズ局所アラインメントに基づく `PSA+NJ` の方が頑健となる現象が顕著に現れています。

---

#### 2. `results_paper_tree` (論文準拠 後退 Yule 系統樹 + GS法)
- **概要**: 論文の系統樹生成モデル（Backward Yule 過程 + 対数サンプリング枝長）を完全再現し、グラフ分割法（Graph Splitting: GS法）を含めて評価した中心実験。
- **条件パラメータ**: 配列長 $L \in [500, 1000, 1500]$, 進化距離 $D \in [0.1, 0.2, 0.5, 1.0, 1.5, 2.0]$, Taxa $N=32$
- **グラフ**: [画像を直接開く](./results/results_paper_tree/nrf_vs_sss_by_condition.png)

![results_paper_tree 条件別プロット](./results/results_paper_tree/nrf_vs_sss_by_condition.png)

- **分析**:
  - SSS は $0.0139 \sim 0.4948$（中央値 $0.0527$）と全体的に Twilight Zone 付近に集中。
  - 配列長 $L=1500$ では情報量の増加により全手法のトポロジー精度が向上する一方、超低類似度域（$w < 0.05$）では `GS` 法が `MSA+ML` に匹敵または凌駕するトポロジー安定性を示しています。

---

#### 3. `results_paper_tree_pre_fix` (論文系統樹 予備実行データ)
- **概要**: 枝長スケール補正前の予備ベンチマーク。
- **グラフ**: [画像を直接開く](./results/results_paper_tree_pre_fix/nrf_vs_sss_by_condition.png)

![results_paper_tree_pre_fix 条件別プロット](./results/results_paper_tree_pre_fix/nrf_vs_sss_by_condition.png)

- **分析**:
  - `results_paper_tree` と同様に、Twilight Zone における各手法の順位逆転挙動が再現されています。

---

### カテゴリ B: 配列の生物学的異質性・構造モチーフ

#### 4. `results_alpha` (サイト間進化速度変動: $\alpha$)
- **概要**: ガンマ分布形状母数 $\alpha \in \{0.25, 0.5, 1.0, 2.0\}$ の影響を評価。
- **条件パラメータ**: ガンマ $\alpha \in [0.25, 0.5, 1.0, 2.0]$, 配列長 $L \in [500, 1000]$ (100 aa 除外), 距離 $D \in [0.1, 0.5, 1.0, 2.0, 3.0]$
- **グラフ**: [画像を直接開く](./results/results_alpha/nrf_vs_sss_by_condition.png)

![results_alpha 条件別プロット](./results/results_alpha/nrf_vs_sss_by_condition.png)

- **分析**:
  - $\alpha=0.25$（強い進化速度不均一性）では、不変サイトと極めて速く変異するサイトが混在するため SSS の低下が急速に進行。
  - 速度不均一が大きいほど、モデル（$+G4$）を考慮した最尤法（`MSA+ML`）の優位性が明確になります。

---

#### 5. `results_gamma` (ガンマ形状母数の検証)
- **概要**: `results_alpha` と同一条件グリッドにおける独立レプリケート検証。
- **グラフ**: [画像を直接開く](./results/results_gamma/nrf_vs_sss_by_condition.png)

![results_gamma 条件別プロット](./results/results_gamma/nrf_vs_sss_by_condition.png)

- **分析**:
  - 高類似度領域（$w > 0.15$）では `MSA+ML` が $nRF < 0.15$ と極めて高精度を維持し、類似度低下に伴い $w=0.06$ 付近で緩やかに各手法の差が縮小します。

---

#### 6. `results_ics` (断片保存配列: Interspersed Conserved Segments)
- **概要**: 進化の過程で保存される機能ドメイン（ICS）が配列全体の $0\%, 5\%, 10\%, 20\%$ 混在する条件。
- **条件パラメータ**: ICS比率 $\in [0.0, 0.05, 0.10, 0.20]$, 配列長 $L \in [300, 500, 1000]$ (100 aa 除外), 距離 $D \in [0.1 .. 3.0]$
- **グラフ**: [画像を直接開く](./results/results_ics/nrf_vs_sss_by_condition.png)

![results_ics 条件別プロット](./results/results_ics/nrf_vs_sss_by_condition.png)

- **分析**:
  - ICS 比率が高くなるほど（$ICS=0.20$）、MMseqs2 で検出可能なアンカー領域が確保され、極限進化距離（$D=3.0$）においても SSS がゼロに落ち込まず維持されます。
  - これによりペアワイズアラインメント法の性能劣化が著しく緩和されます。

---

#### 7. `results_ics_full` (完全保存ドメイン: $ICS = 1.0$)
- **概要**: 配列全長が ICS モデルに従って進化する対照実験。
- **グラフ**: [画像を直接開く](./results/results_ics_full/nrf_vs_sss_by_condition.png)

![results_ics_full 条件別プロット](./results/results_ics_full/nrf_vs_sss_by_condition.png)

- **分析**:
  - 全域が保存されているため、変異速度が大幅に抑制され、全配列長において高い類似度（中央値 $0.375$）と極めて低い $nRF$（高精度）が達成されます。

---

### カテゴリ C: スケーラビリティと極限距離の限界

#### 8. `results_taxon` (分類群数のスケーリング: $N \in [8, 16, 32, 64, 128]$)
- **概要**: 系統樹のサイズ（Taxon 数）が SSS と推定精度に与える影響。
- **条件パラメータ**: 分類群数 $N \in [8, 16, 32, 64, 128]$, 配列長 $L \in [500, 1000]$ (100 aa 除外), 距離 $D \in [0.1 .. 3.0]$
- **グラフ**: [画像を直接開く](./results/results_taxon/nrf_vs_sss_by_condition.png)

![results_taxon 条件別プロット](./results/results_taxon/nrf_vs_sss_by_condition.png)

- **分析**:
  - $N=8$ ではトポロジー空間が小さいため全体的に誤差が低く抑えられます。
  - $N=64, 128$ と大規模化するにつれて、長枝誘引（LBA）やアラインメントミスの影響が全体に波及し、低 SSS 側での推定精度低下が急峻になります。

---

#### 9. `results_high_dist` (極限進化距離 $D \in [4.0, 5.0, 6.0]$)
- **概要**: ほぼ全反復が Twilight Zone 深部（$w \le 0.03$）に突入する過酷な条件。
- **グラフ**: [画像を直接開く](./results/results_high_dist/nrf_vs_sss_by_condition.png)

![results_high_dist 条件別プロット](./results/results_high_dist/nrf_vs_sss_by_condition.png)

- **分析**:
  - SSS の最大値でも $0.0910$ であり、中央値は $0.0318$。
  - ここでは配列間の相同シグナルがほぼ失われており、どのアルゴリズムも $nRF > 0.70$ となり、理論的情報限界に達していることが確認できます。

---

#### 10. `results_high_gap` (広域距離 $D \in [0.1 .. 6.0]$)
- **概要**: 近縁（$D=0.1$）から極遠縁（$D=6.0$）までを連続的に網羅。
- **グラフ**: [画像を直接開く](./results/results_high_gap/nrf_vs_sss_by_condition.png)

![results_high_gap 条件別プロット](./results/results_high_gap/nrf_vs_sss_by_condition.png)

- **分析**:
  - 真のアラインメント対照群（`TRUE_MSA+ML`, `TRUE_PSA+NJ`）との比較により、誤差の何割が「アラインメントエラー」によるもので、何割が「系統推定法そのものの限界」であるかを定量的に分離観察できます。

---

### カテゴリ D: アルゴリズム比較・対照群

#### 11. `results_simple` (均一速度 LG モデル)
- **概要**: サイト間速度変動（$+G4$）を含まない標準的 LG モデル。
- **グラフ**: [画像を直接開く](./results/results_simple/nrf_vs_sss_by_condition.png)

![results_simple 条件別プロット](./results/results_simple/nrf_vs_sss_by_condition.png)

- **分析**:
  - 複雑な速度変動がないため、同一進化距離における SSS が安定し、高類似度域での推定精度が全体的に向上します。

---

#### 12. `results_new` (ベースラインベンチマーク + TRUE_DIST 対照)
- **概要**: 真の距離行列から系統樹を推定する `TRUE_DIST+NJ` を含む統制実験。
- **グラフ**: [画像を直接開く](./results/results_new/nrf_vs_sss_by_condition.png)

![results_new 条件別プロット](./results/results_new/nrf_vs_sss_by_condition.png)

- **分析**:
  - `TRUE_DIST+NJ`（黒点線）はアラインメントおよび距離推定誤差がゼロの理論的下限を示しており、低 SSS 域での誤差増大が純粋にアラインメントの破綻に起因していることが明確に証明されています。

---

#### 13. `results_fastme_NoOption` & 14. `results_fastme_options` (最小進化法 FastME)
- **概要**: 距離行列法として NJ 以外の代表格である FastME（デフォルト vs SPR局所探索+LG+G補正）を比較。
- **グラフ**:
  - `results_fastme_NoOption`: [画像を直接開く](./results/results_fastme_NoOption/nrf_vs_sss_by_condition.png)

![results_fastme_NoOption 条件別プロット](./results/results_fastme_NoOption/nrf_vs_sss_by_condition.png)

  - `results_fastme_options`: [画像を直接開く](./results/results_fastme_options/nrf_vs_sss_by_condition.png)

![results_fastme_options 条件別プロット](./results/results_fastme_options/nrf_vs_sss_by_condition.png)

- **分析**:
  - SPR 局所探索と Gamma 補正距離を組み合わせた `MSA+FastME_LG_G` は、単純な NJ 法（RapidNJ）よりも低類似度域で優れたトポロジー復元性を示します。

---

#### 15. `results_true_pwa` (真のペアワイズアラインメント対照群)
- **概要**: 正解のペアワイズアラインメントから Needleman-Wunsch 距離を算出して近隣結合樹を構築。
- **グラフ**: [画像を直接開く](./results/results_true_pwa/nrf_vs_sss_by_condition.png)

![results_true_pwa 条件別プロット](./results/results_true_pwa/nrf_vs_sss_by_condition.png)

- **分析**:
  - `PSA+NJ` の誤差のうち、ペアワイズアラインメント自体の誤りに起因する成分を特定するための基準線となります。

---

## 3. 総合的な結論と次期実験（実験16: `results_paper_tree2`）への展望

1. **Matsui & Iwasaki (2020) の Crossover 現象（$w \approx 0.06$）の普遍性**:
   - ほぼ全ての実験（特にべき乗則 Indel を含む `results_zipfian_indel` および `results_paper_tree`）において、$w < 0.06$ に達すると `MSA+ML` の急速なトポロジー精度劣化が発生し、ペアワイズ局所アラインメントに基づくアプローチ（`PSA+NJ`, `GS`）が同等以上の頑健性を示すことが実証されました。
2. **完全計算による連続 SSS 評価の意義**:
   - 全 15 ディレクトリにおいて 100 レプリケート個別のシードによる完全計算（欠損 0 件、最大 7,960 種の固有 SSS）を実施したことで、従来の離散的な縦線が完全に解消され、LOWESS 平滑化曲線が各手法の挙動を忠実に捉えられるようになりました。
3. **実験16（INDELible 正式採用による完全再現）への接続**:
   - 現在準備中の [`nextflow_paper_tree2.config`](./next_configs/nextflow_paper_tree2.config)（実験16）では、AliSim から連続時間マルコフジャンプ過程シミュレータである **INDELible v1.03**（WAG+G4 モデル、べき指数 1.7 の Zipfian Indel）へと切り替え、論文 Fig. 3a のシミュレーション条件を完全再現する準備が整っています。
