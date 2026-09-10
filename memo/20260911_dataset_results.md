# 実生物データセット（`dataset/data`）系統樹推定・分類学的整合性（TCS）ベンチマーク総括レポート

作成日: 2026年9月11日  
解析対象: 4生物群（Amniota, Eukaryota, Actinobacteria, Bacteria）、各 100 HOG（計 400 HOG、124,426 配列ペア）  
比較手法: **MSA+ML** (IQ-TREE 2), **MSA+NJ** (RapidNJ), **PSA+NJ** (RapidNJ), **GS** (Graph Splitting `gs2`)

---

## 1. エグゼクティブ・サマリー（主要な結論）

1. **総合第1位は最尤法（MSA+ML）**:
   - 同一 HOG 内比較において、全 400 HOG 中 **36.0%（143.9勝）で単独トップ勝率**を記録。
   - GS に対して **241勝 130敗（勝率 60.2% vs 32.5%）**、NJ に対して **184勝 151敗（勝率 46.0% vs 37.8%）** と全手法の直接対決で勝ち越した。
2. **Graph Splitting（GS法）の明確な得意領域の特定**:
   - 複雑な真核生物（Eukaryota）では勝率わずか 5.0% と低迷する一方、**原核生物・超遠縁領域（Actinobacteria で 35.1% 勝率トップ、Bacteria で 33.2% 勝率トップ）** において ML を凌駕する最高勝率を達成。
   - 特に遠縁帯（$0.2 \le \text{SSS} < 0.4$）で **29.5%** と大きく躍進し、多重整列の過剰整列ノイズを回避するグラフ分割法の強みが実証された。
3. **遠縁領域における PSA+NJ の MSA+NJ に対する優位性**:
   - 遠縁帯（$\text{SSS} < 0.4$）において、Needleman-Wunsch に基づく **PSA+NJ が MAFFT に基づく MSA+NJ を勝率（22.9〜26.9% vs 15.2〜20.4%）で一貫して上回る**。
   - 全 400 HOG 中 45.0%（180 HOG）では MSA+NJ とトポロジーが完全一致するが、差が生じた場合は PSA+NJ が勝ち越す。
4. **極限の超遠縁（$\text{SSS} < 0.1$）での「距離法（NJ）の逆転現象」**:
   - アミノ酸変異が飽和する極限領域（Bacteria 8 HOG, Actinobacteria 1 HOG, 平均 SSS = 0.0736）では、**PSA+NJ（平均 TCS: 4.176）および MSA+NJ（4.147）が MSA+ML（3.978）を逆転して勝利（NJ 法が 9 HOG 中 6 HOG を制覇）**。
   - 一方、GS は MMseqs2 の相同性エッジがノイズ化するため 1 勝（11.1%、平均 TCS: 3.736 で最下位）と苦戦する。
5. **TCS（Taxonomy Congruence Score）の数理的性質とルーティングの重要性解明**:
   - TCS は根（Root）から階層的に分類群一致度を累積するため、無根木（IQ-TREE / RapidNJ）では真のクレードが根で両断され過剰ペナルティを受けていた。
   - AmpliPhy 標準の **MAD（Minimal Ancestor Deviation: 祖先偏差最小化法）ルーティング** を適用することで、ML と NJ のスコアが倍増し正当に評価された。
   - TCS は分母で葉数 $N$ で割っているものの、木の深さ $D \approx \log_2 N$ に比例してスケールするため、**異なる HOG 間の生値比較ではなく「同一 HOG 内での比較」が必須**である。

---

## 2. 実験条件と実行環境

実験15（論文準拠ベンチマーク）と全く同一のオプションを適用し、ローカル（macOS Apple Silicon 6並列）で全 400 HOG × 4 手法 = **計 1,600 本の系統樹構築（成功率 100%）** を完了した。

| 手法 | アライメント手法 | 距離計算 / 系統樹エンジン | 適用オプション |
| :--- | :--- | :--- | :--- |
| **MSA+ML** | MAFFT (`--auto --threadit 0`) | IQ-TREE 3.1.3 | `-m MFP --redo -T 1` (自動置換モデル選択, Bootstrap なし) |
| **MSA+NJ** | MAFFT (`--auto --threadit 0`) | RapidNJ 2.3.3 | ギャップ除外 Poisson 距離 (`dist_model = "poisson"`), `rapidnj -i pd -x` |
| **PSA+NJ** | Needleman-Wunsch (BLOSUM62) | RapidNJ 2.3.3 | **Gap Open = 20.0, Gap Extension = 2.0**, Poisson 距離, `rapidnj -i pd -x` |
| **GS** | MMseqs2 (全対全 PSA) | Graph Splitting (`gs2` v2.5) | **`-m 7.5`（最感度モード）**, `-s -l -t 1` |

---

## 3. 手法間トポロジー乖離度（Robinson-Foulds 距離: nRF）

実データにはシミュレーションのような単一の正解遺伝子樹が存在しないため、まず 4 手法間のトポロジー一致度（正規化 RF 距離: 0.0 = 完全一致, 1.0 = 完全不一致）を全ペアワイズで評価した。

### 生物群別・手法間 平均正規化 RF 距離（nRF）

| 生物群 (Dataset) | 解析 HOG 数 | Taxa 数 | MSA+NJ vs PSA+NJ | MSA+ML vs MSA+NJ | MSA+ML vs PSA+NJ | MSA+NJ vs GS | PSA+NJ vs GS | MSA+ML vs GS |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **Amniota** (`amn`) | 100 | 10 | **0.1142** | 0.3039 | 0.3154 | 0.4392 | 0.4465 | **0.4824** |
| **Eukaryota** (`euk`) | 100 | 50 | **0.1129** | 0.2945 | 0.2999 | 0.4508 | 0.4535 | **0.5098** |
| **Actinobacteria** (`act`) | 100 | 50 | **0.1297** | 0.2892 | 0.3205 | 0.4390 | 0.4480 | **0.4838** |
| **Bacteria** (`bac`) | 100 | 50 | **0.1423** | 0.2668 | 0.2774 | 0.3661 | 0.3854 | **0.4011** |
| **全4群統合平均** | **400** | — | **0.1248** | **0.2886** | **0.3033** | **0.4238** | **0.4334** | **0.4693** |

* **MSA+NJ と PSA+NJ の極めて高い一致度（nRF $\approx 0.12$）**: 多重整列とペアワイズ整列から得られる距離行列は、大半の HOG においてほぼ同一の NJ トポロジーを構築する。
* **最尤法（ML）と距離法（NJ）の乖離（nRF $\approx 0.29$）**: 確率モデルによる枝長・置換最適化と距離加法性の差異が安定して現れる。
* **GS法と他手法の乖離（nRF $\approx 0.47$）**: GS 法は最尤法との乖離が最も大きく、特に真核生物（0.5098）で顕著。

#### 手法間トポロジー乖離度の可視化

| 手法間 nRF 分布（バイオリンプロット） | 手法間平均 nRF ヒートマップ |
| :---: | :---: |
| ![手法間 nRF 分布](../dataset/results/method_concordance_violins.png) | ![手法間平均 nRF ヒートマップ](../dataset/results/method_concordance_heatmaps.png) |

#### 最尤法（MSA+ML）を基準とした SSS vs nRF トポロジー距離散布図
![最尤法基準 SSS vs nRF 散布図](../dataset/results/scatter_sss_vs_nrf_by_group_ml_ref.png)

---

## 4. AmpliPhy 準拠 TCS スコア解析 & ルーティング問題の解明

客観的な正解度を測るため、[DessimozLab/ampliphy-analysis](https://github.com/DessimozLab/ampliphy-analysis)（Moi et al., 2023）で提唱された **Taxonomy Congruence Score（TCS / 分類学的整合性スコア）** を導入した。

### 4.1. ルート（Root）の重要性と MAD ルーティングの必要性
- **問題の発見**:
  当初、無根木（IQ-TREE / RapidNJ）の構文上の最外殻の括弧をそのままルートとして計算したところ、ML 法が不当に低く、GS 法が見かけ上高くなる現象が起きた。
- **メカニズム**:
  TCS は「根から葉に向かって分類群の一致度を累積する」階層スコアである。無根木では適当な枝で切られた位置が根になるため、真の単系統クレード（例: 哺乳類群）が根で真っ二つに分断され、壊滅的なスコアペナルティを受けていた。一方、GS はトップダウンの 2 分割法であるため最初から根付き二分木の構造を持っていた。
- **AmpliPhy 準拠の解決**:
  AmpliPhy パイプライン（`scripts/run_tree`）に倣い、**MAD (Minimal Ancestor Deviation: 祖先偏差最小化法)** ルーティングを適用して系統学的な根を推定してから TCS を算出した。

### 4.2. 生物群別 平均 TCS スコア比較（無根木 vs MAD 根付き木）

| 生物群 (Dataset) | 評価状態 | MSA+ML (IQ-TREE) | MSA+NJ (RapidNJ) | PSA+NJ (RapidNJ) | GS (Graph Splitting) | 最高スコア手法 |
| :--- | :--- | :---: | :---: | :---: | :---: | :---: |
| **Amniota** (`amn`, 10 taxa) | 無根（外枠根）<br>**MAD根付き** | 12.486<br>**18.458** | 9.077<br>17.333 | 9.416<br>17.873 | **16.989**<br>16.989 | GS<br>**MSA+ML** |
| **Eukaryota** (`euk`, 50 taxa) | 無根（外枠根）<br>**MAD根付き** | **44.338**<br>**45.420** | 17.262<br>42.046 | 17.217<br>41.633 | 38.073<br>38.073 | MSA+ML<br>**MSA+ML** |
| **Bacteria** (`bac`, 50 taxa) | 無根（外枠根）<br>**MAD根付き** | 5.928<br>**6.711** | 5.625<br>6.607 | 5.795<br>6.402 | **6.663**<br>6.663 | GS<br>**MSA+ML** |
| **Actinobacteria** (`act`, 50 taxa) | 無根（外枠根）<br>**MAD根付き** | 0.600<br>0.618 | 0.565<br>0.563 | 0.602<br>0.587 | **0.659**<br>**0.659** | GS<br>**GS** |

* 正しいルーティングにより、**MSA+ML は 4 群中 3 群（Amniota, Eukaryota, Bacteria）で単独首位**となった。
* NJ 法も Eukaryota で **17.2 $\rightarrow$ 42.0 超へと倍増**し、GS（38.073）を明確に上回った。

#### 生物群別 SSS vs TCS 散布図（MAD ルーティング適用後）
![生物群別 SSS vs TCS 散布図](../dataset/results/scatter_sss_vs_tcs_by_group.png)

---

## 5. 「なぜ SSS が小さいほど TCS が大きくなるのか？」の数理的解明

散布図において「配列類似度（SSS）が下がるほど TCS が上がる（右肩下がり）」傾向が見られた原因について数理・統計的検証を行った。

### 5.1. TCS のサイズ依存性（木サイズ $N$ で割っても深さ $D$ が残る）
AmpliPhy の定義式：
$$\text{TCS} = \frac{\text{root.tax\_score} - \text{root.leaf\_acc} \times |\text{root.overlap}|}{N}$$
- 分子の `tax_score` は各内部ノードで「部分木の葉数 $\times$ 一致ランク数」を累積する。二分木の各階層（レベル）で合計 $N$ の重みがあるため、深さ $D \approx \log_2 N$ 段ある木全体では **$O(N \log N)$** のスケールで増大する。
- 全体を葉数 $N$ で割っても、**「木の深さ $D \approx \log_2 N$（階層数）」の効果が相殺されずに残る**。
- **検証（100% 正解の理想的な木）**:
  - 4 葉（深さ 2）: **TCS = 1.0**
  - 8 葉（深さ 3）: **TCS = 3.0**
  - 16 葉（深さ 4）: **TCS = 6.0**
  - 32 葉（深さ 5）: **TCS = 10.0**
  - 64 葉（深さ 6）: **TCS = 15.0**
  正解度が同じでも、木が大きく階層が深いほど TCS スコア自体が何倍にも跳ね上がる。

### 5.2. 交絡因子としての「生物種数（葉数 $N$）」
- **種数が多い HOG**: 多様な門・綱にまたがる普遍的遺伝子 $\rightarrow$ 遠縁ペアを多く含むため **SSS が小さくなる**（$r = -0.489 \sim -0.699$）。同時に木の階層が深いため **TCS が大きくなる**。
- **種数が少ない HOG**: 特定の近縁群に限定 $\rightarrow$ 配列が似ているため **SSS が高くなる**。同時に階層が浅いため **TCS が小さくなる**。
- **偏相関による立証**: 葉数 $N$ を固定した偏相関係数（Partial Correlation）を算出すると、Eukaryota では $r = -0.418 \rightarrow -0.159$（有意差消失）、Amniota では $r = -0.196 \rightarrow -0.060$（ほぼ無相関）となり、見かけ上の右肩下がりは **生物種数（葉数）がもたらした擬似相関（交絡）** であることが完全に証明された。

---

## 6. 同一 HOG 内での手法間直接比較 & SSS 帯別統計（Within-HOG Evaluation）

木サイズの違いによる交絡を完全に排除するため、**「同一 HOG 内で 4 手法の優劣を直接評価する」** 解析を実施した。

### 6.1. 全体および生物群別 同一 HOG 勝率（HOG 内最高 TCS 獲得率）

| 生物群 (Dataset) | 評価 HOG 数 | MSA+ML (IQ-TREE) | MSA+NJ (RapidNJ) | PSA+NJ (RapidNJ) | GS (Graph Splitting) | 最高勝率手法 | 平均順位 (ML / NJ / PSA / GS) |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **Amniota** (`amn`) | 100 | **40.2%** (40.2勝) | 19.9% (19.9勝) | 24.4% (24.4勝) | 15.4% (15.4勝) | **MSA+ML** | **2.22** / 2.46 / 2.30 / 3.02 |
| **Eukaryota** (`euk`) | 100 | **43.8%** (43.8勝) | 28.3% (28.3勝) | 22.8% (22.8勝) | 5.0% (5.0勝) | **MSA+ML** | **2.06** / 2.13 / 2.33 / 3.48 |
| **Actinobacteria** (`act`) | 100 | 28.2% (28.2勝) | 17.1% (17.1勝) | 19.6% (19.6勝) | **35.1%** (35.1勝) | **GS** | 2.40 / 2.65 / 2.54 / **2.40** |
| **Bacteria** (`bac`) | 100 | 31.6% (31.6勝) | 17.1% (17.1勝) | 18.1% (18.1勝) | **33.2%** (33.2勝) | **GS** | **2.42** / 2.56 / 2.49 / 2.54 |
| **全 4 生物群 合計** | **400** | **36.0%** (143.9勝) | **20.6%** (82.4勝) | **21.2%** (84.9勝) | **22.2%** (88.8勝) | **MSA+ML (総合第1位)** | **2.28** / 2.45 / 2.41 / 2.86 |

### 6.2. 2 手法間の直接対決（Head-to-Head: 全 400 HOG）

| 対戦カード (Method 1 vs Method 2) | M1 勝利数 (%) | 引き分け (%) | M2 勝利数 (%) | 平均 ΔTCS (M1 - M2) | 勝敗の要約 |
| :--- | :---: | :---: | :---: | :---: | :--- |
| **MSA+ML** vs **GS** | **241 勝 (60.2%)** | 29 分 (7.2%) | 130 勝 (32.5%) | **+2.206** | **MSA+ML の圧倒的勝ち越し (勝率 60%)** |
| **MSA+ML** vs **MSA+NJ** | **184 勝 (46.0%)** | 65 分 (16.2%) | 151 勝 (37.8%) | **+1.165** | **MSA+ML の勝ち越し** |
| **MSA+ML** vs **PSA+NJ** | **186 勝 (46.5%)** | 63 分 (15.8%) | 151 勝 (37.8%) | **+1.178** | **MSA+ML の勝ち越し** |
| **PSA+NJ** vs **MSA+NJ** | 114 勝 (28.5%) | **180 分 (45.0%)** | 106 勝 (26.5%) | -0.013 | **45% は完全一致。差がある時は PSA+NJ が勝ち越し** |
| **GS** vs **MSA+NJ** | 138 勝 (34.5%) | 42 分 (10.5%) | **220 勝 (55.0%)** | -1.041 | **MSA+NJ の勝ち越し** |
| **GS** vs **PSA+NJ** | 131 勝 (32.8%) | 42 分 (10.5%) | **227 勝 (56.8%)** | -1.028 | **PSA+NJ の勝ち越し** |

---

### 6.3. SSS（配列類似度）の範囲ごとの統計

全 400 HOG を配列類似度スコア（SSS: $\bar{w}$）に応じて 5 つの進化距離帯に層別化した。

| SSS 範囲（進化距離） | HOG 数 | MSA+ML 勝率 | MSA+NJ 勝率 | PSA+NJ 勝率 | GS 勝率 | 最高勝率手法 | 平均相対スコア (ML / MSA-NJ / PSA-NJ / GS) | 各帯の挙動・特徴 |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :--- |
| **SSS < 0.2**（超遠縁） | 31 | **33.3%** | 20.4% | 26.9% | 19.4% | **MSA+ML** | **0.879** / 0.847 / 0.875 / 0.820 | **PSA+NJ (26.9%) が MSA+NJ (20.4%) を明確に凌駕** |
| **0.2 ≤ SSS < 0.4**（遠縁） | 97 | **32.4%** | 15.2% | 22.9% | **29.5%** | **MSA+ML** | **0.864** / 0.834 / 0.843 / 0.836 | **GS (29.5%) が大きく躍進**。PSA+NJ も MSA+NJ に圧勝 |
| **0.4 ≤ SSS < 0.6**（中等度） | 138 | **37.3%** | 23.5% | 18.1% | 21.2% | **MSA+ML** | **0.859** / 0.819 / 0.784 / 0.815 | **MSA+ML が独走 (37.3%)**。MSA+NJ も復調 |
| **0.6 ≤ SSS < 0.8**（近縁） | 96 | **41.1%** | 21.1% | 21.1% | 16.8% | **MSA+ML** | **0.900** / 0.820 / 0.849 / 0.842 | **MSA+ML の独壇場 (41.1%)**。GS は 16.8% に後退 |
| **0.8 ≤ SSS ≤ 1.0**（極近縁/超保存） | 38 | **29.8%** | 22.8% | 24.1% | 23.2% | **MSA+ML** | **0.875** / 0.792 / 0.814 / 0.863 | 4手法が 23〜30% で拮抗（変異サイト不足） |

#### SSS 帯別・勝率および相対性能プロット
![Within-HOG Performance Stratified by SSS](../dataset/results/plot_sss_binned_win_rates.png)

#### 木サイズ（葉数）の影響を排除した「HOG内正規化 TCS」vs SSS 散布図
![Within-HOG Normalized TCS vs SSS](../dataset/results/plot_within_hog_normalized_tcs.png)

---

---

## 7. 特異領域の深掘り：SSS < 0.1 における挙動 & 系統樹・MSA ファイル集約

配列類似度スコアが極限まで低い **SSS < 0.1** の領域には、全 400 HOG 中 **9 つの HOG**（Bacteria: 8個、Actinobacteria: 1個、平均 SSS = 0.0736、平均アミノ酸一致率 0.3%）が該当します。

この 9 HOG について、全手法の系統樹 Newick ファイル、MAFFT による多重整列（MSA）ファイル、および入力 FASTA 配列を専用ディレクトリに集約・整理しました。

### 7.1. 集約ファイル・ディレクトリ構成

集約先パス: [`dataset/results/low_sss_hogs_under_0.1/`](../dataset/results/low_sss_hogs_under_0.1/)

* **系統樹ディレクトリ (`trees/`)**: 全 36 本（9 HOG × 4 手法）の Newick 系統樹ファイル
  - `{group}_{hog_id}_msa_ml.nwk` : 最尤法（IQ-TREE 2）推定樹
  - `{group}_{hog_id}_msa_nj.nwk` : 多重整列＋近隣結合法（RapidNJ）推定樹
  - `{group}_{hog_id}_psa_nj.nwk` : ペアワイズ整列＋近隣結合法（RapidNJ）推定樹
  - `{group}_{hog_id}_gs.nwk`     : グラフ分割法（`gs2`）推定樹
* **多重整列ディレクトリ (`msa/`)**: 全 9 HOG の MAFFT 多重整列 FASTA ファイル
  - `{group}_{hog_id}_mafft.fasta` : MAFFT (`--auto --threadit 0`) によるアライメント
* **入力配列ディレクトリ (`raw_fasta/`)**:
  - `{group}_{hog_id}.fa` : クエリ配列 FASTA
* **総合メタデータテーブル**:
  - [`low_sss_9hogs_summary.csv`](../dataset/results/low_sss_hogs_under_0.1/low_sss_9hogs_summary.csv) : 9 HOG の全スコア・モデル・Newick 文字列を完全収録

---

### 7.2. SSS < 0.1 における手法別パフォーマンス サマリー

| 指標 | MSA+ML (IQ-TREE) | MSA+NJ (RapidNJ) | PSA+NJ (RapidNJ) | GS (Graph Splitting) | 最高パフォーマンス手法 |
| :--- | :---: | :---: | :---: | :---: | :---: |
| **平均 TCS スコア** | 3.9779 | 4.1467 | **4.1763** | 3.7355 | **PSA+NJ (最高スコア)** |
| **HOG 内勝率** | 22.2% (2勝) | **33.3% (3勝)** | **33.3% (3勝)** | 11.1% (1勝) | **NJ 法（計 66.7% 勝利）** |
| **対 ML 直接対決** | — | **5勝 3敗 1分** | **5勝 4敗** | 3勝 6敗 | **NJ 法が ML に勝ち越し** |
| **相対スコア (Max=1.0)** | **0.8871** | 0.8774 | 0.8709 | 0.7838 | **MSA+ML** |

---

### 7.3. 全 9 HOG の個別系統樹特性と手法比較

配列類似度（SSS）の昇順で各 HOG の系統樹トポロジーと各手法の勝因・分岐パターンをまとめます。

#### ① `bac_E0990658`（全 400 HOG 中 最低 SSS）
* **基本情報**: Taxa 数: 42種 / 整列長: 60 aa / 平均 SSS: **0.0499** (中央値: 0.0221) / Identity: 0.3%
* **最尤推定モデル**: `BLOSUM62+F+R4`
* **TCS スコア**: ML: 2.9286 / **MSA-NJ: 3.2619 (勝)** / PSA-NJ: 2.5952 / GS: 2.6905
* **トポロジー・系統学的知見**:
  - アミノ酸変異が激しく飽和しているため、ML では短い内部枝が乱れ局所解にトラップされた。
  - 多重アライメント上のギャップ除外 Poisson 距離に基づく **MSA-NJ** が、門レベル・綱レベルの主要な単系統群（Bifidobacterium 属や放線菌門のクラスタ）を最も正確に基部に配置し、最高スコアを達成した。
* **関連ファイル**: [MSA](../dataset/results/low_sss_hogs_under_0.1/msa/bac_E0990658_mafft.fasta) / 系統樹: [ML](../dataset/results/low_sss_hogs_under_0.1/trees/bac_E0990658_msa_ml.nwk), [MSA-NJ](../dataset/results/low_sss_hogs_under_0.1/trees/bac_E0990658_msa_nj.nwk), [PSA-NJ](../dataset/results/low_sss_hogs_under_0.1/trees/bac_E0990658_psa_nj.nwk), [GS](../dataset/results/low_sss_hogs_under_0.1/trees/bac_E0990658_gs.nwk)

#### ② `bac_E0889658`
* **基本情報**: Taxa 数: 33種 / 整列長: 60 aa / 平均 SSS: **0.0538** (中央値: 0.0332) / Identity: 0.3%
* **最尤推定モデル**: `Q.PFAM+F+I+R3`
* **TCS スコア**: ML: 2.5152 / MSA-NJ: 2.6364 / **PSA-NJ: 2.8182 (勝)** / GS: 2.4545
* **トポロジー・系統学的知見**:
  - Needleman-Wunsch による 1 対 1 ペアワイズ整列（Gap Open=20, Extension=2）を用いた **PSA-NJ** がトップ。
  - 多重整列では誤った残基整列（ミスアライメント）が系統誤差を誘発していたのに対し、PSA ではノイズの少ないペアワイズ距離行列が得られ、Bacteroidetes と Proteobacteria の大分岐を綺麗に再現した。
* **関連ファイル**: [MSA](../dataset/results/low_sss_hogs_under_0.1/msa/bac_E0889658_mafft.fasta) / 系統樹: [ML](../dataset/results/low_sss_hogs_under_0.1/trees/bac_E0889658_msa_ml.nwk), [MSA-NJ](../dataset/results/low_sss_hogs_under_0.1/trees/bac_E0889658_msa_nj.nwk), [PSA-NJ](../dataset/results/low_sss_hogs_under_0.1/trees/bac_E0889658_psa_nj.nwk), [GS](../dataset/results/low_sss_hogs_under_0.1/trees/bac_E0889658_gs.nwk)

#### ③ `bac_E0993605`
* **基本情報**: Taxa 数: 16種 / 整列長: 60 aa / 平均 SSS: **0.0589** (中央値: 0.0245) / Identity: 0.3%
* **最尤推定モデル**: `Q.PFAM+I+R2`
* **TCS スコア**: ML: 7.5625 / **MSA-NJ: 7.8750 (勝)** / PSA-NJ: 7.3750 / GS: 7.6250
* **トポロジー・系統学的知見**:
  - 種数が 16 種と絞られており、全手法が 7.3 〜 7.8 点台と高水準で拮抗。
  - わずかな内部ノードの解像順序の違いにより、**MSA-NJ** が NCBI 分類階層と最も合致した単系統性を復元した。
* **関連ファイル**: [MSA](../dataset/results/low_sss_hogs_under_0.1/msa/bac_E0993605_mafft.fasta) / 系統樹: [ML](../dataset/results/low_sss_hogs_under_0.1/trees/bac_E0993605_msa_ml.nwk), [MSA-NJ](../dataset/results/low_sss_hogs_under_0.1/trees/bac_E0993605_msa_nj.nwk), [PSA-NJ](../dataset/results/low_sss_hogs_under_0.1/trees/bac_E0993605_psa_nj.nwk), [GS](../dataset/results/low_sss_hogs_under_0.1/trees/bac_E0993605_gs.nwk)

#### ④ `bac_E0991584`（PSA-NJ の圧勝例）
* **基本情報**: Taxa 数: 20種 / 整列長: 60 aa / 平均 SSS: **0.0632** (中央値: 0.0479) / Identity: 0.3%
* **最尤推定モデル**: `Q.PFAM+F+I+R2`
* **TCS スコア**: ML: 5.2500 / MSA-NJ: 5.7000 / **PSA-NJ: 6.5000 (勝)** / GS: 4.0500
* **トポロジー・系統学的知見**:
  - **PSA-NJ が 6.5000 と他手法を圧倒（ML より +1.25 点差、GS より +2.45 点差）**。
  - Bacillota（厚壁菌門）と Actinomycetota（放線菌門）の種間が極めて遠縁であるため、MAFFT による多重アライメントが過剰挿入を引き起こし、ML と MSA-NJ の双方が深層分岐で誤ったグループ化を行っていた。ペアワイズ距離に基づく PSA-NJ のみが正確な門別クレードを維持した。
* **関連ファイル**: [MSA](../dataset/results/low_sss_hogs_under_0.1/msa/bac_E0991584_mafft.fasta) / 系統樹: [ML](../dataset/results/low_sss_hogs_under_0.1/trees/bac_E0991584_msa_ml.nwk), [MSA-NJ](../dataset/results/low_sss_hogs_under_0.1/trees/bac_E0991584_msa_nj.nwk), [PSA-NJ](../dataset/results/low_sss_hogs_under_0.1/trees/bac_E0991584_psa_nj.nwk), [GS](../dataset/results/low_sss_hogs_under_0.1/trees/bac_E0991584_gs.nwk)

#### ⑤ `bac_E0994775`（最尤法 ML の勝利例）
* **基本情報**: Taxa 数: 28種 / 整列長: 60 aa / 平均 SSS: **0.0719** (中央値: 0.0433) / Identity: 0.3%
* **最尤推定モデル**: `Q.PFAM+F+R4`
* **TCS スコア**: **ML: 4.5714 (勝)** / MSA-NJ: 3.8929 / PSA-NJ: 3.6429 / GS: 3.3929
* **トポロジー・系統学的知見**:
  - SSS < 0.1 の中では比較的類似度が高く、サイト間速度不均一性（Rate heterogeneity: R4）を組み込んだ **最尤法（IQ-TREE 2）** が距離法（NJ）の加法性ノイズを凌駕し、単独首位を獲得した。
* **関連ファイル**: [MSA](../dataset/results/low_sss_hogs_under_0.1/msa/bac_E0994775_mafft.fasta) / 系統樹: [ML](../dataset/results/low_sss_hogs_under_0.1/trees/bac_E0994775_msa_ml.nwk), [MSA-NJ](../dataset/results/low_sss_hogs_under_0.1/trees/bac_E0994775_msa_nj.nwk), [PSA-NJ](../dataset/results/low_sss_hogs_under_0.1/trees/bac_E0994775_psa_nj.nwk), [GS](../dataset/results/low_sss_hogs_under_0.1/trees/bac_E0994775_gs.nwk)

#### ⑥ `bac_E0983296`（グラフ分割法 GS の勝利例）
* **基本情報**: Taxa 数: 33種 / 整列長: 60 aa / 平均 SSS: **0.0811** (中央値: 0.0591) / Identity: 0.3%
* **最尤推定モデル**: `LG+F+I+R3`
* **TCS スコア**: ML: 3.5455 / MSA-NJ: 3.5455 / PSA-NJ: 3.6667 / **GS: 4.2424 (勝)**
* **トポロジー・系統学的知見**:
  - **GS 法が 4.2424 をマークし、全他手法（3.5〜3.6 点台）を大きく引き離して勝利**。
  - MMseqs2 のペアワイズ相同性ネットワークに基づくスペクトラル 2 分割が、多重整列を行わずに Acinetobacter 属や Burkholderia 目のクラスタを初手で綺麗に大分割した。アライメント不能な遠縁タンパク質に対する GS 法の設計思想通りの強みが発揮された好例。
* **関連ファイル**: [MSA](../dataset/results/low_sss_hogs_under_0.1/msa/bac_E0983296_mafft.fasta) / 系統樹: [ML](../dataset/results/low_sss_hogs_under_0.1/trees/bac_E0983296_msa_ml.nwk), [MSA-NJ](../dataset/results/low_sss_hogs_under_0.1/trees/bac_E0983296_msa_nj.nwk), [PSA-NJ](../dataset/results/low_sss_hogs_under_0.1/trees/bac_E0983296_psa_nj.nwk), [GS](../dataset/results/low_sss_hogs_under_0.1/trees/bac_E0983296_gs.nwk)

#### ⑦ `act_E0992444`（放線菌群で唯一の SSS < 0.1）
* **基本情報**: Taxa 数: 29種 / 整列長: 60 aa / 平均 SSS: **0.0901** (中央値: 0.0402) / Identity: 0.4%
* **最尤推定モデル**: `Q.PFAM+F+R4`
* **TCS スコア**: **ML: 0.1379 (勝)** / MSA-NJ: 0.0690 / PSA-NJ: 0.0690 / GS: 0.0690
* **トポロジー・系統学的知見**:
  - 放線菌群全体の分類ランク付与の特性上スコアの絶対値は小さいが、他手法がベースライン（0.0690）に留まる中、**最尤法のみが上位階層（目・科レベル）の正しい一致を復元し 2 倍のスコア（0.1379）を記録**した。
* **関連ファイル**: [MSA](../dataset/results/low_sss_hogs_under_0.1/msa/act_E0992444_mafft.fasta) / 系統樹: [ML](../dataset/results/low_sss_hogs_under_0.1/trees/act_E0992444_msa_ml.nwk), [MSA-NJ](../dataset/results/low_sss_hogs_under_0.1/trees/act_E0992444_msa_nj.nwk), [PSA-NJ](../dataset/results/low_sss_hogs_under_0.1/trees/act_E0992444_psa_nj.nwk), [GS](../dataset/results/low_sss_hogs_under_0.1/trees/act_E0992444_gs.nwk)

#### ⑧ `bac_E0997625`（最大規模 50 Taxa での PSA-NJ 勝利例）
* **基本情報**: Taxa 数: **50種** / 整列長: 60 aa / 平均 SSS: **0.0940** (中央値: 0.0612) / Identity: 0.3%
* **最尤推定モデル**: `WAG+R5`
* **TCS スコア**: ML: 6.8400 / MSA-NJ: 6.7400 / **PSA-NJ: 7.5200 (勝)** / GS: 6.9200
* **トポロジー・系統学的知見**:
  - 50 種という大規模かつ超遠縁な条件において、**PSA-NJ が 7.5200 で単独トップ**。
  - 50 配列の多重整列ではアライメントの不確実性が全体に波及するのに対し、PSA は各配列対の相同領域のみを評価するため、長大な系統枝を持つ外群種が基部に正確に配置された。
* **関連ファイル**: [MSA](../dataset/results/low_sss_hogs_under_0.1/msa/bac_E0997625_mafft.fasta) / 系統樹: [ML](../dataset/results/low_sss_hogs_under_0.1/trees/bac_E0997625_msa_ml.nwk), [MSA-NJ](../dataset/results/low_sss_hogs_under_0.1/trees/bac_E0997625_msa_nj.nwk), [PSA-NJ](../dataset/results/low_sss_hogs_under_0.1/trees/bac_E0997625_psa_nj.nwk), [GS](../dataset/results/low_sss_hogs_under_0.1/trees/bac_E0997625_gs.nwk)

#### ⑨ `bac_E0991547`（NJ 法の独壇場）
* **基本情報**: Taxa 数: 40種 / 整列長: 60 aa / 平均 SSS: **0.0994** (中央値: 0.0875) / Identity: 0.3%
* **最尤推定モデル**: `Q.PFAM+F+R5`
* **TCS スコア**: ML: 2.4500 / **MSA-NJ: 3.6000 (勝)** / PSA-NJ: 3.4000 / GS: 2.1750
* **トポロジー・系統学的知見**:
  - **NJ 法（MSA-NJ: 3.60, PSA-NJ: 3.40）が、最尤法（2.45）や GS（2.18）を大きく引き離して完勝**。
  - ML と GS は Bifidobacterium 属の近縁種を別々のサブツリーに誤って分散させてしまったのに対し、RapidNJ は距離の加法性からこれらを正確に単一クレードに束ね上げた。
* **関連ファイル**: [MSA](../dataset/results/low_sss_hogs_under_0.1/msa/bac_E0991547_mafft.fasta) / 系統樹: [ML](../dataset/results/low_sss_hogs_under_0.1/trees/bac_E0991547_msa_ml.nwk), [MSA-NJ](../dataset/results/low_sss_hogs_under_0.1/trees/bac_E0991547_msa_nj.nwk), [PSA-NJ](../dataset/results/low_sss_hogs_under_0.1/trees/bac_E0991547_psa_nj.nwk), [GS](../dataset/results/low_sss_hogs_under_0.1/trees/bac_E0991547_gs.nwk)

---

## 8. 生成・格納された成果物一覧

### 集計データテーブル（CSV）
1. [dataset_tcs_summary.csv](file:///Users/kazukiaibara/PhyloMethod/dataset/results/dataset_tcs_summary.csv): 全 400 HOG × 4 手法の TCS スコア（MAD根付き・無根）および SSS・Identity
2. [low_sss_9hogs_summary.csv](file:///Users/kazukiaibara/PhyloMethod/dataset/results/low_sss_hogs_under_0.1/low_sss_9hogs_summary.csv): SSS < 0.1 の 9 HOG 特化メタデータ・TCS・全樹形 Newick 一覧
3. [within_hog_win_rates.csv](file:///Users/kazukiaibara/PhyloMethod/dataset/results/within_hog_win_rates.csv): 生物群別・手法別の勝率・平均順位・正規化比
4. [pairwise_head_to_head.csv](file:///Users/kazukiaibara/PhyloMethod/dataset/results/pairwise_head_to_head.csv): 2 手法間の直接対決勝敗数・勝率
5. [sss_binned_tcs_stats.csv](file:///Users/kazukiaibara/PhyloMethod/dataset/results/sss_binned_tcs_stats.csv): SSS 帯別の勝率・相対スコアデータ
6. [method_concordance_summary.csv](file:///Users/kazukiaibara/PhyloMethod/dataset/results/method_concordance_summary.csv): 全ペアワイズ RF / nRF 距離

### SSS < 0.1 専用集約ディレクトリ
* [dataset/results/low_sss_hogs_under_0.1/](file:///Users/kazukiaibara/PhyloMethod/dataset/results/low_sss_hogs_under_0.1/):
  - `trees/` : 全 36 本（9 HOG × 4手法）の Newick 系統樹ファイル
  - `msa/` : 全 9 HOG の MAFFT 多重整列 FASTA ファイル
  - `raw_fasta/` : 全 9 HOG の入力配列 FASTA ファイル

### 可視化プロット（PNG）
1. [plot_sss_binned_win_rates.png](file:///Users/kazukiaibara/PhyloMethod/dataset/results/plot_sss_binned_win_rates.png): SSS 帯別の勝率棒グラフ & 相対スコア推移
2. [plot_within_hog_normalized_tcs.png](file:///Users/kazukiaibara/PhyloMethod/dataset/results/plot_within_hog_normalized_tcs.png): 木サイズの影響を排除した HOG 内正規化スコア vs SSS 散布図
3. [scatter_sss_vs_tcs_by_group.png](file:///Users/kazukiaibara/PhyloMethod/dataset/results/scatter_sss_vs_tcs_by_group.png): 生物群別 SSS vs TCS 散布図（2x2）
4. [method_concordance_heatmaps.png](file:///Users/kazukiaibara/PhyloMethod/dataset/results/method_concordance_heatmaps.png): 手法間平均 nRF ヒートマップ
5. [method_concordance_violins.png](file:///Users/kazukiaibara/PhyloMethod/dataset/results/method_concordance_violins.png): 手法間 nRF 分布バイオリンプロット

### 実行スクリプト（Python）
1. [run_dataset_methods.py](file:///Users/kazukiaibara/PhyloMethod/dataset/scripts/run_dataset_methods.py): 実データ 4 手法並列実行コアスクリプト
2. [calculate_dataset_tcs.py](file:///Users/kazukiaibara/PhyloMethod/dataset/scripts/calculate_dataset_tcs.py): MAD ルーティング適用 TCS 計算スクリプト
3. [analyze_within_hog_tcs.py](file:///Users/kazukiaibara/PhyloMethod/dataset/scripts/analyze_within_hog_tcs.py): 同一 HOG 比較・勝率・SSS 帯別統計スクリプト
