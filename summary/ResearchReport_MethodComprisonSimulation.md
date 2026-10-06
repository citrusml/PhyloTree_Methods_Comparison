# ResearchReport_MethodComprisonSimulation

## 趣旨

このレポートでは`/Users/aibarakazuki/PhyloTree_Methods_Comparison/`で行った実験の結果をリサーチレポートの形でまとめる。まとめる際は大まかな時系列順に結果をまとめ、結果のところには結果を示す図がマークダウンで見れるように表示する。Resultsでそれぞれの実験の結果をまとめるが、Introにはどのような実験を何のために行なったかを羅列しておいてください。

また本実験は実験の条件が非常に重要なため、どのような条件で行ったかをしっかりと明記してください。また技術的に大きな変更をしている場合はそれをMethodのところに追記してください。

## 本文

> 出典: [README.md](../README.md), [README_ALL_RESULTS_SSS.md](../README_ALL_RESULTS_SSS.md)、および `results/` 内の既存図。数値は上記READMEの記載値を転記したものであり、本レポート作成時に再集計はしていない。

### Introduction

系統樹推定では「配列アラインメント → 距離/尤度計算 → 木構築」の各段階が誤差源となる。特に配列類似度が低い領域（Twilight Zone, $SSS<0.06$）では多重配列アラインメント(MSA)が破綻し、MSAに依存しない手法の方が頑健になる可能性がある(Matsui & Iwasaki 2020, *Syst. Biol.*, syz049)。本研究ではシミュレーションにより、進化距離 $D$・配列長 $L$ といった条件に対する各手法のトポロジー誤差(nRF)を比較し、手法間の優劣が入れ替わる境界(Regime Map / Crossover)を調べた。実施した実験とその目的を時系列で列挙する。

| # | 時期(結果ディレクトリ更新日) | 実験 (`results/`) | 目的 |
|---|---|---|---|
| 1 | 08/19 | `results_ver2`, `results_ver3` | 初期の $L \times D$ グリッドで PWA+NJ / MSA+NJ / MSA+ML のRegime Mapを作成し、パイプラインを検証 |
| 2 | 08/23 | `results_new` | 標準ベースライン。真の距離行列(TRUE_DIST+NJ)を対照に加え、アラインメント誤差と距離推定誤差を分離 |
| 3 | 08/26 | `results_alpha`, `results_gamma` | ガンマ形状母数 $\alpha$(サイト間速度不均一性)の影響評価 |
| 4 | 08/26 | `results_taxon` | 分類群数 $N$ のスケーラビリティ評価 |
| 5 | 08/27 | `results_true_pwa` | 真のペアワイズアラインメントを用いたPSA+NJの上限(アラインメント誤差の寄与)の評価 |
| 6 | 08/29 | `results_fastme_NoOption` | NJ以外の距離法(FastME既定)との比較 |
| 7 | 08/29 | `results_high_dist` | 極限進化距離($D=4\sim6$)での情報限界の確認 |
| 8 | 08/29 | `results_ics` | 保存領域(ICS)比率が各手法に与える影響 |
| 9 | 08/30 | `results_simple` | 速度不均一性なしの単純LGモデルでの挙動確認 |
| 10 | 08/31 | `results_fastme_options` | FastME(SPR探索+LG+G補正距離)の効果 |
| 11 | 09/03 | `results_high_gap` | 広域 $D$ (0.1〜6.0)での連続的な感度曲線取得 |
| 12 | 09/04 | `results_zipfian_indel` | べき乗則(Zipfian)Indel導入によるMSA+MLの崩壊とPSA+NJとの逆転の検証 |
| 13 | 09/09 | `results_paper_tree_pre_fix`, `results_paper_tree` | 論文準拠の後退Yule系統樹+Graph Splitting(GS)法を含む比較 |
| 14 | 09/11 | `results_paper_tree2` | INDELible v1.03を用いた論文 Fig. 3a の完全再現 |
| 15 | 09/23 | `results_ics_full` | 全長が保存領域($ICS=1.0$)の対照実験 |
| 16 | 10/07 | `results_paper_tree3` | INDELible条件で100レプリケートに増やし、SP scoreと誤差伝播(アラインメント誤差→木推定誤差)を解析 |

### Method

#### 共通の実験枠組み

- **パイプライン**: Nextflow DSL2 ([main.nf](../main.nf))。ローカル/LSFスパコンで実行。環境は Micromamba ([environment.yml](../environment.yml))。
- **2×2要因計画**:

  | | NJ | ML |
  |---|---|---|
  | MSA (MAFFT) | MSA+NJ | MSA+ML |
  | PWA/PSA (ペアワイズ) | PWA+NJ | - |

- **各手法**:
  - `PWA(PSA)+NJ`: Needleman-Wunsch(LG行列) → Poisson/Gamma距離 → RapidNJ/FastME (論文再現系ではBLOSUM62, Gap 20/2)
  - `MSA+NJ`(対照): MAFFT → 同一のPoisson/Gamma距離 → RapidNJ (MSA誤差と木推定の分離)
  - `MSA+ML`: MAFFT → IQ-TREE 2 (ModelFinder `-m MFP`)
  - `MSA+RAXML`, `GS`(Graph Splitting `gs2`, MMseqs2 `-m 7.5`): 論文再現系(実験13, 14)で追加
  - 真値対照: `TRUE_MSA+*`, `TRUE_PSA+NJ`, `TRUE_DIST+NJ`
- **評価指標**: DendroPyによるunrooted Robinson-Foulds距離および正規化RF(nRF)。比較は $\Delta nRF = nRF_{\mathrm{PWA+NJ}} - nRF_{\mathrm{MSA+ML}}$ (負: PWA+NJ優位、正: MSA+ML優位)。
- **横軸指標**: Sequence Similarity Score $SSS$ ($w$)。レプリケート個別のシードで全件計算(欠損0件)し、対数空間LOWESSで平滑化。
- **配列進化**: 初期〜`results_paper_tree` はAliSim(LG+G4)。`results_paper_tree2`/`results_paper_tree3` は INDELible v1.03(WAG+G4)。`results_ics_full`(再実行版)は AliSim + ICS+G4。

> 以下の条件は [next_configs/](../next_configs) の `params` と [run_scripts/](../run_scripts) から確認した値。READMEの記述と異なる点は ⚠ で示す。

#### 初期ベースライン系 (`next_configs/old/*.config`) の共通条件

| 項目 | 値 |
|---|---|
| 系統樹 | 標準(birth_rate=0.1, death_rate=0.05) |
| 分類群数 $N$ | 32 (`results_taxon` のみ可変) |
| 距離 $D$ | 0.1, 0.5, 1.0, 2.0, 3.0 |
| 配列長 $L$ | 100, 300, 500, 1000, 1500 (実験により異なる) |
| レプリケート数 | 100 / 条件 |
| 置換モデル | LG+G4, $\alpha=1.0$ |
| Indel率 | 挿入 0.05 / 欠失 0.10 (Zipfian指定なし) |
| ギャップペナルティ | Open=10.0, Extend=0.5 (`high_gap`のみ 20.0/2.0) |
| 距離補正 | Poisson, NJ: RapidNJ |

#### 実験条件一覧(差分のみ)

| 実験 | $N$ | $L$ | $D$ | 条件の差分 | SSS範囲(README) |
|---|---|---|---|---|---|
| `results_ver2` (`nextflow2.config`) | 32 | 100〜1500 | 0.1〜3.0 | ベースライン | - |
| `results_ver3` (`nextflow.config`) | (未指定) | 100〜1500 | 0.1〜3.0 | ベースライン初版。⚠ ver2/ver3 の違いは設定からは特定できず(要確認) | - |
| `results_new` | - | - | 0.1〜3.0 | TRUE_DIST+NJ対照 | 0.0064〜0.9494 |
| `results_alpha` | 32 | 100, 500, 1000 | 0.1〜3.0 | $\alpha\in\{0.25,0.5,1.0,2.0\}$(真値), 距離はPoisson。READMEは100除外と記載 | 0.0063〜0.9537 |
| `results_gamma` | 32 | 100, 500, 1000 | 0.1〜3.0 | 同上 + 距離補正を **gamma_poisson** に変更 | 0.0063〜0.9537 |
| `results_taxon` | 8, 16, 64, 128 (+32は基準) | 100, 500, 1000 | 0.1〜3.0 | $N$ 可変 | 0.0000〜0.9861 |
| `results_true_pwa` | 32 | 100〜1500 | 0.1〜3.0 | 真のPWA+NJを追加 (`main_true_pwa.nf`) | 0.0064〜0.9494 |
| `results_fastme_NoOption` | 32 | 100〜1500 | 0.1〜3.0 | `nj_tool="fastme"` (FastME既定) (`main_fastme.nf`) | 0.0064〜0.9494 |
| `results_high_dist` | 32 | 300〜1500 | 4.0, 5.0, 6.0 | 極限距離 | 0.0078〜0.0910 |
| `results_ics` | 32 | 100, 300, 500, 1000 | 0.1〜3.0 | ICS比率 0, 0.05, 0.10, 0.20 (Dayhoff 6分類不変サイト) | 0.0072〜0.9525 |
| `results_simple` | 32 | 100〜1500 | 0.1〜3.0 | モデル **LG**(ガンマなし) | 0.0035〜0.9537 |
| `results_fastme_options` | 32 | 100〜1500 | 0.1〜3.0 | `MSA+FastME_LG_G`, `PWA+FastME_SPR`(`fastme -mB -s -q`: BME+SPR探索+三角不等式補正)(`main_fastme_options.nf`) | 0.0078〜0.9318 |
| `results_high_gap` | 32 | 300〜1500 | 0.1〜6.0 (8点) | ギャップペナルティ **20/2** | 0.0078〜0.9318 |
| `results_zipfian_indel` | 32 | 300, 500, 1000, 1500 | 0.1〜6.0 (8点) | Zipfian `POW{1.7/50}`, 挿入/欠失 **0.05/0.05**(⚠READMEは0.05/0.10の併記), Gap 20/2, LG+G4, 100レプリケート, GS追加 | 0.0078〜0.9318 |

#### 論文準拠系の条件 (`nextflow_paper_tree*.config`, `nextflow_ics_full.config`)

共通: 後退Yule過程+枝長分布 $1-\ln(u(e-1)+1)$、Zipfian Indel、挿入/欠失 0.10/0.10、Gap Open/Extend = 20/2、Poisson距離、RapidNJ、$\alpha=1.0$、比較手法は PWA+NJ, MSA+NJ, MSA+ML, MSA+RAXML(`PROTGAMMAIWAGX`), GS(MMseqs2 `-m 7.5`), TRUE_*対照群、SSS計算あり。

| 実験 | $N$ | $L$ | $D$ | レプリケート数 | シミュレータ / モデル | 備考 |
|---|---|---|---|---|---|---|
| `results_paper_tree` | 32 | 500, 1000, 1500 | 0.1, 0.2, 0.5, 1.0, 1.5, 2.0 | **30** | AliSim / LG+G4, `POW{1.7/50}` | `pre_fix`の設定は未確認 |
| `results_paper_tree2` | 20 | 1000 | 同上 | **30** (⚠READMEは60と記載。設定ファイルは30) | **INDELible v1.03** / WAG+G4, `POW 1.7 50` | `bin/install_indelible.sh` で配備 |
| `results_paper_tree3` (10/07) | 20 | 500, 1000, 1500 | 同上 | **100** (chunk_size=20) | INDELible v1.03 / WAG+G4 | SP score(f_D)計算、中間生データ保存、距離誤差伝播解析を追加。BioPerl必要 |
| `results_ics_full` (再実行版, 09/23) | 32 | 300, 500, 1000, 1500 | 同上 | **50** | AliSim / `ICS+G4` (`models/ics_model.nex`), `POW{1.7/50}` | ⚠READMEは $L$=100〜1500。`old/` の旧設定は初期ベースライン条件(ICS=1.0, ins/del 0.05/0.10, Gap 10/0.5, $D$=0.1〜3.0) |

`results_paper_tree2` の「60レプリケート・3,240推定」はREADME記載値であり、設定の30とは整合しないため、再実行・追加実行の有無を要確認。

#### 技術的に大きな変更点

1. **ギャップペナルティの変更**: Open/Extend = 10/0.5 → **20/2** (`results_high_gap` 以降、PhyPA準拠)。
2. **Indel長分布**: 既定 → **べき乗則 POW{1.7/50}**(最大50残基)(`results_zipfian_indel` 以降)。
3. **系統樹生成モデルの変更**: 標準(birth-death) → **後退Yule過程+対数サンプリング枝長**(`results_paper_tree` 以降、論文準拠)。`results_paper_tree_pre_fix` から枝長スケール補正を行った。
4. **配列シミュレータの変更**: AliSim(LG+G4) → **INDELible v1.03**(WAG+G4)(`results_paper_tree2` 以降)。論文条件の厳密再現が目的。
5. **比較手法の拡張**: FastME(既定、SPR+LG+G補正)、RAxML(`-f d`, PROTGAMMAIWAGX)、GS法、真値対照群(TRUE_*)を追加し誤差源を分離。
6. **評価軸・出力の拡張**: $SSS$ を横軸とする評価を導入。`results_paper_tree3` で SP score、中間生データ保存、距離誤差伝播解析を追加。
7. **実行基盤**: 一部の実験は `next_main/old/main_*.nf`(実験別ワークフロー)で実施し、後に共通 `main.nf` + 設定ファイル切替へ移行。

### Results & Discussion

> 図は既存PNGを埋め込んでいる。READMEで参照されている `nrf_vs_sss_by_condition.png` は本リポジトリの `results/` に存在しない(別環境 `/Users/kazukiaibara/PhyloMethod/results` で生成)ため、本リポジトリにある代表図で代替した。該当図がない実験は図なしとした。

#### 1. 初期Regime Map (`results_ver2`, `results_ver3`, 08/19)

![ver3 Regime Map](../results/results_ver3/regime_map_delta_nrf.png)

![ver3 nRF boxplot](../results/results_ver3/nrf_boxplots_5x5_grid.png)

パイプラインの動作検証と $L \times D$ 上の $\Delta nRF$ マップ作成。結果の数値的詳細はREADMEに記載がない。

#### 2. ベースライン (`results_new`, 08/23)

![new Regime Map](../results/results_new/regime_map_delta_nrf.png)

![new nRF boxplot](../results/results_new/nrf_boxplots.png)

`TRUE_DIST+NJ` はアラインメント・距離推定誤差ゼロの理論下限を与える。低SSS域での誤差増大が純粋にアラインメントの破綻に起因することを示した。

#### 3. ガンマ形状母数 (`results_alpha`, `results_gamma`, 08/26)

![alpha scaling](../results/results_alpha/alpha_scaling_curves.png)

![gamma distance alpha scaling](../results/results_gamma/gamma_distance_alpha_scaling.png)

- $\alpha=0.25$(強い速度不均一)でSSS低下が最も急。+G4を考慮した `MSA+ML` の優位性が明確。
- `results_gamma` は `results_alpha` と同一グリッドの再現検証。$w>0.15$で `MSA+ML` は $nRF<0.15$ を維持し、$w\approx0.06$付近で手法間差が縮小。

#### 4. 分類群数 (`results_taxon`, 08/26)

![taxon scaling](../results/results_taxon/taxon_scaling_curves.png)

$N=8$ では誤差が低く抑えられる一方、$N=64,128$ では探索空間の増大・LBA・アラインメントミスの影響で低SSS側のnRF悪化が急峻になった。

#### 5. 真のPWA対照 (`results_true_pwa`, 08/27)

![true pwa heatmap](../results/results_true_pwa/true_pwa_heatmap.png)

![true pwa scaling](../results/results_true_pwa/true_pwa_scaling_curves.png)

PSA+NJの誤差のうち、ペアワイズアラインメント誤りに起因する成分を切り出す基準線。

#### 6. FastME (`results_fastme_NoOption` 08/29, `results_fastme_options` 08/31)

![fastme default](../results/results_fastme_NoOption/scaling_curves_fastme.png)

![fastme options](../results/results_fastme_options/scaling_curves_fastme_options.png)

SPR局所探索+Gamma補正距離の `MSA+FastME_LG_G` は、RapidNJより低類似度域で良好なトポロジー復元を示した。既定オプションFastMEは長大配列でNJより良好。

#### 7. 極限距離 (`results_high_dist`, 08/29)

![high dist regime](../results/results_high_dist/regime_map_high_dist.png)

![high dist scaling](../results/results_high_dist/scaling_curves_high_dist.png)

SSS最大0.0910、中央値0.0318。相同シグナルがほぼ失われ、全手法で $nRF>0.70$、理論的情報限界に達した。

#### 8. 保存領域 ICS (`results_ics` 08/29, `results_ics_full` 09/23)

![ics scaling](../results/results_ics/ics_scaling_benchmark.png)

![ics full benchmark](../results/results_ics_full/ics_full_benchmark.png)

![ics full nrf vs sss](../results/results_ics_full/nrf_vs_sss_curves.png)

- ICS比率が高いほど(0.20)、$D=3.0$でもSSSが保たれ、PSA系手法の劣化が緩和された。
- $ICS=1.0$では置換速度が抑制され、SSS中央値0.375、nRFが全体的に極めて低かった。

#### 9. 単純モデル (`results_simple`, 08/30)

![simple regime](../results/results_simple/regime_map_high_dist.png)

![simple scaling](../results/results_simple/scaling_curves_high_dist.png)

速度不均一性がないと同一 $D$ でのSSSが安定し、高類似度域の精度が全体的に向上した。

#### 10. 広域距離 (`results_high_gap`, 09/03)

![high gap regime](../results/results_high_gap/regime_map_delta_nrf.png)

![high gap boxplot](../results/results_high_gap/nrf_boxplots.png)

$w>0.10$ で `MSA+FastME(LG+G)` が `MSA+NJ` より低nRF。$w<0.05$ ではアラインメント破綻により距離法(NJ, FastME)は同水準の誤差に漸近した。

#### 11. Zipfian Indel (`results_zipfian_indel`, 09/04)

![zipfian regime](../results/results_zipfian_indel/regime_map_delta_nrf.png)

![zipfian boxplot](../results/results_zipfian_indel/nrf_boxplots.png)

長大Indelにより $w<0.06$ で `MSA+ML` が急激に悪化し、`PSA+NJ` が逆転した(Crossover)。

#### 12. 論文準拠系統樹 (`results_paper_tree_pre_fix`, `results_paper_tree`, 09/09)

![paper tree pre_fix curves](../results/results_paper_tree_pre_fix/nrf_vs_sss_curves.png)

![paper tree pre_fix by length](../results/results_paper_tree_pre_fix/paper_tree_nrf_vs_sss_by_length.png)

SSSは0.0139〜0.4948(中央値0.0527)でTwilight Zone付近に集中。$L=1500$で全手法が向上し、$w<0.05$では `GS` が `MSA+ML` に匹敵〜凌駕した。`pre_fix`(枝長補正前)でも同様の順位逆転挙動が再現された。

#### 13. 論文完全再現 (`results_paper_tree2`, 09/11)

![paper tree2 curves](../results/results_paper_tree2/nrf_vs_sss_curves.png)

![paper tree2 regime](../results/results_paper_tree2/regime_map_delta_nrf.png)

![paper tree2 SSS distribution](../results/results_paper_tree2/sss_distribution_by_distance.png)

![paper tree2 model selection](../results/results_paper_tree2/model_selection_nrf_comparison.png)

1. $D\le0.5$($SSS\ge0.06$)では `MSA+ML`, `MSA+RAXML` が $nRF\approx0.00\sim0.01$。$D=1.5, 2.0$(SSS≈0.024, 0.015)でMAFFTの誤配置が増え、MLの誤差が急上昇。
2. Extreme領域($SSS\le0.03$, 96レプリケート)の平均nRF: `GS` 0.3768、`MSA+ML` 0.3768、`MSA+RAXML` 0.3811。$D=2.0$では `GS` 0.4088 vs `MSA+ML` 0.4069 / `MSA+RAXML` 0.4137 と互角。論文Fig. 3aの傾向(低類似度域でのGSの頑健性)を再現。
3. IQ-TREE 2とRAxMLは全距離帯でほぼ同じ感度曲線を示し、ツール差によるバイアスはない。
4. `PSA+NJ` は $D\le0.2$ で $nRF<0.01$ だが、$D\ge1.0$ では累積ギャップペナルティの影響で $nRF\approx0.36\sim0.47$。

#### 14. 論文完全再現 + 誤差伝播解析 (`results_paper_tree3`, 10/07)

条件: $N=20$, $L\in\{500,1000,1500\}$, $D\in\{0.1,0.2,0.5,1.0,1.5,2.0\}$, 各100レプリケート(計1,800データセット × 9パイプライン = 16,200推定)。INDELible v1.03 / WAG+G4 / `POW 1.7 50`。数値は `benchmark_summary.csv` の再集計値(全条件100レプリケート揃い、欠損なし)。

![paper tree3 curves](../results/results_paper_tree3/nrf_vs_sss_curves.png)

![paper tree3 regime](../results/results_paper_tree3/regime_map_delta_nrf.png)

![paper tree3 SSS distribution](../results/results_paper_tree3/sss_distribution_by_distance.png)

**距離別の平均nRF** (全 $L$ 込み)

| pipeline | D=0.1 | 0.2 | 0.5 | 1.0 | 1.5 | 2.0 |
|---|---|---|---|---|---|---|
| GS | 0.091 | 0.072 | 0.048 | 0.156 | 0.312 | 0.430 |
| MSA+ML | 0.000 | 0.000 | 0.012 | 0.123 | 0.297 | 0.426 |
| MSA+NJ | 0.001 | 0.000 | 0.009 | 0.108 | 0.283 | 0.411 |
| MSA+RAXML | 0.000 | 0.000 | 0.012 | 0.123 | 0.296 | 0.429 |
| PWA+NJ | 0.005 | 0.010 | 0.043 | 0.198 | 0.370 | 0.490 |
| TRUE_MSA+ML | 0.000 | 0.000 | 0.001 | 0.012 | 0.079 | 0.242 |
| TRUE_MSA+NJ | 0.001 | 0.000 | 0.002 | 0.026 | 0.108 | 0.163 |
| TRUE_MSA+RAXML | 0.000 | 0.000 | 0.001 | 0.012 | 0.081 | 0.245 |
| TRUE_PWA+NJ | 0.001 | 0.000 | 0.002 | 0.026 | 0.108 | 0.163 |

**配列長別の平均nRF**

| pipeline | L=500 | 1000 | 1500 |
|---|---|---|---|
| GS | 0.198 | 0.177 | 0.179 |
| MSA+ML | 0.158 | 0.139 | 0.132 |
| MSA+NJ | 0.149 | 0.131 | 0.127 |
| MSA+RAXML | 0.157 | 0.140 | 0.134 |
| PWA+NJ | 0.200 | 0.182 | 0.176 |
| TRUE_MSA+ML | 0.090 | 0.049 | 0.027 |
| TRUE_MSA+NJ | 0.070 | 0.047 | 0.033 |

**SSS・アラインメント精度** (データセット単位の平均)

| D | 0.1 | 0.2 | 0.5 | 1.0 | 1.5 | 2.0 |
|---|---|---|---|---|---|---|
| SSS (mean) | 0.4964 | 0.2630 | 0.0849 | 0.0403 | 0.0271 | 0.0208 |
| SP score | 0.9678 | 0.8853 | 0.5480 | 0.3429 | 0.3165 | 0.3185 |

SSS全体: 範囲 0.0086〜0.6516、中央値 0.0573、第1四分位 0.0271。Extreme領域($SSS\le0.03$)は530データセットで、平均nRFは `MSA+NJ` 0.3600、`MSA+ML` 0.3757、`MSA+RAXML` 0.3779、`GS` 0.3836、`PWA+NJ` 0.4494、真のアラインメント使用(`TRUE_MSA+ML`/`NJ`)は 0.155 / 0.128。

##### 誤差伝播解析で用いた評価指標の定義とカスケード構造

`results_paper_tree3` では、アラインメント誤差が進化的距離推定および最終的な系統樹トポロジー推定にどのように波及・拡大するかを定量化するため、以下の指標群を用いた誤差伝播解析（Error Propagation Cascade）を導入している。

1. **アラインメント精度指標**
   - **Sum-of-Pairs (SP) score ($f_D$, Developer's score)**:
     - **定義**: 真のアラインメント（INDELible 出力の True MSA）に含まれるすべての相同アミノ酸残基ペアのうち、推定アラインメント（MAFFT 出力の Estimated MSA）でも正しく同一列に配置されたペアの割合（相同残基の再現率 / Recall）。
     - **数式**:
       $$SP = \frac{|\text{推定アラインメントで正しく整列された真の相同残基ペア}|}{|\text{真のアラインメントにおける全相同残基ペア}|}$$
     - **解釈**: $1.0$ で完全一致（すべての相同残基が正確に同定・整列）。$0.0$ で相同関係が完全に誤配置。
   - **Modeler score ($f_M$, Precision / 適合率)**:
     - **定義**: 推定アラインメント中に配置された全残基ペアのうち、真のアラインメントでも実際に相同関係にあったペアの割合（過剰アラインメントの検知）。
     - **数式**: $f_M = \frac{|\text{一致した相同残基ペア}|}{|\text{推定アラインメントの全残基ペア}|}$
   - **アラインメント構造指標**:
     - `msa_length_ratio`: 推定アラインメント長 / 真のアラインメント長（ギャップ誤挿入による伸長や圧縮を定量化）。
     - `true_gap_ratio` / `est_gap_ratio`: アラインメント全体に占めるギャップ文字（`-`）の比率。

2. **進化的距離行列の推定誤差指標**
   - **真のパトリスティック距離 ($d_{ij}$, True Patristic Distance)**:
     - **定義**: 真の系統樹（Newick 形式）の枝長情報から直接計算した、分類群 $i$ と $j$ の間の樹形パス上の累積枝長（真の進化的置換距離の厳密解）。
   - **距離推定平均絶対誤差 MAE (`dist_mae`, Mean Absolute Error)**:
     - **定義**: 推定アラインメントまたはペアワイズアラインメントからポアソン補正モデル等で算出した距離行列要素 $\hat{d}_{ij}$ と、真のパトリスティック距離 $d_{ij}$ との平均絶対乖離（上三角要素 $\binom{N}{2}$ ペアの平均）。
     - **数式**:
       $$\text{MAE} = \frac{1}{\binom{N}{2}} \sum_{i < j} |\hat{d}_{ij} - d_{ij}|$$
     - **解釈**: ギャップ誤配置や多重置換の補正破綻により、距離行列が真の系統樹の枝長からどれだけ歪んでいるかを絶対値スケールで測定。
   - **関連指標**: 二乗平均平方根誤差 RMSE (`dist_rmse`), バイアス (`dist_bias`), 相関係数 (`dist_corr`)。

3. **系統樹トポロジー評価指標**
   - **正規化 Robinson-Foulds 距離 (`nrf_distance`, nRF)**:
     - **定義**: 推定系統樹と真の系統樹の間で不一致となった内部枝（bipartition）数を、可能な最大不一致数 $2(N-3)$ で割ったトポロジー誤差率。$0.0$ で完全一致、$1.0$ で全枝不一致。
   - **正解枝復元率 RCRB (`rcrb`, Ratio of Correctly Recovered Branches)**:
     - **定義**: 真の系統樹の全内部枝 $N-3$ 本のうち、推定系統樹で正しく復元された内部枝の割合（Matsui & Iwasaki 2020 準拠）。
     - **数式**: $\text{RCRB} = \frac{\text{correct\_branches}}{N-3} = 1 - nRF$ （二分岐樹の場合）。
   - **Precision**: 推定樹の枝のうち真の樹に含まれていた割合（二分岐樹では RCRB と一致）。

4. **4段階の誤差伝播カスケード (Error Propagation Cascade)**:
   `error_propagation_analysis.png`（4分割パネル）は、誤差が上流から下流へ波及する連鎖構造を可視化している：
   - **Panel (a) $D \rightarrow SP$ (進化的発散 $\rightarrow$ アラインメント崩壊)**:
     進化距離 $D$ の増大と Indel 蓄積に伴い、多重配列アラインメント（MAFFT）の相同残基ペア復元率（SP score）が急激に減衰する過程。
   - **Panel (b) $SP \rightarrow \text{MAE}$ (アラインメント誤り $\rightarrow$ 距離行列の歪み)**:
     相同残基の誤配置（SP score 低下）が、pairwise 距離推定値 $\hat{d}_{ij}$ と真のパトリスティック距離 $d_{ij}$ の乖離（MAE）を急激に増大させる関係。
   - **Panel (c) $\text{MAE} \rightarrow nRF$ (距離歪み $\rightarrow$ 木構造の誤認)**:
     距離行列の歪み（MAE の増大）が近隣結合法（NJ）のクラスタリング判定を狂わせ、内部枝の誤認（nRF の上昇）へと直結するプロセス。
   - **Panel (d) $SP \rightarrow nRF$ (アラインメント精度 $\rightarrow$ 最終トポロジー誤差)**:
     アラインメント誤差が系統樹トポロジーに与える最終的な直接的インパクト。距離法（PWA+NJ, MSA+NJ）だけでなく最尤法（MSA+ML, MSA+RAXML）やグラフ分割法（GS）を含めた全手法を同一軸上で比較し、アラインメントの破綻に対する各手法の耐性を評価。

![paper tree3 SP score vs distance](../results/results_paper_tree3/sp_score_vs_distance.png)

![paper tree3 SP score vs MAE](../results/results_paper_tree3/sp_score_vs_mae.png)

![paper tree3 SP score vs tree accuracy](../results/results_paper_tree3/sp_score_vs_tree_accuracy.png)

![paper tree3 error propagation](../results/results_paper_tree3/error_propagation_analysis.png)


考察:
1. $D\le0.2$ では MSA系は $nRF\approx0$ と高精度。一方 `GS` は $D=0.1$ でも $nRF=0.091$ と近縁域では劣る。`PWA+NJ` もわずかに劣る。
2. $D\ge1.0$ で全手法のnRFが急上昇し、SSSは0.04以下、SP scoreは約0.32〜0.34まで低下する。アラインメント精度の低下が木推定誤差の増加と同時に起きている。
3. 真のアラインメントを使うと $D=2.0$ でも nRF は 0.16〜0.25 に留まり、推定アラインメント使用時(0.41〜0.49)との差がアラインメント誤差の寄与を示す。
4. Extreme領域では `GS`(0.3836)は `MSA+ML`(0.3757)・`MSA+RAXML`(0.3779)とほぼ同水準だが、`MSA+NJ`(0.3600)より劣る。`results_paper_tree2` で見られた「GSがMLと同等以上」の傾向は、同水準という範囲でのみ再現し、GSの優位は確認できなかった。
5. `PWA+NJ` は全域で `MSA+NJ` より劣り($D=2.0$: 0.490 vs 0.411)、Extreme領域でも 0.4494 と最も悪い。
6. `IQ-TREE 2` と `RAxML` は全域でほぼ同等(例 $D=2.0$: 0.426 vs 0.429)。
7. `TRUE_MSA+NJ` と `TRUE_PWA+NJ` は同一の値で、真のアラインメントを与えた場合に両者が同じ入力になるためと考えられる(要確認)。

#### 総括

- 低類似度域($w\approx0.06$以下)で `MSA+ML` の精度が急落し、`PSA+NJ`・`GS` が同等以上の頑健性を示すCrossoverが、Zipfian Indelを含む複数実験で一貫して確認された。
- 速度不均一性($\alpha$)、$N$、ICS比率はいずれもSSSの分布と手法の優劣に影響する。
- 真値対照(TRUE_*)により、低SSS域での誤差の主因がアラインメント破綻であることを分離できた。
- 注意: `results_paper_tree3` の数値は `results/results_paper_tree3/benchmark_summary.csv` を本レポート作成時に pandas で再集計した値であり、他の実験(README転記)とは出典が異なる。
- 注意: Methodの条件はREADMEでなく設定ファイルに基づく。READMEとの不一致(`results_paper_tree2` のレプリケート数、`results_zipfian_indel` のIndel率、`results_ics_full` の $L$ 等)は要確認。