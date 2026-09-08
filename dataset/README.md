# Dataset: HOG Ortholog Groups & Similarity Analysis

本ディレクトリには、Matsui & Iwasaki (Syst. Biol., 2020) の比較検証に用いるオーソログデータセット（HOGs）、生物種リスト（Taxa）、およびそれらに対する配列類似度スコア（SSS: Sequence Similarity Score）と一致率（Identity %）の解析スクリプト・結果が格納されています。

---

## ディレクトリ構成

```text
dataset/
├── act_hogs.tsv          # Actinobacteria 100 HOGs
├── act_taxa.tsv          # Actinobacteria 50 taxa
├── amn_hogs.tsv          # Amniota 100 HOGs
├── amn_taxa.tsv          # Amniota 10 taxa
├── bac_hogs.tsv          # Bacteria 100 HOGs
├── bac_taxa.tsv          # Bacteria 50 taxa
├── euk_hogs.tsv          # Eukaryota 100 HOGs
├── euk_taxa.tsv          # Eukaryota 50 taxa
├── scripts/
│   ├── similarity_calculator.py       # Needleman-Wunsch BLOSUM62 による SSS / Identity 計算コア
│   ├── analyze_dataset_similarity.py  # 全対全ペアワイズ並列バッチ解析スクリプト
│   └── plot_dataset_similarity.py     # 論文基準（300 DPI）可視化プロット生成スクリプト
└── results/
    ├── dataset_level_summary.csv      # 4生物群のマクロ統計サマリー
    ├── hog_similarity_summary.csv     # 全400 HOGの個別統計
    ├── pairwise_similarity_details.csv# 全124,426配列ペアの詳細テーブル
    ├── dataset_similarity_multiplot.png # SSS & Identity 分布バイオリンプロット（4群統合）
    ├── sss_distribution_comparison.png  # SSS 分布プロット
    ├── identity_distribution_comparison.png # 一致率分布プロット
    └── sss_vs_identity_correlation.png  # SSS vs 一致率 相関散布図
```

---

## 実行方法

リポジトリルートから以下のコマンドで解析・プロットを再実行できます：

```bash
# 1. 類似度（SSS）および一致率（Identity）の並列計算（400 HOGs / 6 workers）
micromamba run -n phylomethod_env python dataset/scripts/analyze_dataset_similarity.py --workers 6

# 2. グラフ描画（バイオリンプロット・相関散布図）
micromamba run -n phylomethod_env python dataset/scripts/plot_dataset_similarity.py
```

---

## 解析結果サマリー

| データセット | 分類群 | 解析 HOG 数 | 計算ペア数 | 平均 SSS ($w$) | 中央値 SSS | 平均 Identity (%) | 中央値 Identity (%) | SSS $\le 0.10$ の割合 |
| :--- | :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **`amn`** | **Amniota** (10種) | 100 | 4,500 | **0.7371** | **0.7688** | **79.55%** | **83.37%** | **0% (0 / 100)** |
| **`euk`** | **Eukaryota** (50種) | 100 | 38,292 | **0.5448** | **0.5362** | **63.43%** | **60.57%** | **0% (0 / 100)** |
| **`act`** | **Actinobacteria** (50種) | 100 | 41,048 | **0.4263** | **0.4044** | **53.26%** | **48.34%** | **1% (1 / 100)** |
| **`bac`** | **Bacteria** (50種) | 100 | 40,586 | **0.3070** | **0.2406** | **42.60%** | **36.89%** | **8% (8 / 100)** |
| **全体** | **全4群統合** | **400** | **124,426** | **0.5038** | **0.4875** | **59.71%** | **57.29%** | **2.25% (9 / 400)** |
