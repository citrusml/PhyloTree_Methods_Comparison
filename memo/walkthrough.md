# 30-Condition Evolutionary Simulation, MSA & Homology Benchmark

## 1. Overview & Experimental Matrix

This comprehensive benchmark tracks **30 unique evolutionary simulation conditions** across **5 evolutionary distances** ($D \in \{0.1, 0.5, 1.0, 2.0, 3.0\}$), **3 Indel rate regimes**, and **2 initial sequence lengths** ($L = 100\text{ aa}$ and $L = 500\text{ aa}$) with $N = 32\text{ taxa}$, using `LG+G` ($\alpha=1.0$) and a log-normal relaxed clock ($\sigma=0.5$).

### The 30-Condition Matrix
- **Initial Sequence Lengths**: $L = 100\text{ aa}$ (Short) and $L = 500\text{ aa}$ (Long)
- **Evolutionary Distances**: $D = 0.1$ (Ultra-shallow), $D = 0.5$ (Shallow), $D = 1.0$ (Moderate), $D = 2.0$ (Deep), $D = 3.0$ (Ultra-deep)
- **Indel Regimes**:
  1. **Conservative / Balanced (`indel0.01`)**: $\text{Ins} = 0.01, \text{Del} = 0.01$
  2. **Moderate / High Indel (`indel0.05`)**: $\text{Ins} = 0.05, \text{Del} = 0.05$
  3. **Asymmetric Expansion (`ins0.1_del0.05`)**: $\text{Ins} = 0.10, \text{Del} = 0.05$

---

## 2. Master Statistical Summary (All 30 Conditions)

### A. $L = 100\text{ aa}$ Benchmark Results

| Indel Regime | Distance $D$ | Mean Length (aa) | Std Dev | Min - Max | MSA Length (cols) | Total Gaps | Gap % | Expansion Ratio |
|:---|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| **Round 1** (`Ins=0.01, Del=0.01`) | **0.1** | 100.0 | $\pm 0.39$ | 98 - 101 | 101 | 33 | **1.02%** | 1.01x |
| | **0.5** | 101.0 | $\pm 2.14$ | 99 - 108 | 111 | 321 | **9.04%** | 1.11x |
| | **1.0** | 99.9 | $\pm 1.90$ | 95 - 105 | 116 | 514 | **13.85%** | 1.16x |
| | **2.0** | 99.4 | $\pm 21.42$ | 22 - 159 | 192 | 2,964 | **48.24%** | 1.92x |
| | **3.0** | 95.3 | $\pm 7.16$ | 75 - 108 | 125 | 950 | **23.75%** | 1.25x |
| **Round 2** (`Ins=0.05, Del=0.05`) | **0.1** | 100.9 | $\pm 2.15$ | 99 - 106 | 119 | 578 | **15.18%** | 1.19x |
| | **0.5** | 98.8 | $\pm 9.31$ | 58 - 119 | 158 | 1,896 | **37.50%** | 1.58x |
| | **1.0** | 87.8 | $\pm 22.78$ | 63 - 170 | 270 | 5,829 | **67.47%** | 2.70x |
| | **2.0** | 67.5 | $\pm 28.93$ | 20 - 140 | 185 | 3,759 | **63.50%** | 1.85x |
| | **3.0** | 101.1 | $\pm 38.52$ | 13 - 201 | 308 | 6,622 | **67.19%** | 3.08x |
| **Round 3** (`Ins=0.10, Del=0.05`) | **0.1** | 98.2 | $\pm 1.24$ | 95 - 100 | 105 | 217 | **6.46%** | 1.05x |
| | **0.5** | 107.9 | $\pm 20.89$ | 15 - 167 | 227 | 3,810 | **52.45%** | 2.27x |
| | **1.0** | 119.0 | $\pm 32.74$ | 55 - 188 | 347 | 7,296 | **65.71%** | 3.47x |
| | **2.0** | 202.3 | $\pm 68.72$ | 90 - 359 | 508 | 9,782 | **60.17%** | 5.08x |
| | **3.0** | 200.3 | $\pm 63.18$ | 80 - 330 | 597 | 12,694 | **66.45%** | 5.97x |

---

### B. $L = 500\text{ aa}$ Benchmark Results

| Indel Regime | Distance $D$ | Mean Length (aa) | Std Dev | Min - Max | MSA Length (cols) | Total Gaps | Gap % | Expansion Ratio |
|:---|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| **Round 1** (`Ins=0.01, Del=0.01`) | **0.1** | 500.1 | $\pm 1.67$ | 496 - 502 | 505 | 157 | **0.97%** | 1.01x |
| | **0.5** | 438.8 | $\pm 54.41$ | 387 - 575 | 598 | 5,096 | **26.63%** | 1.20x |
| | **1.0** | 518.6 | $\pm 24.79$ | 475 - 560 | 709 | 6,094 | **26.86%** | 1.42x |
| | **2.0** | 487.9 | $\pm 36.81$ | 428 - 613 | 776 | 9,220 | **37.13%** | 1.55x |
| | **3.0** | 535.0 | $\pm 57.14$ | 438 - 641 | 967 | 13,824 | **44.67%** | 1.93x |
| **Round 2** (`Ins=0.05, Del=0.05`) | **0.1** | 499.2 | $\pm 8.33$ | 467 - 514 | 588 | 2,842 | **15.10%** | 1.18x |
| | **0.5** | 468.7 | $\pm 61.43$ | 340 - 563 | 857 | 12,426 | **45.31%** | 1.71x |
| | **1.0** | 493.6 | $\pm 39.80$ | 408 - 575 | 1,077 | 18,668 | **54.17%** | 2.15x |
| | **2.0** | 520.9 | $\pm 130.70$ | 299 - 792 | 1,392 | 27,876 | **62.58%** | 2.78x |
| | **3.0** | 607.6 | $\pm 176.46$ | 166 - 851 | 1,534 | 29,646 | **60.39%** | 3.07x |
| **Round 3** (`Ins=0.10, Del=0.05`) | **0.1** | 514.3 | $\pm 25.15$ | 483 - 581 | 804 | 9,272 | **36.04%** | 1.61x |
| | **0.5** | 510.4 | $\pm 49.87$ | 408 - 618 | 1,278 | 24,563 | **60.06%** | 2.56x |
| | **1.0** | 651.4 | $\pm 92.78$ | 497 - 877 | 1,912 | 40,340 | **65.93%** | 3.82x |
| | **2.0** | 740.2 | $\pm 251.05$ | 275 - 1,142 | 2,638 | 60,731 | **71.94%** | 5.28x |
| | **3.0** | 1,457.1 | $\pm 209.28$ | 1,104 - 1,918 | 4,245 | 89,213 | **65.68%** | 8.49x |

---

## 3. Visualizations: 15-Condition Master MSA Raster Grids (3 $\times$ 5)

### A. $L = 100\text{ aa}$ Master Raster Overview (All 15 Conditions)
![L=100 All 15 Conditions MSA Raster Grid](/Users/kazukiaibara/.gemini/antigravity-ide/brain/3991a44c-0669-480a-b064-0ea1599a516f/figures/l100_all_15conditions_msa_grid.png)

### B. $L = 500\text{ aa}$ Master Raster Overview (All 15 Conditions)
![L=500 All 15 Conditions MSA Raster Grid](/Users/kazukiaibara/.gemini/antigravity-ide/brain/3991a44c-0669-480a-b064-0ea1599a516f/figures_L500/l500_all_15conditions_msa_grid.png)

---

## 4. Individual Regime 5-Panel Detailed Raster Plots

````carousel
![L=100 Ins=0.01, Del=0.01 5 Distances](/Users/kazukiaibara/.gemini/antigravity-ide/brain/3991a44c-0669-480a-b064-0ea1599a516f/figures/msa_raster_indel0.01.png)
<!-- slide -->
![L=100 Ins=0.05, Del=0.05 5 Distances](/Users/kazukiaibara/.gemini/antigravity-ide/brain/3991a44c-0669-480a-b064-0ea1599a516f/figures/msa_raster_indel0.05.png)
<!-- slide -->
![L=100 Ins=0.10, Del=0.05 5 Distances](/Users/kazukiaibara/.gemini/antigravity-ide/brain/3991a44c-0669-480a-b064-0ea1599a516f/figures/msa_raster_ins0.1_del0.05.png)
<!-- slide -->
![L=500 Ins=0.01, Del=0.01 5 Distances](/Users/kazukiaibara/.gemini/antigravity-ide/brain/3991a44c-0669-480a-b064-0ea1599a516f/figures_L500/msa_raster_indel0.01.png)
<!-- slide -->
![L=500 Ins=0.05, Del=0.05 5 Distances](/Users/kazukiaibara/.gemini/antigravity-ide/brain/3991a44c-0669-480a-b064-0ea1599a516f/figures_L500/msa_raster_indel0.05.png)
<!-- slide -->
![L=500 Ins=0.10, Del=0.05 5 Distances](/Users/kazukiaibara/.gemini/antigravity-ide/brain/3991a44c-0669-480a-b064-0ea1599a516f/figures_L500/msa_raster_ins0.1_del0.05.png)
````

---

## 5. Sequence Homology vs Evolutionary Distance ($L = 500\text{ aa}, \text{Ins} = 0.10, \text{Del} = 0.05$)

We computed **2,480 pairwise sequence alignments** across all 5 evolutionary distances using both Needleman-Wunsch Pairwise Alignment (PWA) and MAFFT (MSA), comparing them against the true phylogenetic tree patristic distance.

### A. Distance vs Sequence Identity & Similarity Decay
![Homology Decay vs True Tree Distance](/Users/kazukiaibara/.gemini/antigravity-ide/brain/3991a44c-0669-480a-b064-0ea1599a516f/figures_homology/identity_vs_true_distance_scatter.png)

- **$D = 0.1$**: Patristic Dist $0.05 \sim 0.18$, Identity **$85.4\% \sim 99.8\%$**, Similarity **$92.1\% \sim 100.0\%$**. Alignment is nearly gap-free.
- **$D = 0.5$**: Patristic Dist $0.20 \sim 0.85$, Identity **$58.2\% \sim 82.5\%$**, Similarity **$72.0\% \sim 90.4\%$**. Smooth exponential decay following theoretical curve $e^{-d}$.
- **$D = 1.0$**: Patristic Dist $0.40 \sim 1.80$, Identity **$42.0\% \sim 65.0\%$**, Similarity **$60.0\% \sim 78.0\%$**.
- **$D = 2.0 \sim 3.0$**: Patristic Dist $1.5 \sim 6.4$, Identity hits the **Twilight Zone** plateau at **$25\% \sim 35\%$** (due to DP forced alignment artifact), while BLOSUM62 similarity retains phylogenetic signal at **$45\% \sim 55\%$**.

### B. Pairwise Alignment (PWA) vs Multiple Sequence Alignment (MSA)
![PWA vs MSA Homology Comparison](/Users/kazukiaibara/.gemini/antigravity-ide/brain/3991a44c-0669-480a-b064-0ea1599a516f/figures_homology/pwa_vs_msa_homology_comparison.png)

### C. Distance Saturation & Correction
![Distance Saturation Curves](/Users/kazukiaibara/.gemini/antigravity-ide/brain/3991a44c-0669-480a-b064-0ea1599a516f/figures_homology/poisson_distance_estimation_accuracy.png)
- Standard Poisson distance caps at $\approx 1.3 \sim 1.5$ when $d_{\text{true}} > 2.0$.
- **Gamma-Poisson distance ($\alpha = 1.0$)** accurately recovers true tree distance up to $d > 4.5$ without saturation.
