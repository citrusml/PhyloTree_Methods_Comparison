# PhyloTree Methods Comparison (Regime Map Benchmark)

[![Nextflow](https://img.shields.io/badge/Nextflow-DSL2-brightgreen)](https://www.nextflow.io/)
[![Environment](https://img.shields.io/badge/Package%20Manager-Micromamba-blue)](https://mamba.readthedocs.io/en/latest/user_guide/micromamba.html)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

A high-performance, automated benchmark pipeline to systematically map the **Regime Map (相図)** of phylogenetic reconstruction methods across evolutionary distance ($D$) and sequence length ($L$).

This study investigates the performance trade-off and regime boundary ($\Delta RF = 0$) between:
1. **`PWA + NJ`**: Pairwise Sequence Alignment (Needleman-Wunsch with LG matrix) $\rightarrow$ Poisson/Gamma Distance $\rightarrow$ Neighbor-Joining (RapidNJ / FastME).
2. **`MSA + NJ` (Control)**: Multiple Sequence Alignment (MAFFT) $\rightarrow$ Identical Poisson/Gamma Distance $\rightarrow$ Neighbor-Joining (isolating MSA noise from tree inference).
3. **`MSA + ML`**: Multiple Sequence Alignment (MAFFT) $\rightarrow$ Maximum Likelihood with automated model selection (IQ-TREE 2 ModelFinder).

---

## 2x2 Factorial Experimental Design

```text
                  Tree Inference Algorithm
                 Neighbor-Joining (NJ)    Maximum Likelihood (ML)
Alignment
  MSA (MAFFT)          MSA + NJ                  MSA + ML
  PWA (Pairwise)       PWA + NJ                     -
```

---

## Directory Structure

```text
PhyloTree_Methods_Comparison/
├── main.nf                      # Nextflow pipeline orchestration workflow (DSL2)
├── nextflow.config              # Cluster resource configuration (LSF & Local)
├── run_nextflow_job.sh          # LSF batch submission wrapper script
├── environment.yml              # Micromamba / Conda environment definition
├── bin/
│   ├── simulate_data.py         # Tree generation (relaxed clock) & AliSim (LG+G) sequence simulation
│   ├── run_pwa_nj.py            # Needleman-Wunsch (LG matrix) + Poisson distance + RapidNJ
│   ├── run_msa_nj.py            # MAFFT MSA + Poisson distance + RapidNJ (Control group)
│   ├── run_msa_ml.py            # MAFFT MSA + IQ-TREE 2 ModelFinder (-m MFP -B 1000)
│   ├── evaluate_trees.py        # DendroPy unrooted RF / nRF metric computation
│   ├── plot_regime_map.py       # 2D Regime Map Heatmap generation (Seaborn)
│   ├── run_local_test.py        # Local Phase 0 mini benchmark runner
│   └── aggregate_work_results.py# Fast recovery & aggregation script from Nextflow work directory
├── results/                     # Output directory for benchmark artifacts
│   ├── benchmark_summary.csv    # Aggregated performance metric table
│   └── regime_map_delta_nrf.png # 2D phase diagram heatmap (ΔnRF)
└── plan/
    └── design2.md               # Detailed research plan and parameter definitions
```

---

## Installation & Setup

### 1. Clone the Repository

```bash
git clone https://github.com/citrusml/PhyloTree_Methods_Comparison.git
cd PhyloTree_Methods_Comparison
```

### 2. Install Micromamba (if not already installed)

```bash
curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj bin/micromamba
mkdir -p ~/bin
mv bin/micromamba ~/bin/micromamba
export PATH="$HOME/bin:$PATH"

# Add to ~/.bashrc for persistence
echo 'export PATH="$HOME/bin:$PATH"' >> ~/.bashrc
source ~/.bashrc
```

### 3. Create the Environment

```bash
micromamba create -f environment.yml -y
micromamba activate phylomethod_env
```

---

## How to Run

### Option A: Supercomputer (LSF Cluster Batch Job - Recommended)

To run the full-scale benchmark (2,500 parameter combinations $\times$ replicates) stably on an LSF supercomputer without being affected by SSH disconnects:

```bash
# Submit the Nextflow master job to LSF
bsub < run_nextflow_job.sh
```

**Monitor job progress:**
```bash
# Check running LSF jobs
bjobs

# Follow real-time execution log
tail -f logs/nextflow_lsf_*.log
```

---

### Option B: Local / Small-scale Verification Test (Phase 0)

To quickly test the end-to-end pipeline locally (or on a login node) with a mini parameter grid:

```bash
micromamba activate phylomethod_env

# Run via Nextflow local profile
nextflow run main.nf -profile local
```

Or run the standalone Python test runner:
```bash
python bin/run_local_test.py
```

---

### Option C: Direct CLI Execution with Resume Support

```bash
micromamba activate phylomethod_env
mkdir -p logs

# Run in background with input detached
nohup nextflow run main.nf -profile supercomputer -resume -ansi-log false > logs/run_nextflow.log 2>&1 < /dev/null &

# Monitor logs
tail -f logs/run_nextflow.log
```

---

## Result Aggregation & Instant Recovery

If a pipeline run was interrupted or you want to immediately aggregate all completed evaluation records directly from the `work/` directory without re-running any computations:

```bash
python bin/aggregate_work_results.py
```

This will automatically search all `work/` subdirectories, deduplicate results, and produce:
- `results/benchmark_summary.csv`
- `results/regime_map_delta_nrf.png`

---

## Output Metrics & Interpretation

| Output File | Description |
|---|---|
| `results/benchmark_summary.csv` | Full table of distance, length, replicate, pipeline, RF distance, Normalized RF ($nRF$), and selected evolutionary models. |
| `results/regime_map_delta_nrf.png` | 2D Phase Diagram displaying $\Delta nRF = nRF_{\text{PWA+NJ}} - nRF_{\text{MSA+ML}}$. |

### Interpreting the Regime Map ($\Delta nRF$)
- **Blue ($\Delta nRF < 0$)**: **`PWA+NJ` outperforms `MSA+ML`** (typically in short sequences / large divergence where MSA noise degrades ML).
- **White ($\Delta nRF \approx 0$)**: **Regime Boundary** (the critical transition curve where both methods perform equally).
- **Red ($\Delta nRF > 0$)**: **`MSA+ML` outperforms `PWA+NJ`** (typically in long alignments / close relatives where ML statistical power dominates).

---

## License

This project is licensed under the MIT License.