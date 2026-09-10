nextflow.enable.dsl=2

/*
 * Unified Comparative Benchmark Pipeline: AliSim vs INDELible (research2)
 *
 * Evaluates:
 *   - INDELible (rate=0.10/0.10)
 *   - INDELible (rate=0.05/0.05)
 *   - AliSim (rate=0.10/0.10)
 *   - AliSim (rate=0.05/0.05)
 * Across:
 *   - Lengths: [500, 1000, 1500]
 *   - Distances: [0.1, 0.2, 0.5, 1.0, 1.5, 2.0]
 *   - Replicates: 10
 * Under:
 *   - Taxa: 20
 *   - Tree model: paper_yule (backward Yule + log branch distribution)
 *   - Indel size: POW 1.7 50 (Zipfian power-law)
 *   - Substitution model: LG + Gamma4 (alpha=1.0)
 */

include { SIMULATE_COMPARATIVE_DATA } from './modules/simulate_comparative_data'
include {
    RUN_PWA_NJ;
    RUN_MAFFT;
    RUN_MSA_NJ;
    RUN_MSA_ML;
    RUN_TRUE_MSA_NJ;
    RUN_TRUE_MSA_ML;
    RUN_SSS
} from './modules/reconstruct_trees'
include { CALCULATE_PROPERTIES } from './modules/calculate_properties'
include { COLLECT_AND_PLOT_RESEARCH2 } from './modules/collect_and_plot'

workflow {
    ch_conditions = Channel.fromList([
        ['indelible_rate0.10', 'indelible', 0.10, 0.10],
        ['indelible_rate0.05', 'indelible', 0.05, 0.05],
        ['alisim_rate0.10',    'alisim',    0.10, 0.10],
        ['alisim_rate0.05',    'alisim',    0.05, 0.05]
    ])

    def dist_list = (params.distances instanceof Collection) ? params.distances.flatten() : [params.distances]
    def len_list  = (params.lengths instanceof Collection) ? params.lengths.flatten() : [params.lengths]
    def rep_count = (params.replicates instanceof List) ? params.replicates[0] as int : params.replicates as int

    ch_distances = Channel.fromList(dist_list)
    ch_lengths   = Channel.fromList(len_list)
    ch_reps      = Channel.fromList(1..rep_count)

    ch_sim_input = ch_conditions
        .combine(ch_distances)
        .combine(ch_lengths)
        .combine(ch_reps)

    // 1. シミュレーション実行
    SIMULATE_COMPARATIVE_DATA(ch_sim_input)
    ch_sim = SIMULATE_COMPARATIVE_DATA.out.sim_data

    // 2. 入力チャネルの準備
    // ch_sim structure: [key(0), cond(1), sim(2), ins(3), del(4), dist(5), len(6), rep(7), tree(8), seqs(9), true_msa(10)]
    ch_pwa_input     = ch_sim.map { it -> [it[0], it[9]] }
    ch_mafft_input   = ch_sim.map { it -> [it[0], it[9]] }
    ch_true_nj_input = ch_sim.map { it -> [it[0], it[10]] }
    ch_true_ml_input = ch_sim.map { it -> [it[0], it[10]] }
    ch_sss_input     = ch_sim.map { it -> [it[0], it[5], it[6], it[7], it[9]] }

    // 3. 各推定・解析プロセスの実行
    RUN_PWA_NJ(ch_pwa_input)
    RUN_MAFFT(ch_mafft_input)
    RUN_MSA_NJ(RUN_MAFFT.out.msa)
    RUN_MSA_ML(RUN_MAFFT.out.msa)
    RUN_TRUE_MSA_NJ(ch_true_nj_input)
    RUN_TRUE_MSA_ML(ch_true_ml_input)
    RUN_SSS(ch_sss_input)

    // 4. キーによる統合
    ch_joined = ch_sim
        .join(RUN_MAFFT.out.msa)
        .join(RUN_PWA_NJ.out.tree)
        .join(RUN_MSA_NJ.out.tree)
        .join(RUN_MSA_ML.out.tree)
        .join(RUN_TRUE_MSA_NJ.out.tree)
        .join(RUN_TRUE_MSA_ML.out.tree)
        .join(RUN_SSS.out.csv)

    // 5. 各レプリケートの MSA 性質および系統樹精度を算出
    CALCULATE_PROPERTIES(ch_joined)

    // 6. 全結果の集約・プロット
    COLLECT_AND_PLOT_RESEARCH2(CALCULATE_PROPERTIES.out.csv.collect())
}
