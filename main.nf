nextflow.enable.dsl=2

/*
 * Unified Phylogenetic Benchmark Pipeline (main.nf)
 *
 * 全実験共通のコアワークフロー（プロジェクトルート用エントリポイント）。
 * 各プロセスは `modules/` 配下の独立したモジュールファイルからインポートされます。
 * 実験ごとの条件（Indelモデル、置換モデル、ギャップペナルティ、距離・配列長・反復数など）は、
 * すべて `next_configs/*.config` の `params` ブロックで指定します。
 *
 * Execution:
 *   nextflow run main.nf -c next_configs/your_experiment.config -profile supercomputer
 */

// モジュールのインポート
include { SIMULATE_DATA }    from './modules/simulate_data'
include { RUN_PWA_NJ }       from './modules/run_pwa_nj'
include { RUN_MAFFT }        from './modules/run_mafft'
include { RUN_MSA_NJ }       from './modules/run_msa_nj'
include { RUN_MSA_ML }       from './modules/run_msa_ml'
include { RUN_TRUE_PWA_NJ }  from './modules/run_true_pwa_nj'
include { RUN_TRUE_MSA_NJ }  from './modules/run_true_msa_nj'
include { RUN_TRUE_MSA_ML }  from './modules/run_true_msa_ml'
include { COLLECT_AND_PLOT } from './modules/collect_and_plot'

workflow {
    def dist_list = (params.distances instanceof Collection) ? params.distances.flatten() : [params.distances]
    def len_list  = (params.lengths instanceof Collection) ? params.lengths.flatten() : [params.lengths]
    ch_distances = Channel.fromList(dist_list)
    ch_lengths   = Channel.fromList(len_list)

    def reps = (params.replicates instanceof List) ? params.replicates[0] as int : params.replicates as int
    def chk_size = (params.chunk_size instanceof List) ? params.chunk_size[0] as int : params.chunk_size as int
    def num_chunks = Math.max(1, Math.ceil(reps / chk_size) as int)
    ch_chunks = Channel.fromList( (1..num_chunks).collect { c ->
        def r_start = (c - 1) * chk_size + 1
        def r_end   = Math.min(reps, c * chk_size)
        [c, r_start, r_end]
    } )

    ch_params = ch_distances.combine(ch_lengths).combine(ch_chunks)

    SIMULATE_DATA(ch_params)

    ch_sim_data = SIMULATE_DATA.out.sim_data

    // 1. PWA+NJ pipeline (Gap Open 20, Extension 2)
    RUN_PWA_NJ(ch_sim_data)

    // 2. Shared MAFFT calculation
    RUN_MAFFT(ch_sim_data)
    ch_msa_data = RUN_MAFFT.out.msa_data

    // 3. MSA+NJ and MSA+ML pipelines
    RUN_MSA_NJ(ch_msa_data)
    RUN_MSA_ML(ch_msa_data)

    ch_all_csvs = RUN_PWA_NJ.out.csv.mix(RUN_MSA_NJ.out.csv).mix(RUN_MSA_ML.out.csv)

    // 4. True PWA + NJ
    if (params.run_true_pwa) {
        ch_true_msa_data = SIMULATE_DATA.out.true_msa_data
        RUN_TRUE_PWA_NJ(ch_true_msa_data)
        ch_all_csvs = ch_all_csvs.mix(RUN_TRUE_PWA_NJ.out.csv)
    }

    // 5. True MSA + NJ / True MSA + ML
    if (params.run_true_msa) {
        ch_true_msa_data = SIMULATE_DATA.out.true_msa_data
        RUN_TRUE_MSA_NJ(ch_true_msa_data)
        RUN_TRUE_MSA_ML(ch_true_msa_data)
        ch_all_csvs = ch_all_csvs.mix(RUN_TRUE_MSA_NJ.out.csv).mix(RUN_TRUE_MSA_ML.out.csv)
    }

    COLLECT_AND_PLOT(ch_all_csvs.collect())
}
