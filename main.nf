nextflow.enable.dsl=2

/*
 * Unified Phylogenetic Benchmark Pipeline (main.nf)
 *
 * 全実験共通のコアワークフロー（プロジェクトルート用エントリポイント）。
 * 各プロセスは `modules/` 配下の独立したモジュールファイルからインポートされます。
 * 実験ごとの条件（Indelモデル、置換モデル、ギャップペナルティ、距離・配列長・反復数など）は、
 * すべて `next_configs/*.config` の `params` ブロックで指定します。
 *
 * 各手法のモジュール実行は config または実行時引数のスイッチ（run_*）で個別に ON / OFF 可能です。
 * ただし Nextflow DSL2 の制約上、プロセス呼び出しは if ブロックの1段目のみ有効です。
 * 個別スイッチは Channel.filter() によって実装しています。
 *
 * Pipelines supported:
 *   - PWA+NJ      (params.run_pwa_nj)
 *   - MSA+NJ      (params.run_msa_nj)
 *   - MSA+ML      (params.run_msa_ml)
 *   - MSA+BI      (params.run_msa_bi / params.run_bi)
 *   - GS          (params.run_gs)
 *   - TRUE_PWA+NJ (params.run_true_pwa)
 *   - TRUE_MSA+NJ (params.run_true_msa / params.run_true_msa_nj)
 *   - TRUE_MSA+ML (params.run_true_msa / params.run_true_msa_ml)
 *   - TRUE_MSA+BI (params.run_true_bi / params.run_true_msa_bi)
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
include { RUN_MSA_RAXML }    from './modules/run_msa_raxml'
include { RUN_MSA_BI }       from './modules/run_msa_bi'
include { RUN_GS }           from './modules/run_gs'
include { RUN_TRUE_PWA_NJ }  from './modules/run_true_pwa_nj'
include { RUN_TRUE_MSA_NJ }  from './modules/run_true_msa_nj'
include { RUN_TRUE_MSA_ML }  from './modules/run_true_msa_ml'
include { RUN_TRUE_MSA_RAXML } from './modules/run_true_msa_raxml'
include { RUN_TRUE_MSA_BI }  from './modules/run_true_msa_bi'
include { CALCULATE_SSS }    from './modules/calculate_sss'
include { COLLECT_AND_PLOT } from './modules/collect_and_plot'

// 型安全なブール値パース関数（CLI文字列 "false"/"true" と Boolean 型の双方に対応）
def asBool(val, default_val = false) {
    if (val == null) return default_val
    if (val instanceof Boolean) return val
    def s = val.toString().trim().toLowerCase()
    return (s == 'true' || s == '1' || s == 'yes') ? true
         : (s == 'false' || s == '0' || s == 'no') ? false
         : default_val
}

workflow {
    def dist_list = (params.distances instanceof Collection) ? params.distances.flatten() : [params.distances]
    def len_list  = (params.lengths instanceof Collection) ? params.lengths.flatten() : [params.lengths]
    ch_distances = Channel.fromList(dist_list)
    ch_lengths   = Channel.fromList(len_list)

    def reps     = (params.replicates instanceof List) ? params.replicates[0] as int : params.replicates as int
    def chk_size = (params.chunk_size  instanceof List) ? params.chunk_size[0]  as int : params.chunk_size  as int
    def num_chunks = Math.max(1, Math.ceil(reps / chk_size) as int)
    ch_chunks = Channel.fromList( (1..num_chunks).collect { c ->
        [c, (c - 1) * chk_size + 1, Math.min(reps, c * chk_size)]
    } )

    ch_params = ch_distances.combine(ch_lengths).combine(ch_chunks)

    SIMULATE_DATA(ch_params)
    ch_sim_data      = SIMULATE_DATA.out.sim_data
    ch_true_msa_data = SIMULATE_DATA.out.true_msa_data   // 常に取得（filter で制御）

    // --- フラグ評価（全スイッチを事前に確定） ---
    def do_pwa_nj      = asBool(params.containsKey('run_pwa_nj') ? params.run_pwa_nj : true, true)
    def do_msa_nj      = asBool(params.containsKey('run_msa_nj') ? params.run_msa_nj : true, true)
    def do_msa_ml      = asBool(params.containsKey('run_msa_ml') ? params.run_msa_ml : true, true)
    def do_msa_raxml   = asBool(params.containsKey('run_msa_raxml') ? params.run_msa_raxml : false, false)
    // run_msa_bi と run_bi（エイリアス）のどちらかが true なら有効（OR ロジック）
    def do_msa_bi      = asBool(params.containsKey('run_msa_bi') ? params.run_msa_bi : false, false) ||
                         asBool(params.containsKey('run_bi')     ? params.run_bi     : false, false)
    def do_gs          = asBool(params.containsKey('run_gs') ? params.run_gs : true, true)
    def do_true_pwa    = asBool(params.containsKey('run_true_pwa') ? params.run_true_pwa : true, true)
    def master_tmsa    = asBool(params.containsKey('run_true_msa') ? params.run_true_msa : true, true)
    def do_true_msa_nj = master_tmsa && asBool(params.containsKey('run_true_msa_nj') ? params.run_true_msa_nj : true, true)
    def do_true_msa_ml = master_tmsa && asBool(params.containsKey('run_true_msa_ml') ? params.run_true_msa_ml : true, true)
    def do_true_msa_raxml = master_tmsa && asBool(params.containsKey('run_true_msa_raxml') ? params.run_true_msa_raxml : false, false)
    // run_true_msa_bi と run_true_bi（エイリアス）のどちらかが true なら有効（OR ロジック）
    def do_true_msa_bi = asBool(params.containsKey('run_true_msa_bi') ? params.run_true_msa_bi : false, false) ||
                         asBool(params.containsKey('run_true_bi')     ? params.run_true_bi     : false, false)
    def do_calc_sss    = asBool(params.containsKey('calc_sss') ? params.calc_sss : true, true)

    ch_all_csvs = Channel.empty()

    // 0. SSS (Sequence Similarity Score) 計算 (Matsui & Iwasaki 2020 準拠)
    if (do_calc_sss) {
        CALCULATE_SSS(ch_sim_data)
        ch_all_csvs = ch_all_csvs.mix(CALCULATE_SSS.out.csv)
    }

    // 1. PWA+NJ
    if (do_pwa_nj) {
        RUN_PWA_NJ(ch_sim_data)
        ch_all_csvs = ch_all_csvs.mix(RUN_PWA_NJ.out.csv)
    }

    // 2. MAFFT + MSA系（需要があるときのみ MAFFT を実行）
    def need_mafft = do_msa_nj || do_msa_ml || do_msa_raxml || do_msa_bi
    if (need_mafft) {
        RUN_MAFFT(ch_sim_data)
        ch_msa_data = RUN_MAFFT.out.msa_data
    }

    // 3. MSA+NJ, MSA+ML, MSA+RAXML, MSA+BI（各々独立した if、ネストなし）
    if (need_mafft && do_msa_nj) {
        RUN_MSA_NJ(RUN_MAFFT.out.msa_data)
        ch_all_csvs = ch_all_csvs.mix(RUN_MSA_NJ.out.csv)
    }
    if (need_mafft && do_msa_ml) {
        RUN_MSA_ML(RUN_MAFFT.out.msa_data)
        ch_all_csvs = ch_all_csvs.mix(RUN_MSA_ML.out.csv)
    }
    if (need_mafft && do_msa_raxml) {
        RUN_MSA_RAXML(RUN_MAFFT.out.msa_data)
        ch_all_csvs = ch_all_csvs.mix(RUN_MSA_RAXML.out.csv)
    }
    if (need_mafft && do_msa_bi) {
        RUN_MSA_BI(RUN_MAFFT.out.msa_data)
        ch_all_csvs = ch_all_csvs.mix(RUN_MSA_BI.out.csv)
    }

    // 4. GS（アライメントフリー）
    if (do_gs) {
        RUN_GS(ch_sim_data)
        ch_all_csvs = ch_all_csvs.mix(RUN_GS.out.csv)
    }

    // 5-9. 正解アライメント系（各々独立した if、ネストなし）
    if (do_true_pwa) {
        RUN_TRUE_PWA_NJ(ch_true_msa_data)
        ch_all_csvs = ch_all_csvs.mix(RUN_TRUE_PWA_NJ.out.csv)
    }
    if (do_true_msa_nj) {
        RUN_TRUE_MSA_NJ(ch_true_msa_data)
        ch_all_csvs = ch_all_csvs.mix(RUN_TRUE_MSA_NJ.out.csv)
    }
    if (do_true_msa_ml) {
        RUN_TRUE_MSA_ML(ch_true_msa_data)
        ch_all_csvs = ch_all_csvs.mix(RUN_TRUE_MSA_ML.out.csv)
    }
    if (do_true_msa_raxml) {
        RUN_TRUE_MSA_RAXML(ch_true_msa_data)
        ch_all_csvs = ch_all_csvs.mix(RUN_TRUE_MSA_RAXML.out.csv)
    }
    if (do_true_msa_bi) {
        RUN_TRUE_MSA_BI(ch_true_msa_data)
        ch_all_csvs = ch_all_csvs.mix(RUN_TRUE_MSA_BI.out.csv)
    }

    COLLECT_AND_PLOT(ch_all_csvs.collect())
}
