nextflow.enable.dsl=2

/*
 * FastME Options Benchmark Pipeline (Experiment 10)
 * Evaluates:
 * 1. MSA + FastME_LG_G (MAFFT -> FastME with LG+Gamma model + SPR + TI correction)
 * 2. PWA + FastME_SPR  (PWA Poisson Matrix -> FastME with SPR + TI correction)
 * across Distance D x Length L.
 * Chunked / Batched Execution: Groups 10 replicates into a single task.
 */

params.taxa        = 32
params.distances   = [0.1, 0.5, 1.0, 2.0, 3.0]
params.lengths     = [100, 300, 500, 1000, 1500]
params.replicates  = 100
params.chunk_size  = 10
params.birth_rate  = 0.1
params.death_rate  = 0.05
params.insert_rate = 0.05
params.delete_rate = 0.10
params.model       = "LG+G4"
params.alpha       = 1.0
params.dist_model  = "poisson"
params.gap_open    = 10.0
params.gap_extend  = 0.5
params.outdir      = "results/results_fastme_options"

process SIMULATE_DATA {
    tag "D=${dist}_L=${len}_chk=${chunk_id}[${rep_start}..${rep_end}]"

    input:
    tuple val(dist), val(len), val(chunk_id), val(rep_start), val(rep_end)

    output:
    tuple val(dist), val(len), val(chunk_id), val(rep_start), val(rep_end), path("true_tree_*.nwk"), path("seqs_*.fasta"), emit: sim_data

    script:
    """
    for rep in \$(seq ${rep_start} ${rep_end}); do
        python3 ${projectDir}/../bin/simulate_data.py \\
            --num_taxa ${params.taxa} \\
            --distance ${dist} \\
            --length ${len} \\
            --birth_rate ${params.birth_rate} \\
            --death_rate ${params.death_rate} \\
            --insert_rate ${params.insert_rate} \\
            --delete_rate ${params.delete_rate} \\
            --model "${params.model}" \\
            --alpha ${params.alpha} \\
            --outtree true_tree_\${rep}.nwk \\
            --outfasta seqs_\${rep}.fasta \\
            --seed \${rep}
    done
    """
}

process RUN_PWA_FASTME_SPR {
    tag "D=${dist}_L=${len}_chk=${chunk_id}[${rep_start}..${rep_end}]"

    input:
    tuple val(dist), val(len), val(chunk_id), val(rep_start), val(rep_end), path(true_trees), path(fastas)

    output:
    path("chunk_pwa_fastme_spr_D${dist}_L${len}_chk${chunk_id}.csv"), emit: csv
    path("pwa_fastme_spr_*.nwk")
    path("pwa_matrix_*.phylip")

    script:
    """
    for rep in \$(seq ${rep_start} ${rep_end}); do
        # 1. PWA によるペアワイズ距離行列の計算
        python3 ${projectDir}/../bin/run_pwa_nj.py \\
            --fasta seqs_\${rep}.fasta \\
            --outmatrix pwa_matrix_\${rep}.phylip \\
            --outtree temp_\${rep}.nwk \\
            --gap_open ${params.gap_open} \\
            --gap_extend ${params.gap_extend} \\
            --dist_model ${params.dist_model} \\
            --threads ${task.cpus}

        # 2. FastME を直接実行 (Balanced Minimum Evolution + SPR探索 + 三角不等式補正)
        fastme -i pwa_matrix_\${rep}.phylip -o pwa_fastme_spr_\${rep}.nwk -mB -s -q
        rm -f temp_\${rep}.nwk

        # 3. 系統樹精度の評価
        python3 ${projectDir}/../bin/evaluate_trees.py \\
            --truetree true_tree_\${rep}.nwk \\
            --esttree pwa_fastme_spr_\${rep}.nwk \\
            --pipeline PWA+FastME_SPR \\
            --distance ${dist} \\
            --length ${len} \\
            --alpha ${params.alpha} \\
            --replicate \${rep} \\
            --outcsv chunk_pwa_fastme_spr_D${dist}_L${len}_chk${chunk_id}.csv
    done
    """
}

process RUN_MAFFT {
    tag "D=${dist}_L=${len}_chk=${chunk_id}[${rep_start}..${rep_end}]"

    input:
    tuple val(dist), val(len), val(chunk_id), val(rep_start), val(rep_end), path(true_trees), path(fastas)

    output:
    tuple val(dist), val(len), val(chunk_id), val(rep_start), val(rep_end), path(true_trees), path("msa_*.fasta"), emit: msa_data

    script:
    """
    for rep in \$(seq ${rep_start} ${rep_end}); do
        python3 ${projectDir}/../bin/run_msa.py \\
            --fasta seqs_\${rep}.fasta \\
            --outmsa msa_\${rep}.fasta \\
            --threads ${task.cpus}
    done
    """
}

process RUN_MSA_FASTME_LG_G {
    tag "D=${dist}_L=${len}_chk=${chunk_id}[${rep_start}..${rep_end}]"

    input:
    tuple val(dist), val(len), val(chunk_id), val(rep_start), val(rep_end), path(true_trees), path(msas)

    output:
    path("chunk_msa_fastme_lgg_D${dist}_L${len}_chk${chunk_id}.csv"), emit: csv
    path("msa_fastme_lgg_*.nwk")

    script:
    """
    for rep in \$(seq ${rep_start} ${rep_end}); do
        # 1. FASTA形式 MSA を PHYLIP形式に変換
        python3 -c "from Bio import SeqIO; SeqIO.write(SeqIO.parse('msa_\${rep}.fasta', 'fasta'), 'msa_\${rep}.phy', 'phylip-relaxed')"

        # 2. FastME を直接実行 (LG置換モデル + ガンマalpha + SPR探索 + 三角不等式補正)
        fastme -i msa_\${rep}.phy -o msa_fastme_lgg_\${rep}.nwk -pL -g${params.alpha} -s -q
        rm -f msa_\${rep}.phy

        # 3. 系統樹精度の評価
        python3 ${projectDir}/../bin/evaluate_trees.py \\
            --truetree true_tree_\${rep}.nwk \\
            --esttree msa_fastme_lgg_\${rep}.nwk \\
            --pipeline MSA+FastME_LG_G \\
            --distance ${dist} \\
            --length ${len} \\
            --alpha ${params.alpha} \\
            --replicate \${rep} \\
            --outcsv chunk_msa_fastme_lgg_D${dist}_L${len}_chk${chunk_id}.csv
    done
    """
}

process COLLECT_AND_PLOT {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path csv_files

    output:
    path("benchmark_fastme_options_summary.csv")
    path("scaling_curves_fastme_options.png"), optional: true
    path("scaling_curves_fastme_options.pdf"), optional: true
    path("regime_map_fastme_options.png"), optional: true
    path("regime_map_fastme_options.pdf"), optional: true
    path("benchmark_fastme_options_statistics.csv"), optional: true

    script:
    """
    cat *.csv | awk 'NR==1 || \$0 !~ /^alpha/' > benchmark_fastme_options_summary.csv
    python3 ${projectDir}/../bin/plot_fastme_options.py --csv benchmark_fastme_options_summary.csv --outdir .
    """
}

workflow {
    ch_distances = Channel.fromList(params.distances)
    ch_lengths   = Channel.fromList(params.lengths)

    def num_chunks = Math.max(1, (params.replicates / params.chunk_size) as int)
    ch_chunks = Channel.fromList( (1..num_chunks).collect { c ->
        def r_start = (c - 1) * params.chunk_size + 1
        def r_end   = Math.min(params.replicates, c * params.chunk_size)
        [c, r_start, r_end]
    } )

    ch_params = ch_distances.combine(ch_lengths).combine(ch_chunks)

    SIMULATE_DATA(ch_params)
    ch_sim_data = SIMULATE_DATA.out.sim_data

    // 1. PWA + FastME_SPR
    RUN_PWA_FASTME_SPR(ch_sim_data)

    // 2. Shared MAFFT MSA
    RUN_MAFFT(ch_sim_data)
    ch_msa_data = RUN_MAFFT.out.msa_data

    // 3. MSA + FastME_LG_G
    RUN_MSA_FASTME_LG_G(ch_msa_data)

    ch_all_csvs = RUN_PWA_FASTME_SPR.out.csv
        .mix(RUN_MSA_FASTME_LG_G.out.csv)
        .collect()

    COLLECT_AND_PLOT(ch_all_csvs)
}
