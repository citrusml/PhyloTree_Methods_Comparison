nextflow.enable.dsl=2

/*
 * Unified Phylogenetic Benchmark Pipeline (next_main/main.nf)
 *
 * 全実験共通のコアワークフロー。
 * 実験ごとの条件（Indelモデル、置換モデル、ギャップペナルティ、距離・配列長・反復数など）は、
 * すべて `next_configs/*.config` の `params` ブロックで指定します。
 *
 * Pipelines supported:
 *   - PWA+NJ (BioPython Needleman-Wunsch + Poisson distance + RapidNJ)
 *   - MSA+NJ (MAFFT + Poisson distance + RapidNJ)
 *   - MSA+ML (MAFFT + IQ-TREE 2 ModelFinder)
 *   - TRUE_PWA+NJ (True alignment + Poisson distance + RapidNJ)
 *   - TRUE_MSA+NJ (True alignment + Poisson distance + RapidNJ)
 *   - TRUE_MSA+ML (True alignment + IQ-TREE 2 ModelFinder)
 *
 * Execution:
 *   nextflow run next_main/main.nf -c next_configs/your_experiment.config -profile supercomputer
 */

process SIMULATE_DATA {
    tag "D=${dist}_L=${len}_chk=${chunk_id}[${rep_start}..${rep_end}]"

    input:
    tuple val(dist), val(len), val(chunk_id), val(rep_start), val(rep_end)

    output:
    tuple val(dist), val(len), val(chunk_id), val(rep_start), val(rep_end), path("true_tree_*.nwk"), path("seqs_*.fasta"), emit: sim_data
    tuple val(dist), val(len), val(chunk_id), val(rep_start), val(rep_end), path("true_tree_*.nwk"), path("true_msa_*.fasta"), emit: true_msa_data

    script:
    def model_str = params.alpha ? "${params.model}{${params.alpha}}" : "${params.model}"
    def indel_size_arg = params.indel_size ? "--indel-size \"${params.indel_size}\"" : ""
    """
    IQTREE_BIN=\$(which iqtree2 2>/dev/null || which iqtree3 2>/dev/null || which iqtree 2>/dev/null)
    if [ -z "\$IQTREE_BIN" ]; then
        echo "Error: IQ-TREE / AliSim not found in PATH" >&2
        exit 1
    fi
    for rep in \$(seq ${rep_start} ${rep_end}); do
        \${IQTREE_BIN} --alisim sim_\${rep} \\
            -m "${model_str}" \\
            --length ${len} \\
            -t "RANDOM{bd{${params.birth_rate}/${params.death_rate}}/${params.taxa}}" \\
            --indel ${params.insert_rate},${params.delete_rate} \\
            ${indel_size_arg} \\
            --branch-scale ${dist} \\
            -af fasta \\
            -seed \${rep} \\
            --redo

        if [ -f sim_\${rep}.treefile ]; then
            mv sim_\${rep}.treefile true_tree_\${rep}.nwk
        else
            mv sim_\${rep}.tree true_tree_\${rep}.nwk
        fi
        mv sim_\${rep}.unaligned.fa seqs_\${rep}.fasta
        mv sim_\${rep}.fa true_msa_\${rep}.fasta
    done
    """
}

process RUN_PWA_NJ {
    tag "D=${dist}_L=${len}_chk=${chunk_id}[${rep_start}..${rep_end}]"

    input:
    tuple val(dist), val(len), val(chunk_id), val(rep_start), val(rep_end), path(true_trees), path(fastas)

    output:
    path("chunk_pwa_nj_D${dist}_L${len}_chk${chunk_id}.csv"), emit: csv
    path("pwa_nj_*.nwk")
    path("pwa_matrix_*.phylip")

    script:
    """
    for rep in \$(seq ${rep_start} ${rep_end}); do
        python3 ${projectDir}/../bin/run_pwa_nj.py \\
            --fasta seqs_\${rep}.fasta \\
            --outtree pwa_nj_\${rep}.nwk \\
            --outmatrix pwa_matrix_\${rep}.phylip \\
            --gap_open ${params.gap_open} \\
            --gap_extend ${params.gap_extend} \\
            --dist_model ${params.dist_model} \\
            --tool ${params.nj_tool} \\
            --threads ${task.cpus}

        python3 ${projectDir}/../bin/evaluate_trees.py \\
            --truetree true_tree_\${rep}.nwk \\
            --esttree pwa_nj_\${rep}.nwk \\
            --pipeline PWA+NJ \\
            --distance ${dist} \\
            --length ${len} \\
            --alpha ${params.alpha} \\
            --replicate \${rep} \\
            --outcsv chunk_pwa_nj_D${dist}_L${len}_chk${chunk_id}.csv
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
        mafft --threadit 0 --auto --thread ${task.cpus} --quiet seqs_\${rep}.fasta > msa_\${rep}.fasta
    done
    """
}

process RUN_MSA_NJ {
    tag "D=${dist}_L=${len}_chk=${chunk_id}[${rep_start}..${rep_end}]"

    input:
    tuple val(dist), val(len), val(chunk_id), val(rep_start), val(rep_end), path(true_trees), path(msas)

    output:
    path("chunk_msa_nj_D${dist}_L${len}_chk${chunk_id}.csv"), emit: csv
    path("msa_nj_*.nwk")
    path("msa_matrix_*.phylip")

    script:
    """
    for rep in \$(seq ${rep_start} ${rep_end}); do
        python3 ${projectDir}/../bin/run_msa_nj.py \\
            --msa msa_\${rep}.fasta \\
            --outtree msa_nj_\${rep}.nwk \\
            --outmatrix msa_matrix_\${rep}.phylip \\
            --dist_model ${params.dist_model} \\
            --tool ${params.nj_tool}

        python3 ${projectDir}/../bin/evaluate_trees.py \\
            --truetree true_tree_\${rep}.nwk \\
            --esttree msa_nj_\${rep}.nwk \\
            --pipeline MSA+NJ \\
            --distance ${dist} \\
            --length ${len} \\
            --alpha ${params.alpha} \\
            --replicate \${rep} \\
            --outcsv chunk_msa_nj_D${dist}_L${len}_chk${chunk_id}.csv
    done
    """
}

process RUN_MSA_ML {
    tag "D=${dist}_L=${len}_chk=${chunk_id}[${rep_start}..${rep_end}]"

    input:
    tuple val(dist), val(len), val(chunk_id), val(rep_start), val(rep_end), path(true_trees), path(msas)

    output:
    path("chunk_msa_ml_D${dist}_L${len}_chk${chunk_id}.csv"), emit: csv
    path("msa_ml_*.nwk")
    path("msa_ml_meta_*.json")

    script:
    """
    for rep in \$(seq ${rep_start} ${rep_end}); do
        python3 ${projectDir}/../bin/run_msa_ml.py \\
            --msa msa_\${rep}.fasta \\
            --outtree msa_ml_\${rep}.nwk \\
            --outjson msa_ml_meta_\${rep}.json \\
            --threads ${task.cpus}

        python3 ${projectDir}/../bin/evaluate_trees.py \\
            --truetree true_tree_\${rep}.nwk \\
            --esttree msa_ml_\${rep}.nwk \\
            --pipeline MSA+ML \\
            --distance ${dist} \\
            --length ${len} \\
            --alpha ${params.alpha} \\
            --replicate \${rep} \\
            --json msa_ml_meta_\${rep}.json \\
            --outcsv chunk_msa_ml_D${dist}_L${len}_chk${chunk_id}.csv
    done
    """
}

process RUN_TRUE_PWA_NJ {
    tag "D=${dist}_L=${len}_chk=${chunk_id}[${rep_start}..${rep_end}]"

    input:
    tuple val(dist), val(len), val(chunk_id), val(rep_start), val(rep_end), path(true_trees), path(true_msas)

    output:
    path("chunk_true_pwa_nj_D${dist}_L${len}_chk${chunk_id}.csv"), emit: csv
    path("true_pwa_nj_*.nwk")
    path("true_pwa_matrix_*.phylip")

    script:
    """
    for rep in \$(seq ${rep_start} ${rep_end}); do
        python3 ${projectDir}/../bin/run_msa_nj.py \\
            --msa true_msa_\${rep}.fasta \\
            --outtree true_pwa_nj_\${rep}.nwk \\
            --outmatrix true_pwa_matrix_\${rep}.phylip \\
            --dist_model ${params.dist_model} \\
            --tool ${params.nj_tool}

        python3 ${projectDir}/../bin/evaluate_trees.py \\
            --truetree true_tree_\${rep}.nwk \\
            --esttree true_pwa_nj_\${rep}.nwk \\
            --pipeline TRUE_PWA+NJ \\
            --distance ${dist} \\
            --length ${len} \\
            --alpha ${params.alpha} \\
            --replicate \${rep} \\
            --outcsv chunk_true_pwa_nj_D${dist}_L${len}_chk${chunk_id}.csv
    done
    """
}

process RUN_TRUE_MSA_NJ {
    tag "D=${dist}_L=${len}_chk=${chunk_id}[${rep_start}..${rep_end}]"

    input:
    tuple val(dist), val(len), val(chunk_id), val(rep_start), val(rep_end), path(true_trees), path(true_msas)

    output:
    path("chunk_true_msa_nj_D${dist}_L${len}_chk${chunk_id}.csv"), emit: csv
    path("true_msa_nj_*.nwk")
    path("true_msa_matrix_*.phylip")

    script:
    """
    for rep in \$(seq ${rep_start} ${rep_end}); do
        python3 ${projectDir}/../bin/run_msa_nj.py \\
            --msa true_msa_\${rep}.fasta \\
            --outtree true_msa_nj_\${rep}.nwk \\
            --outmatrix true_msa_matrix_\${rep}.phylip \\
            --dist_model ${params.dist_model} \\
            --tool ${params.nj_tool}

        python3 ${projectDir}/../bin/evaluate_trees.py \\
            --truetree true_tree_\${rep}.nwk \\
            --esttree true_msa_nj_\${rep}.nwk \\
            --pipeline TRUE_MSA+NJ \\
            --distance ${dist} \\
            --length ${len} \\
            --alpha ${params.alpha} \\
            --replicate \${rep} \\
            --outcsv chunk_true_msa_nj_D${dist}_L${len}_chk${chunk_id}.csv
    done
    """
}

process RUN_TRUE_MSA_ML {
    tag "D=${dist}_L=${len}_chk=${chunk_id}[${rep_start}..${rep_end}]"

    input:
    tuple val(dist), val(len), val(chunk_id), val(rep_start), val(rep_end), path(true_trees), path(true_msas)

    output:
    path("chunk_true_msa_ml_D${dist}_L${len}_chk${chunk_id}.csv"), emit: csv
    path("true_msa_ml_*.nwk")
    path("true_msa_ml_meta_*.json")

    script:
    """
    for rep in \$(seq ${rep_start} ${rep_end}); do
        python3 ${projectDir}/../bin/run_msa_ml.py \\
            --msa true_msa_\${rep}.fasta \\
            --outtree true_msa_ml_\${rep}.nwk \\
            --outjson true_msa_ml_meta_\${rep}.json \\
            --threads ${task.cpus}

        python3 ${projectDir}/../bin/evaluate_trees.py \\
            --truetree true_tree_\${rep}.nwk \\
            --esttree true_msa_ml_\${rep}.nwk \\
            --pipeline TRUE_MSA+ML \\
            --distance ${dist} \\
            --length ${len} \\
            --alpha ${params.alpha} \\
            --replicate \${rep} \\
            --json true_msa_ml_meta_\${rep}.json \\
            --outcsv chunk_true_msa_ml_D${dist}_L${len}_chk${chunk_id}.csv
    done
    """
}

process COLLECT_AND_PLOT {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path csv_files

    output:
    path("benchmark_summary.csv")
    path("regime_map_delta_nrf.png"), optional: true
    path("nrf_boxplots.png"), optional: true
    path("summary_statistics.csv"), optional: true
    path("method_comparisons"), optional: true
    path("length_analysis"), optional: true

    script:
    """
    head -n 1 \$(ls *.csv | head -1) > benchmark_summary.csv
    for f in *.csv; do tail -n +2 "\$f"; done >> benchmark_summary.csv
    python3 ${projectDir}/../bin/plot/generate_all_reports.py --csv benchmark_summary.csv --outdir .
    """
}

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
