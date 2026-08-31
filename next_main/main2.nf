nextflow.enable.dsl=2

/*
 * Method Comparison Benchmark Pipeline 2 (main2.nf)
 * Evaluates PWA+NJ, MSA+NJ, MSA+ML, TRUE_MSA+NJ, TRUE_MSA+ML, and TRUE_DIST+NJ
 * under the updated MIN_OVERLAP = 10 pairwise distance calculation.
 * Chunked / Batched Execution: Groups 50 replicates into a single task to maximize HPC throughput.
 */
params.taxa        = 32
params.distances   = [0.1, 0.5, 1.0, 2.0, 3.0]
params.lengths     = [100, 300, 500, 1000, 1500]
params.replicates  = 100
params.chunk_size  = 50
params.birth_rate  = 0.1
params.death_rate  = 0.05
params.insert_rate = 0.05
params.delete_rate = 0.10
params.model       = "LG+G4"
params.alpha       = 1.0
params.dist_model  = "poisson"
params.gap_open    = 10.0
params.gap_extend  = 0.5
params.run_true_msa  = true
params.run_true_dist = true
params.outdir      = "results/results2"
params.nj_tool     = "rapidnj"

process SIMULATE_DATA {
    tag "D=${dist}_L=${len}_chk=${chunk_id}[${rep_start}..${rep_end}]"

    input:
    tuple val(dist), val(len), val(chunk_id), val(rep_start), val(rep_end)

    output:
    tuple val(dist), val(len), val(chunk_id), val(rep_start), val(rep_end), path("true_tree_*.nwk"), path("seqs_*.fasta"), emit: sim_data
    tuple val(dist), val(len), val(chunk_id), val(rep_start), val(rep_end), path("true_tree_*.nwk"), path("true_msa_*.fasta"), emit: true_msa_data
    tuple val(dist), val(len), val(chunk_id), val(rep_start), val(rep_end), path("true_tree_*.nwk"), path("true_matrix_*.phylip"), emit: true_dist_data

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
            --outtrue_msa true_msa_\${rep}.fasta \\
            --outtrue_matrix true_matrix_\${rep}.phylip \\
            --seed \${rep}
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
        python3 ${projectDir}/../bin/run_msa.py \\
            --fasta seqs_\${rep}.fasta \\
            --outmsa msa_\${rep}.fasta \\
            --threads ${task.cpus}
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

process RUN_TRUE_DIST_NJ {
    tag "D=${dist}_L=${len}_chk=${chunk_id}[${rep_start}..${rep_end}]"

    input:
    tuple val(dist), val(len), val(chunk_id), val(rep_start), val(rep_end), path(true_trees), path(true_matrices)

    output:
    path("chunk_true_dist_nj_D${dist}_L${len}_chk${chunk_id}.csv"), emit: csv
    path("true_dist_nj_*.nwk")

    script:
    """
    for rep in \$(seq ${rep_start} ${rep_end}); do
        python3 ${projectDir}/../bin/run_pwa_nj.py \\
            --matrix true_matrix_\${rep}.phylip \\
            --outtree true_dist_nj_\${rep}.nwk \\
            --tool ${params.nj_tool}

        python3 ${projectDir}/../bin/evaluate_trees.py \\
            --truetree true_tree_\${rep}.nwk \\
            --esttree true_dist_nj_\${rep}.nwk \\
            --pipeline TRUE_DIST+NJ \\
            --distance ${dist} \\
            --length ${len} \\
            --alpha ${params.alpha} \\
            --replicate \${rep} \\
            --outcsv chunk_true_dist_nj_D${dist}_L${len}_chk${chunk_id}.csv
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
    # Concatenate CSVs: keep header from first file only
    head -n 1 \$(ls *.csv | head -1) > benchmark_summary.csv
    for f in *.csv; do tail -n +2 "\$f"; done >> benchmark_summary.csv
    python3 ${projectDir}/../bin/plot/generate_all_reports.py --csv benchmark_summary.csv --outdir .
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

    // 1. PWA+NJ pipeline
    RUN_PWA_NJ(ch_sim_data)

    // 2. Shared MAFFT calculation
    RUN_MAFFT(ch_sim_data)
    ch_msa_data = RUN_MAFFT.out.msa_data

    // 3. MSA+NJ and MSA+ML pipelines
    RUN_MSA_NJ(ch_msa_data)
    RUN_MSA_ML(ch_msa_data)

    ch_all_csvs = RUN_PWA_NJ.out.csv.mix(RUN_MSA_NJ.out.csv).mix(RUN_MSA_ML.out.csv)

    // 4. True MSA + NJ / True MSA + ML evaluation
    if (params.run_true_msa) {
        ch_true_msa_data = SIMULATE_DATA.out.true_msa_data
        RUN_TRUE_MSA_NJ(ch_true_msa_data)
        RUN_TRUE_MSA_ML(ch_true_msa_data)
        ch_all_csvs = ch_all_csvs.mix(RUN_TRUE_MSA_NJ.out.csv).mix(RUN_TRUE_MSA_ML.out.csv)
    }

    // 5. True Patristic Distance + NJ evaluation
    if (params.run_true_dist) {
        ch_true_dist_data = SIMULATE_DATA.out.true_dist_data
        RUN_TRUE_DIST_NJ(ch_true_dist_data)
        ch_all_csvs = ch_all_csvs.mix(RUN_TRUE_DIST_NJ.out.csv)
    }

    COLLECT_AND_PLOT(ch_all_csvs.collect())
}
