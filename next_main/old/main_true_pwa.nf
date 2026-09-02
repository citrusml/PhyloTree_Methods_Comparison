nextflow.enable.dsl=2

/*
 * True Pairwise Alignment (TRUE_PWA+NJ) Benchmark (Chunked / Batched Execution)
 * Groups 10 replicates into a single task to maximize HPC throughput.
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
params.gap_open    = 10.0
params.gap_extend  = 0.5
params.dist_model  = "poisson"
params.outdir      = "results/results_true_pwa"
params.nj_tool     = "rapidnj"

process SIMULATE_DATA {
    tag "D=${dist}_L=${len}_chk=${chunk_id}[${rep_start}..${rep_end}]"

    input:
    tuple val(dist), val(len), val(chunk_id), val(rep_start), val(rep_end)

    output:
    tuple val(dist), val(len), val(chunk_id), val(rep_start), val(rep_end), path("true_tree_*.nwk"), path("true_msa_*.fasta"), emit: true_msa_data

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
            --seed \${rep}
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
            --alpha ${params.alpha} \\
            --tool ${params.nj_tool}

        python3 ${projectDir}/../bin/evaluate_trees.py \\
            --truetree true_tree_\${rep}.nwk \\
            --esttree true_pwa_nj_\${rep}.nwk \\
            --pipeline TRUE_PWA+NJ \\
            --distance ${dist} \\
            --length ${len} \\
            --replicate \${rep} \\
            --outcsv chunk_true_pwa_nj_D${dist}_L${len}_chk${chunk_id}.csv
    done
    """
}

process COLLECT_AND_PLOT {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path csv_files

    output:
    path("benchmark_true_pwa_summary.csv")
    path("true_pwa_heatmap.png")
    path("true_pwa_heatmap.pdf"), optional: true
    path("true_pwa_scaling_curves.png")
    path("true_pwa_scaling_curves.pdf"), optional: true
    path("benchmark_true_pwa_statistics.csv")

    script:
    """
    cat *.csv | awk 'NR==1 || \$0 !~ /^distance/' > benchmark_true_pwa_summary.csv
    python3 ${projectDir}/../bin/plot_true_pwa_benchmark.py --csv benchmark_true_pwa_summary.csv --outdir .
    """
}

workflow {
    ch_distances = Channel.fromList(params.distances)
    ch_lengths   = Channel.fromList(params.lengths)

    // Calculate chunks based on replicates and chunk_size
    def num_chunks = Math.max(1, (params.replicates / params.chunk_size) as int)
    ch_chunks = Channel.fromList( (1..num_chunks).collect { c ->
        def r_start = (c - 1) * params.chunk_size + 1
        def r_end   = Math.min(params.replicates, c * params.chunk_size)
        [c, r_start, r_end]
    } )

    ch_params = ch_distances.combine(ch_lengths).combine(ch_chunks)

    SIMULATE_DATA(ch_params)

    RUN_TRUE_PWA_NJ(SIMULATE_DATA.out.true_msa_data)

    COLLECT_AND_PLOT(RUN_TRUE_PWA_NJ.out.csv.collect())
}
