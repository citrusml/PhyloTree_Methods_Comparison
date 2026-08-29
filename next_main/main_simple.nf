nextflow.enable.dsl=2

/*
 * Simple Condition(LG model) Benchmark Pipeiline (Experiment 9)
 * Evaluates PWA+NJ, MSA+NJ, and MSA+ML.
 * Chunked / Batched Execution: Groups 20 replicates into a single task to maximize HPC throughput.
 */
params.taxa        = 32
params.distances   = [0.1, 0.5, 1.0, 2.0, 3.0]
params.lengths     = [100, 300, 500, 1000, 1500]
params.replicates  = 100
params.chunk_size  = 20
params.birth_rate  = 0.1
params.death_rate  = 0.05
params.insert_rate = 0.05
params.delete_rate = 0.10
params.model       = "LG"
alpha              = 1.0  
params.dist_model  = "poisson"
params.gap_open    = 10.0
params.gap_extend  = 0.5
params.outdir      = "results/results_simple"
params.nj_tool     = "rapidnj"

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

process RUN_PWA_NJ {
    tag "D=${dist}_L=${len}_chk=${chunk_id}[${rep_start}..${rep_end}]"

    input:
    tuple val(dist), val(len), val(chunk_id), val(rep_start), val(rep_end), path(true_trees), path(fastas)

    output:
    path("chunk_pwa_nj_high_dist_D${dist}_L${len}_chk${chunk_id}.csv"), emit: csv
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
            --outcsv chunk_pwa_nj_high_dist_D${dist}_L${len}_chk${chunk_id}.csv
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
    path("chunk_msa_nj_high_dist_D${dist}_L${len}_chk${chunk_id}.csv"), emit: csv
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
            --outcsv chunk_msa_nj_high_dist_D${dist}_L${len}_chk${chunk_id}.csv
    done
    """
}

process RUN_MSA_ML {
    tag "D=${dist}_L=${len}_chk=${chunk_id}[${rep_start}..${rep_end}]"

    input:
    tuple val(dist), val(len), val(chunk_id), val(rep_start), val(rep_end), path(true_trees), path(msas)

    output:
    path("chunk_msa_ml_high_dist_D${dist}_L${len}_chk${chunk_id}.csv"), emit: csv
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
            --outcsv chunk_msa_ml_high_dist_D${dist}_L${len}_chk${chunk_id}.csv
    done
    """
}

process COLLECT_AND_PLOT {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path csv_files

    output:
    path("benchmark_high_dist_summary.csv")
    path("scaling_curves_high_dist.png")
    path("scaling_curves_high_dist.pdf"), optional: true
    path("regime_map_high_dist.png"), optional: true
    path("regime_map_high_dist.pdf"), optional: true
    path("benchmark_high_dist_statistics.csv")

    script:
    """
    cat *.csv | awk 'NR==1 || \$0 !~ /^alpha/' > benchmark_high_dist_summary.csv
    python3 ${projectDir}/../bin/plot_high_dist_benchmark.py --csv benchmark_high_dist_summary.csv --outdir .
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

    ch_all_csvs = RUN_PWA_NJ.out.csv.mix(RUN_MSA_NJ.out.csv).mix(RUN_MSA_ML.out.csv).collect()

    COLLECT_AND_PLOT(ch_all_csvs)
}
