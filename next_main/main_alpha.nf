nextflow.enable.dsl=2

/*
 * PWA+NJ vs MSA+ML Rate Heterogeneity (Alpha) Benchmark (alpha = 0.25, 0.5, 1.0, 2.0)
 * Default parameter fallback definitions
 */
params.taxa        = 32
params.distances   = [0.1, 0.5, 1.0, 2.0, 3.0]
params.lengths     = [100, 500, 1000]
params.replicates  = 100
params.birth_rate  = 0.1
params.death_rate  = 0.05
params.insert_rate = 0.05
params.delete_rate = 0.10
params.model       = "LG+G4"
params.alpha       = [0.25, 0.5, 1.0, 2.0]
params.gap_open    = 10.0
params.gap_extend  = 0.5
params.dist_model  = "poisson"
params.outdir      = "results/results_alpha"
params.nj_tool     = "rapidnj"

process SIMULATE_DATA {
    tag "A=${alpha}_D=${dist}_L=${len}_rep=${rep}"

    input:
    tuple val(alpha), val(dist), val(len), val(rep)

    output:
    tuple val(alpha), val(dist), val(len), val(rep), path("true_tree.nwk"), path("seqs.fasta"), emit: sim_data

    script:
    """
    python3 ${projectDir}/../bin/simulate_data.py \\
        --num_taxa ${params.taxa} \\
        --distance ${dist} \\
        --length ${len} \\
        --birth_rate ${params.birth_rate} \\
        --death_rate ${params.death_rate} \\
        --insert_rate ${params.insert_rate} \\
        --delete_rate ${params.delete_rate} \\
        --model "${params.model}" \\
        --alpha ${alpha} \\
        --outtree true_tree.nwk \\
        --outfasta seqs.fasta \\
        --seed ${rep}
    """
}

process RUN_PWA_NJ {
    tag "A=${alpha}_D=${dist}_L=${len}_rep=${rep}"

    input:
    tuple val(alpha), val(dist), val(len), val(rep), path(true_tree), path(fasta)

    output:
    path("pwa_nj_A${alpha}_D${dist}_L${len}_rep${rep}.csv"), emit: csv
    path("pwa_nj.nwk")
    path("pwa_matrix.phylip")

    script:
    """
    python3 ${projectDir}/../bin/run_pwa_nj.py \\
        --fasta ${fasta} \\
        --outtree pwa_nj.nwk \\
        --outmatrix pwa_matrix.phylip \\
        --gap_open ${params.gap_open} \\
        --gap_extend ${params.gap_extend} \\
        --dist_model ${params.dist_model} \\
        --alpha ${alpha} \\
        --tool ${params.nj_tool} \\
        --threads ${task.cpus}

    python3 ${projectDir}/../bin/evaluate_trees.py \\
        --truetree ${true_tree} \\
        --esttree pwa_nj.nwk \\
        --pipeline PWA+NJ \\
        --distance ${dist} \\
        --length ${len} \\
        --alpha ${alpha} \\
        --replicate ${rep} \\
        --outcsv pwa_nj_A${alpha}_D${dist}_L${len}_rep${rep}.csv
    """
}

process RUN_MAFFT {
    tag "A=${alpha}_D=${dist}_L=${len}_rep=${rep}"

    input:
    tuple val(alpha), val(dist), val(len), val(rep), path(true_tree), path(fasta)

    output:
    tuple val(alpha), val(dist), val(len), val(rep), path(true_tree), path("msa.fasta"), emit: msa_data

    script:
    """
    python3 ${projectDir}/../bin/run_msa.py \\
        --fasta ${fasta} \\
        --outmsa msa.fasta \\
        --threads ${task.cpus}
    """
}

process RUN_MSA_NJ {
    tag "A=${alpha}_D=${dist}_L=${len}_rep=${rep}"

    input:
    tuple val(alpha), val(dist), val(len), val(rep), path(true_tree), path(msa)

    output:
    path("msa_nj_A${alpha}_D${dist}_L${len}_rep${rep}.csv"), emit: csv
    path("msa_nj.nwk")
    path("msa_matrix.phylip")

    script:
    """
    python3 ${projectDir}/../bin/run_msa_nj.py \\
        --msa ${msa} \\
        --outtree msa_nj.nwk \\
        --outmatrix msa_matrix.phylip \\
        --dist_model ${params.dist_model} \\
        --alpha ${alpha} \\
        --tool ${params.nj_tool}

    python3 ${projectDir}/../bin/evaluate_trees.py \\
        --truetree ${true_tree} \\
        --esttree msa_nj.nwk \\
        --pipeline MSA+NJ \\
        --distance ${dist} \\
        --length ${len} \\
        --alpha ${alpha} \\
        --replicate ${rep} \\
        --outcsv msa_nj_A${alpha}_D${dist}_L${len}_rep${rep}.csv
    """
}

process RUN_MSA_ML {
    tag "A=${alpha}_D=${dist}_L=${len}_rep=${rep}"

    input:
    tuple val(alpha), val(dist), val(len), val(rep), path(true_tree), path(msa)

    output:
    path("msa_ml_A${alpha}_D${dist}_L${len}_rep${rep}.csv"), emit: csv
    path("msa_ml.nwk")
    path("msa_ml_meta.json")

    script:
    """
    python3 ${projectDir}/../bin/run_msa_ml.py \\
        --msa ${msa} \\
        --outtree msa_ml.nwk \\
        --outjson msa_ml_meta.json \\
        --threads ${task.cpus}

    python3 ${projectDir}/../bin/evaluate_trees.py \\
        --truetree ${true_tree} \\
        --esttree msa_ml.nwk \\
        --pipeline MSA+ML \\
        --distance ${dist} \\
        --length ${len} \\
        --alpha ${alpha} \\
        --replicate ${rep} \\
        --json msa_ml_meta.json \\
        --outcsv msa_ml_A${alpha}_D${dist}_L${len}_rep${rep}.csv
    """
}

process COLLECT_AND_PLOT {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path csv_files

    output:
    path("benchmark_alpha_summary.csv")
    path("alpha_scaling_curves.png")
    path("alpha_scaling_curves.pdf"), optional: true
    path("benchmark_alpha_statistics.csv")

    script:
    """
    cat *.csv | awk 'NR==1 || \$0 !~ /^alpha/' > benchmark_alpha_summary.csv
    python3 ${projectDir}/../bin/plot_alpha_benchmark.py --csv benchmark_alpha_summary.csv --outdir .
    """
}

workflow {
    ch_alphas     = Channel.fromList(params.alpha instanceof List ? params.alpha : [params.alpha])
    ch_distances  = Channel.fromList(params.distances)
    ch_lengths    = Channel.fromList(params.lengths)
    ch_replicates = Channel.from(1..params.replicates)

    ch_params = ch_alphas.combine(ch_distances).combine(ch_lengths).combine(ch_replicates)

    SIMULATE_DATA(ch_params)

    ch_sim_data = SIMULATE_DATA.out.sim_data

    // 1. PWA+NJ pipeline
    RUN_PWA_NJ(ch_sim_data)

    // 2. Shared MAFFT calculation
    RUN_MAFFT(ch_sim_data)
    ch_msa_data = RUN_MAFFT.out.msa_data

    // 3. MSA+NJ and MSA+ML pipelines (reuse shared MSA)
    RUN_MSA_NJ(ch_msa_data)
    RUN_MSA_ML(ch_msa_data)

    ch_all_csvs = RUN_PWA_NJ.out.csv.mix(RUN_MSA_NJ.out.csv).mix(RUN_MSA_ML.out.csv).collect()

    COLLECT_AND_PLOT(ch_all_csvs)
}
