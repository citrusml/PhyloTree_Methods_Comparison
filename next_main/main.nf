nextflow.enable.dsl=2

/*
 * PWA+NJ vs MSA+ML Regime Map Benchmark Main Workflow
 * Default parameter fallback definitions
 */
params.distances     = [0.1, 0.5, 1.0, 2.0, 3.0]
params.lengths       = [100, 300, 500, 1000, 1500]
params.replicates    = 100
params.num_taxa      = 32
params.birth_rate    = 0.1
params.death_rate    = 0.05
params.insert_rate   = 0.05
params.delete_rate   = 0.10
params.model         = "LG+G4"
params.alpha         = 1.0
params.gap_open      = 10.0
params.gap_extend    = 0.5
params.dist_model    = "poisson"
params.run_true_msa  = true
params.run_true_dist = true
params.outdir        = "results"
params.nj_tool       = "rapidnj"

process SIMULATE_DATA {
    tag "D=${dist}_L=${len}_rep=${rep}"

    input:
    tuple val(dist), val(len), val(rep)

    output:
    tuple val(dist), val(len), val(rep), path("true_tree.nwk"), path("seqs.fasta"), emit: sim_data
    tuple val(dist), val(len), val(rep), path("true_tree.nwk"), path("true_msa.fasta"), emit: true_msa_data
    tuple val(dist), val(len), val(rep), path("true_tree.nwk"), path("true_matrix.phylip"), emit: true_dist_data

    script:
    """
    python3 ${projectDir}/../bin/simulate_data.py \\
        --distance ${dist} \\
        --length ${len} \\
        --num_taxa ${params.num_taxa} \\
        --birth_rate ${params.birth_rate} \\
        --death_rate ${params.death_rate} \\
        --insert_rate ${params.insert_rate} \\
        --delete_rate ${params.delete_rate} \\
        --model "${params.model}" \\
        --alpha ${params.alpha} \\
        --outtree true_tree.nwk \\
        --outfasta seqs.fasta \\
        --outtrue_msa true_msa.fasta \\
        --outtrue_matrix true_matrix.phylip \\
        --seed ${rep}
    """
}

process RUN_PWA_NJ {
    tag "D=${dist}_L=${len}_rep=${rep}"

    input:
    tuple val(dist), val(len), val(rep), path(true_tree), path(fasta)

    output:
    path("pwa_nj_D${dist}_L${len}_rep${rep}.csv"), emit: csv
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
        --alpha ${params.alpha} \\
        --tool ${params.nj_tool}

    python3 ${projectDir}/../bin/evaluate_trees.py \\
        --truetree ${true_tree} \\
        --esttree pwa_nj.nwk \\
        --pipeline PWA+NJ \\
        --distance ${dist} \\
        --length ${len} \\
        --replicate ${rep} \\
        --outcsv pwa_nj_D${dist}_L${len}_rep${rep}.csv
    """
}

process RUN_MAFFT {
    tag "D=${dist}_L=${len}_rep=${rep}"

    input:
    tuple val(dist), val(len), val(rep), path(true_tree), path(fasta)

    output:
    tuple val(dist), val(len), val(rep), path(true_tree), path("msa.fasta"), emit: msa_data

    script:
    """
    python3 ${projectDir}/../bin/run_msa.py \\
        --fasta ${fasta} \\
        --outmsa msa.fasta \\
        --threads ${task.cpus}
    """
}

process RUN_MSA_NJ {
    tag "D=${dist}_L=${len}_rep=${rep}"

    input:
    tuple val(dist), val(len), val(rep), path(true_tree), path(msa)

    output:
    path("msa_nj_D${dist}_L${len}_rep${rep}.csv"), emit: csv
    path("msa_nj.nwk")
    path("msa_matrix.phylip")

    script:
    """
    python3 ${projectDir}/../bin/run_msa_nj.py \\
        --msa ${msa} \\
        --outtree msa_nj.nwk \\
        --outmatrix msa_matrix.phylip \\
        --dist_model ${params.dist_model} \\
        --alpha ${params.alpha} \\
        --tool ${params.nj_tool}

    python3 ${projectDir}/../bin/evaluate_trees.py \\
        --truetree ${true_tree} \\
        --esttree msa_nj.nwk \\
        --pipeline MSA+NJ \\
        --distance ${dist} \\
        --length ${len} \\
        --replicate ${rep} \\
        --outcsv msa_nj_D${dist}_L${len}_rep${rep}.csv
    """
}

process RUN_MSA_ML {
    tag "D=${dist}_L=${len}_rep=${rep}"

    input:
    tuple val(dist), val(len), val(rep), path(true_tree), path(msa)

    output:
    path("msa_ml_D${dist}_L${len}_rep${rep}.csv"), emit: csv
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
        --replicate ${rep} \\
        --json msa_ml_meta.json \\
        --outcsv msa_ml_D${dist}_L${len}_rep${rep}.csv
    """
}

process RUN_TRUE_MSA_NJ {
    tag "D=${dist}_L=${len}_rep=${rep}"

    input:
    tuple val(dist), val(len), val(rep), path(true_tree), path(true_msa)

    output:
    path("true_msa_nj_D${dist}_L${len}_rep${rep}.csv"), emit: csv
    path("true_msa_nj.nwk")
    path("true_msa_matrix.phylip")

    script:
    """
    python3 ${projectDir}/../bin/run_msa_nj.py \\
        --msa ${true_msa} \\
        --outtree true_msa_nj.nwk \\
        --outmatrix true_msa_matrix.phylip \\
        --dist_model ${params.dist_model} \\
        --alpha ${params.alpha} \\
        --tool ${params.nj_tool}

    python3 ${projectDir}/../bin/evaluate_trees.py \\
        --truetree ${true_tree} \\
        --esttree true_msa_nj.nwk \\
        --pipeline TRUE_MSA+NJ \\
        --distance ${dist} \\
        --length ${len} \\
        --replicate ${rep} \\
        --outcsv true_msa_nj_D${dist}_L${len}_rep${rep}.csv
    """
}

process RUN_TRUE_MSA_ML {
    tag "D=${dist}_L=${len}_rep=${rep}"

    input:
    tuple val(dist), val(len), val(rep), path(true_tree), path(true_msa)

    output:
    path("true_msa_ml_D${dist}_L${len}_rep${rep}.csv"), emit: csv
    path("true_msa_ml.nwk")
    path("true_msa_ml_meta.json")

    script:
    """
    python3 ${projectDir}/../bin/run_msa_ml.py \\
        --msa ${true_msa} \\
        --outtree true_msa_ml.nwk \\
        --outjson true_msa_ml_meta.json \\
        --threads ${task.cpus}

    python3 ${projectDir}/../bin/evaluate_trees.py \\
        --truetree ${true_tree} \\
        --esttree true_msa_ml.nwk \\
        --pipeline TRUE_MSA+ML \\
        --distance ${dist} \\
        --length ${len} \\
        --replicate ${rep} \\
        --json true_msa_ml_meta.json \\
        --outcsv true_msa_ml_D${dist}_L${len}_rep${rep}.csv
    """
}

process RUN_TRUE_DIST_NJ {
    tag "D=${dist}_L=${len}_rep=${rep}"

    input:
    tuple val(dist), val(len), val(rep), path(true_tree), path(true_matrix)

    output:
    path("true_dist_nj_D${dist}_L${len}_rep${rep}.csv"), emit: csv
    path("true_dist_nj.nwk")

    script:
    """
    python3 ${projectDir}/../bin/run_pwa_nj.py \\
        --matrix ${true_matrix} \\
        --outtree true_dist_nj.nwk \\
        --tool ${params.nj_tool}

    python3 ${projectDir}/../bin/evaluate_trees.py \\
        --truetree ${true_tree} \\
        --esttree true_dist_nj.nwk \\
        --pipeline TRUE_DIST+NJ \\
        --distance ${dist} \\
        --length ${len} \\
        --replicate ${rep} \\
        --outcsv true_dist_nj_D${dist}_L${len}_rep${rep}.csv
    """
}

process COLLECT_AND_PLOT {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path csv_files

    output:
    path("benchmark_summary.csv")
    path("regime_map_delta_nrf.png")
    path("nrf_boxplots.png")
    path("summary_statistics.csv")
    path("method_comparisons/*") optional true

    script:
    """
    cat *.csv | awk 'NR==1 || \$0 !~ /^distance/' > benchmark_summary.csv
    python3 ${projectDir}/../bin/generate_all_reports.py --csv benchmark_summary.csv --outdir .
    """
}

workflow {
    // Generate combination channel for D x L x replicate
    ch_distances = Channel.fromList(params.distances)
    ch_lengths   = Channel.fromList(params.lengths)
    ch_replicates = Channel.from(1..params.replicates)

    ch_params = ch_distances.combine(ch_lengths).combine(ch_replicates)

    SIMULATE_DATA(ch_params)

    ch_sim_data = SIMULATE_DATA.out.sim_data

    // 1. PWA+NJ pipeline
    RUN_PWA_NJ(ch_sim_data)

    // 2. Shared MAFFT MSA calculation
    RUN_MAFFT(ch_sim_data)
    ch_msa_data = RUN_MAFFT.out.msa_data

    // 3. MSA+NJ and MSA+ML pipelines (reuse shared MSA)
    RUN_MSA_NJ(ch_msa_data)
    RUN_MSA_ML(ch_msa_data)

    ch_csvs = RUN_PWA_NJ.out.csv.mix(RUN_MSA_NJ.out.csv).mix(RUN_MSA_ML.out.csv)

    // 4. Optional: True MSA + NJ / True MSA + ML evaluation
    if (params.run_true_msa) {
        ch_true_msa_data = SIMULATE_DATA.out.true_msa_data
        RUN_TRUE_MSA_NJ(ch_true_msa_data)
        RUN_TRUE_MSA_ML(ch_true_msa_data)
        ch_csvs = ch_csvs.mix(RUN_TRUE_MSA_NJ.out.csv).mix(RUN_TRUE_MSA_ML.out.csv)
    }

    // 5. Optional: True Patristic Distance + NJ evaluation
    if (params.run_true_dist) {
        ch_true_dist_data = SIMULATE_DATA.out.true_dist_data
        RUN_TRUE_DIST_NJ(ch_true_dist_data)
        ch_csvs = ch_csvs.mix(RUN_TRUE_DIST_NJ.out.csv)
    }

    COLLECT_AND_PLOT(ch_csvs.collect())
}
