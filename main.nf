nextflow.enable.dsl=2

/*
 * PWA+NJ vs MSA+ML Regime Map Benchmark Main Workflow
 */

process SIMULATE_DATA {
    tag "D=${dist}_L=${len}_rep=${rep}"
    publishDir "${params.outdir}/replications/D${dist}_L${len}_rep${rep}", mode: 'copy'

    input:
    tuple val(dist), val(len), val(rep)

    output:
    tuple val(dist), val(len), val(rep), path("true_tree.nwk"), path("seqs.fasta"), emit: sim_data
    path("true_tree.nwk")
    path("seqs.fasta")

    script:
    """
    python3 ${projectDir}/bin/simulate_data.py \\
        --distance ${dist} \\
        --length ${len} \\
        --num_taxa ${params.num_taxa} \\
        --sigma ${params.sigma} \\
        --model "${params.model}" \\
        --alpha ${params.alpha} \\
        --indel_rate ${params.indel_rate} \\
        --outtree true_tree.nwk \\
        --outfasta seqs.fasta \\
        --seed ${rep}
    """
}

process RUN_PWA_NJ {
    tag "D=${dist}_L=${len}_rep=${rep}"
    publishDir "${params.outdir}/replications/D${dist}_L${len}_rep${rep}", mode: 'copy'

    input:
    tuple val(dist), val(len), val(rep), path(true_tree), path(fasta)

    output:
    path("pwa_nj_D${dist}_L${len}_rep${rep}.csv"), emit: csv
    path("pwa_nj.nwk")
    path("pwa_matrix.phylip")

    script:
    """
    python3 ${projectDir}/bin/run_pwa_nj.py \\
        --fasta ${fasta} \\
        --outtree pwa_nj.nwk \\
        --outmatrix pwa_matrix.phylip \\
        --gap_open ${params.gap_open} \\
        --gap_extend ${params.gap_extend} \\
        --dist_model ${params.dist_model} \\
        --alpha ${params.alpha} \\
        --tool ${params.nj_tool}

    python3 ${projectDir}/bin/evaluate_trees.py \\
        --truetree ${true_tree} \\
        --esttree pwa_nj.nwk \\
        --pipeline PWA+NJ \\
        --distance ${dist} \\
        --length ${len} \\
        --replicate ${rep} \\
        --outcsv pwa_nj_D${dist}_L${len}_rep${rep}.csv
    """
}

process RUN_MSA_NJ {
    tag "D=${dist}_L=${len}_rep=${rep}"
    publishDir "${params.outdir}/replications/D${dist}_L${len}_rep${rep}", mode: 'copy'

    input:
    tuple val(dist), val(len), val(rep), path(true_tree), path(fasta)

    output:
    path("msa_nj_D${dist}_L${len}_rep${rep}.csv"), emit: csv
    path("msa_nj.nwk")
    path("msa.fasta")
    path("msa_matrix.phylip")

    script:
    """
    python3 ${projectDir}/bin/run_msa_nj.py \\
        --fasta ${fasta} \\
        --outtree msa_nj.nwk \\
        --outmsa msa.fasta \\
        --outmatrix msa_matrix.phylip \\
        --dist_model ${params.dist_model} \\
        --alpha ${params.alpha} \\
        --tool ${params.nj_tool}

    python3 ${projectDir}/bin/evaluate_trees.py \\
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
    publishDir "${params.outdir}/replications/D${dist}_L${len}_rep${rep}", mode: 'copy'

    input:
    tuple val(dist), val(len), val(rep), path(true_tree), path(fasta)

    output:
    path("msa_ml_D${dist}_L${len}_rep${rep}.csv"), emit: csv
    path("msa_ml.nwk")
    path("msa_ml_meta.json")

    script:
    """
    python3 ${projectDir}/bin/run_msa_ml.py \\
        --fasta ${fasta} \\
        --outtree msa_ml.nwk \\
        --outjson msa_ml_meta.json \\
        --threads ${task.cpus}

    python3 ${projectDir}/bin/evaluate_trees.py \\
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

process COLLECT_AND_PLOT {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path csv_files

    output:
    path("benchmark_summary.csv")
    path("regime_map_delta_nrf.png")

    script:
    """
    cat *.csv | awk 'NR==1 || \$0 !~ /^distance/' > benchmark_summary.csv
    python3 ${projectDir}/bin/plot_regime_map.py --csv benchmark_summary.csv --outdir .
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

    RUN_PWA_NJ(ch_sim_data)
    RUN_MSA_NJ(ch_sim_data)
    RUN_MSA_ML(ch_sim_data)

    ch_all_csvs = RUN_PWA_NJ.out.csv.mix(RUN_MSA_NJ.out.csv).mix(RUN_MSA_ML.out.csv).collect()

    COLLECT_AND_PLOT(ch_all_csvs)
}
