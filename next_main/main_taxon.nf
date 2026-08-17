nextflow.enable.dsl=2

/*
 * PWA+NJ vs MSA+ML Taxon Scaling Benchmark (N = 8, 32, 64, 128)
 * Default parameter fallback definitions
 */
params.taxa       = [8, 32, 64, 128]
params.distances  = [0.1, 0.5, 1.0, 2.0, 3.0]
params.lengths    = [100, 500, 1000]
params.replicates = 100
params.sigma      = 0.5
params.model      = "LG+G"
params.alpha      = 1.0
params.indel_rate = 0.01
params.gap_open   = 10.0
params.gap_extend = 0.5
params.dist_model = "poisson"
params.outdir     = "results/results_taxon"
params.nj_tool    = "rapidnj"

process SIMULATE_DATA {
    tag "N=${taxa}_D=${dist}_L=${len}_rep=${rep}"

    input:
    tuple val(taxa), val(dist), val(len), val(rep)

    output:
    tuple val(taxa), val(dist), val(len), val(rep), path("true_tree.nwk"), path("seqs.fasta"), emit: sim_data

    script:
    """
    python3 ${projectDir}/../bin/simulate_data.py \\
        --num_taxa ${taxa} \\
        --distance ${dist} \\
        --length ${len} \\
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
    tag "N=${taxa}_D=${dist}_L=${len}_rep=${rep}"

    input:
    tuple val(taxa), val(dist), val(len), val(rep), path(true_tree), path(fasta)

    output:
    path("pwa_nj_N${taxa}_D${dist}_L${len}_rep${rep}.csv"), emit: csv
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
        --tool ${params.nj_tool} \\
        --threads ${task.cpus}

    python3 ${projectDir}/../bin/evaluate_trees.py \\
        --truetree ${true_tree} \\
        --esttree pwa_nj.nwk \\
        --pipeline PWA+NJ \\
        --distance ${dist} \\
        --length ${len} \\
        --replicate ${rep} \\
        --outcsv pwa_nj_N${taxa}_D${dist}_L${len}_rep${rep}.csv
    """
}

process RUN_MSA_NJ {
    tag "N=${taxa}_D=${dist}_L=${len}_rep=${rep}"

    input:
    tuple val(taxa), val(dist), val(len), val(rep), path(true_tree), path(fasta)

    output:
    path("msa_nj_N${taxa}_D${dist}_L${len}_rep${rep}.csv"), emit: csv
    path("msa_nj.nwk")
    path("msa.fasta")
    path("msa_matrix.phylip")

    script:
    """
    python3 ${projectDir}/../bin/run_msa_nj.py \\
        --fasta ${fasta} \\
        --outtree msa_nj.nwk \\
        --outmsa msa.fasta \\
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
        --outcsv msa_nj_N${taxa}_D${dist}_L${len}_rep${rep}.csv
    """
}

process RUN_MSA_ML {
    tag "N=${taxa}_D=${dist}_L=${len}_rep=${rep}"

    input:
    tuple val(taxa), val(dist), val(len), val(rep), path(true_tree), path(fasta)

    output:
    path("msa_ml_N${taxa}_D${dist}_L${len}_rep${rep}.csv"), emit: csv
    path("msa_ml.nwk")
    path("msa_ml_meta.json")

    script:
    """
    python3 ${projectDir}/../bin/run_msa_ml.py \\
        --fasta ${fasta} \\
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
        --outcsv msa_ml_N${taxa}_D${dist}_L${len}_rep${rep}.csv
    """
}

process COLLECT_AND_PLOT {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path csv_files

    output:
    path("benchmark_taxon_summary.csv")
    path("taxon_scaling_curves.png")
    path("benchmark_taxon_statistics.csv")

    script:
    """
    cat *.csv | awk 'NR==1 || \$0 !~ /^distance/' > benchmark_taxon_summary.csv
    python3 ${projectDir}/../bin/plot_taxon_benchmark.py --csv benchmark_taxon_summary.csv --outdir .
    """
}

workflow {
    ch_taxa       = Channel.fromList(params.taxa)
    ch_distances  = Channel.fromList(params.distances)
    ch_lengths    = Channel.fromList(params.lengths)
    ch_replicates = Channel.from(1..params.replicates)

    ch_params = ch_taxa.combine(ch_distances).combine(ch_lengths).combine(ch_replicates)

    SIMULATE_DATA(ch_params)

    ch_sim_data = SIMULATE_DATA.out.sim_data

    RUN_PWA_NJ(ch_sim_data)
    RUN_MSA_NJ(ch_sim_data)
    RUN_MSA_ML(ch_sim_data)

    ch_all_csvs = RUN_PWA_NJ.out.csv.mix(RUN_MSA_NJ.out.csv).mix(RUN_MSA_ML.out.csv).collect()

    COLLECT_AND_PLOT(ch_all_csvs)
}
