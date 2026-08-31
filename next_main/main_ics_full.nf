nextflow.enable.dsl=2

/*
 * Full-length Invariant Category Sites (ICS Full) Benchmark Pipeline
 * Simulates sequences under 100% ICS model with Dayhoff 6 classes and realistic Indels.
 * Compares PWA+NJ, MSA+NJ, and MSA+ML pipelines across Distance D x Sequence Length L.
 */
params.taxa        = 32
params.distances   = [0.1, 0.5, 1.0, 2.0, 3.0]
params.lengths     = [100, 300, 500, 1000, 1500]
params.ics_prop    = 1.0                         // 100% full-length ICS
params.replicates  = 100
params.chunk_size  = 10
params.birth_rate  = 0.1
params.death_rate  = 0.05
params.insert_rate = 0.05
params.delete_rate = 0.10
params.model       = "LG+G4"
params.alpha       = 1.0                         // Fixed alpha for ICS+G4
params.dist_model  = "poisson"                   // Distance formula
params.gap_open    = 10.0
params.gap_extend  = 0.5
params.outdir      = "results/results_ics_full"
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
            --ics_prop ${params.ics_prop} \\
            --ics_model_file ${projectDir}/../models/ics_model.nex \\
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
    path("chunk_pwa_nj_ics_full_D${dist}_L${len}_chk${chunk_id}.csv"), emit: csv
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
            --ics_prop ${params.ics_prop} \\
            --replicate \${rep} \\
            --outcsv chunk_pwa_nj_ics_full_D${dist}_L${len}_chk${chunk_id}.csv
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
    path("chunk_msa_nj_ics_full_D${dist}_L${len}_chk${chunk_id}.csv"), emit: csv
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
            --ics_prop ${params.ics_prop} \\
            --replicate \${rep} \\
            --outcsv chunk_msa_nj_ics_full_D${dist}_L${len}_chk${chunk_id}.csv
    done
    """
}

process RUN_MSA_ML {
    tag "D=${dist}_L=${len}_chk=${chunk_id}[${rep_start}..${rep_end}]"

    input:
    tuple val(dist), val(len), val(chunk_id), val(rep_start), val(rep_end), path(true_trees), path(msas)

    output:
    path("chunk_msa_ml_ics_full_D${dist}_L${len}_chk${chunk_id}.csv"), emit: csv
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
            --ics_prop ${params.ics_prop} \\
            --replicate \${rep} \\
            --json msa_ml_meta_\${rep}.json \\
            --outcsv chunk_msa_ml_ics_full_D${dist}_L${len}_chk${chunk_id}.csv
    done
    """
}

process COLLECT_AND_PLOT {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path csv_files

    output:
    path("benchmark_ics_full_summary.csv")
    path("ics_full_benchmark.png")
    path("ics_full_benchmark.pdf"), optional: true
    path("benchmark_ics_full_statistics.csv")

    script:
    """
    # Concatenate CSVs: keep header from first file only
    head -n 1 \$(ls *.csv | head -1) > benchmark_ics_full_summary.csv
    for f in *.csv; do tail -n +2 "\$f"; done >> benchmark_ics_full_summary.csv
    python3 ${projectDir}/../bin/plot/plot_ics_full_benchmark.py --csv benchmark_ics_full_summary.csv --outdir .
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
