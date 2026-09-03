nextflow.enable.dsl=2

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
        python3 ${moduleDir}/../bin/run_pwa_nj.py \\
            --fasta seqs_\${rep}.fasta \\
            --outtree pwa_nj_\${rep}.nwk \\
            --outmatrix pwa_matrix_\${rep}.phylip \\
            --gap_open ${params.gap_open} \\
            --gap_extend ${params.gap_extend} \\
            --dist_model ${params.dist_model} \\
            --tool ${params.nj_tool} \\
            --threads ${task.cpus}

        python3 ${moduleDir}/../bin/evaluate_trees.py \\
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
