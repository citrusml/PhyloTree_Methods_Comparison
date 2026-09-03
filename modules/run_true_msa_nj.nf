nextflow.enable.dsl=2

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
        python3 ${moduleDir}/../bin/run_msa_nj.py \\
            --msa true_msa_\${rep}.fasta \\
            --outtree true_msa_nj_\${rep}.nwk \\
            --outmatrix true_msa_matrix_\${rep}.phylip \\
            --dist_model ${params.dist_model} \\
            --tool ${params.nj_tool}

        python3 ${moduleDir}/../bin/evaluate_trees.py \\
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
