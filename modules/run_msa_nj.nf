nextflow.enable.dsl=2

process RUN_MSA_NJ {
    tag "D=${dist}_L=${len}_chk=${chunk_id}[${rep_start}..${rep_end}]"

    publishDir "${params.outdir}/replications", mode: 'copy',
        enabled: (params.containsKey('save_replications') ? params.save_replications : true),
        saveAs: { filename ->
            def m_tree = filename =~ /^msa_nj_(\d+)\.nwk$/
            if (m_tree) return "D${dist}_L${len}_rep${m_tree[0][1]}/msa_nj.nwk"
            def m_mat = filename =~ /^msa_matrix_(\d+)\.phylip$/
            if (m_mat) return "D${dist}_L${len}_rep${m_mat[0][1]}/msa_matrix.phylip"
            return null
        }

    input:
    tuple val(dist), val(len), val(chunk_id), val(rep_start), val(rep_end), path(true_trees), path(msas)

    output:
    path("chunk_msa_nj_D${dist}_L${len}_chk${chunk_id}.csv"), emit: csv
    path("msa_nj_*.nwk")
    path("msa_matrix_*.phylip")

    script:
    """
    for rep in \$(seq ${rep_start} ${rep_end}); do
        python3 ${moduleDir}/../bin/run_msa_nj.py \\
            --msa msa_\${rep}.fasta \\
            --outtree msa_nj_\${rep}.nwk \\
            --outmatrix msa_matrix_\${rep}.phylip \\
            --dist_model ${params.dist_model} \\
            --tool ${params.nj_tool}

        python3 ${moduleDir}/../bin/evaluate_trees.py \\
            --truetree true_tree_\${rep}.nwk \\
            --esttree msa_nj_\${rep}.nwk \\
            --matrix msa_matrix_\${rep}.phylip \\
            --pipeline MSA+NJ \\
            --distance ${dist} \\
            --length ${len} \\
            --alpha ${params.alpha} \\
            --replicate \${rep} \\
            --outcsv chunk_msa_nj_D${dist}_L${len}_chk${chunk_id}.csv
    done
    """
}
