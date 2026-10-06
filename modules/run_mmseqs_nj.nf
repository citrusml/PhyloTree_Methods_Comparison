nextflow.enable.dsl=2

process RUN_MMSEQS_NJ {
    tag "D=${dist}_L=${len}_chk=${chunk_id}[${rep_start}..${rep_end}]"

    publishDir "${params.outdir}/replications", mode: 'copy',
        enabled: (params.containsKey('save_replications') ? params.save_replications : true),
        saveAs: { filename ->
            def m_tree = filename =~ /^mmseqs_nj_(\d+)\.nwk$/
            if (m_tree) return "D${dist}_L${len}_rep${m_tree[0][1]}/mmseqs_nj.nwk"
            def m_mat = filename =~ /^mmseqs_matrix_(\d+)\.phylip$/
            if (m_mat) return "D${dist}_L${len}_rep${m_mat[0][1]}/mmseqs_matrix.phylip"
            def m_json = filename =~ /^mmseqs_meta_(\d+)\.json$/
            if (m_json) return "D${dist}_L${len}_rep${m_json[0][1]}/mmseqs_meta.json"
            return null
        }

    input:
    tuple val(dist), val(len), val(chunk_id), val(rep_start), val(rep_end), path(true_trees), path(fastas)

    output:
    path("chunk_mmseqs_nj_D${dist}_L${len}_chk${chunk_id}.csv"), emit: csv
    path("mmseqs_nj_*.nwk")
    path("mmseqs_matrix_*.phylip")
    path("mmseqs_meta_*.json")

    script:
    def nj_tool = params.containsKey('nj_tool') ? params.nj_tool : 'rapidnj'
    def sensitivity = params.containsKey('gs_sensitivity') ? params.gs_sensitivity : 7.5
    """
    for rep in \$(seq ${rep_start} ${rep_end}); do
        python3 ${moduleDir}/../bin/run_mmseqs_nj.py \\
            --fasta seqs_\${rep}.fasta \\
            --outtree mmseqs_nj_\${rep}.nwk \\
            --outmatrix mmseqs_matrix_\${rep}.phylip \\
            --outjson mmseqs_meta_\${rep}.json \\
            --dist_model log \\
            --tool ${nj_tool} \\
            --sensitivity ${sensitivity} \\
            --threads ${task.cpus}

        python3 ${moduleDir}/../bin/evaluate_trees.py \\
            --truetree true_tree_\${rep}.nwk \\
            --esttree mmseqs_nj_\${rep}.nwk \\
            --matrix mmseqs_matrix_\${rep}.phylip \\
            --pipeline MMSEQS+NJ \\
            --distance ${dist} \\
            --length ${len} \\
            --alpha ${params.alpha} \\
            --replicate \${rep} \\
            --json mmseqs_meta_\${rep}.json \\
            --outcsv chunk_mmseqs_nj_D${dist}_L${len}_chk${chunk_id}.csv
    done
    """
}
