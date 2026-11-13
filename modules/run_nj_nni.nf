nextflow.enable.dsl=2

process RUN_NJ_NNI {
    tag "${pipeline_name} D=${dist}_L=${len}_chk=${chunk_id}[${rep_start}..${rep_end}]"

    publishDir "${params.outdir}/replications", mode: 'copy',
        enabled: (params.containsKey('save_replications') ? params.save_replications : true),
        saveAs: { filename ->
            def m_tree = filename =~ /^(.+_nj_nni)_(\d+)\.nwk$/
            if (m_tree) return "D${dist}_L${len}_rep${m_tree[0][2]}/${m_tree[0][1]}.nwk"
            return null
        }

    input:
    tuple val(dist), val(len), val(chunk_id), val(rep_start), val(rep_end), path(true_trees), path(nj_trees), path(matrices), val(base_prefix), val(pipeline_name)

    output:
    path("chunk_${base_prefix}_nj_nni_D${dist}_L${len}_chk${chunk_id}.csv"), emit: csv
    path("${base_prefix}_nj_nni_*.nwk")

    script:
    def nni_type = params.containsKey('nj_nni_type') ? params.nj_nni_type : 'B'
    """
    for rep in \$(seq ${rep_start} ${rep_end}); do
        python3 ${moduleDir}/../bin/run_nj_nni.py \\
            --matrix ${base_prefix}_matrix_\${rep}.phylip \\
            --nj_tree ${base_prefix}_nj_\${rep}.nwk \\
            --outtree ${base_prefix}_nj_nni_\${rep}.nwk \\
            --nni_type ${nni_type}

        python3 ${moduleDir}/../bin/evaluate_trees.py \\
            --truetree true_tree_\${rep}.nwk \\
            --esttree ${base_prefix}_nj_nni_\${rep}.nwk \\
            --matrix ${base_prefix}_matrix_\${rep}.phylip \\
            --pipeline ${pipeline_name} \\
            --distance ${dist} \\
            --length ${len} \\
            --alpha ${params.alpha} \\
            --replicate \${rep} \\
            --outcsv chunk_${base_prefix}_nj_nni_D${dist}_L${len}_chk${chunk_id}.csv
    done
    """
}
