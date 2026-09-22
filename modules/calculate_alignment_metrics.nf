nextflow.enable.dsl=2

process CALCULATE_ALIGNMENT_METRICS {
    tag "D=${dist}_L=${len}_chk=${chunk_id}[${rep_start}..${rep_end}]"

    publishDir "${params.outdir}/replications", mode: 'copy',
        enabled: (params.containsKey('save_replications') ? params.save_replications : true),
        saveAs: { filename ->
            def m_mat = filename =~ /^true_matrix_(\d+)\.phylip$/
            if (m_mat) return "D${dist}_L${len}_rep${m_mat[0][1]}/true_matrix.phylip"
            return null
        }

    input:
    tuple val(dist), val(len), val(chunk_id), val(rep_start), val(rep_end), path(true_trees), path(true_msas), path(est_msas)

    output:
    path("chunk_alignment_D${dist}_L${len}_chk${chunk_id}.csv"), emit: csv
    path("true_matrix_*.phylip"), optional: true

    script:
    """
    for rep in \$(seq ${rep_start} ${rep_end}); do
        python3 ${moduleDir}/../bin/calculate_alignment_metrics.py \\
            --truemsa true_msa_\${rep}.fasta \\
            --estmsa msa_\${rep}.fasta \\
            --truetree true_tree_\${rep}.nwk \\
            --outmatrix true_matrix_\${rep}.phylip \\
            --distance ${dist} \\
            --length ${len} \\
            --replicate \${rep} \\
            --outcsv chunk_alignment_D${dist}_L${len}_chk${chunk_id}.csv
    done
    """
}
