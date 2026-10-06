nextflow.enable.dsl=2

process RUN_MSA_RAXML_LINSI {
    tag "D=${dist}_L=${len}_chk=${chunk_id}[${rep_start}..${rep_end}]"

    publishDir "${params.outdir}/replications", mode: 'copy',
        enabled: (params.containsKey('save_replications') ? params.save_replications : true),
        saveAs: { filename ->
            def m_tree = filename =~ /^msa_linsi_raxml_(\d+)\.nwk$/
            if (m_tree) return "D${dist}_L${len}_rep${m_tree[0][1]}/msa_linsi_raxml.nwk"
            def m_json = filename =~ /^msa_linsi_raxml_meta_(\d+)\.json$/
            if (m_json) return "D${dist}_L${len}_rep${m_json[0][1]}/msa_linsi_raxml_meta.json"
            return null
        }

    input:
    tuple val(dist), val(len), val(chunk_id), val(rep_start), val(rep_end), path(true_trees), path(msas)

    output:
    path("chunk_msa_linsi_raxml_D${dist}_L${len}_chk${chunk_id}.csv"), emit: csv
    path("msa_linsi_raxml_*.nwk")
    path("msa_linsi_raxml_meta_*.json")

    script:
    def raxml_model = params.containsKey('raxml_model') ? params.raxml_model : 'PROTGAMMAIWAGX'
    """
    for rep in \$(seq ${rep_start} ${rep_end}); do
        python3 ${moduleDir}/../bin/run_msa_raxml.py \\
            --msa msa_linsi_\${rep}.fasta \\
            --outtree msa_linsi_raxml_\${rep}.nwk \\
            --outjson msa_linsi_raxml_meta_\${rep}.json \\
            --model "${raxml_model}" \\
            --seed \${rep} \\
            --threads ${task.cpus}

        python3 ${moduleDir}/../bin/evaluate_trees.py \\
            --truetree true_tree_\${rep}.nwk \\
            --esttree msa_linsi_raxml_\${rep}.nwk \\
            --pipeline MSA_LINSI+RAXML \\
            --distance ${dist} \\
            --length ${len} \\
            --alpha ${params.alpha} \\
            --replicate \${rep} \\
            --json msa_linsi_raxml_meta_\${rep}.json \\
            --outcsv chunk_msa_linsi_raxml_D${dist}_L${len}_chk${chunk_id}.csv
    done
    """
}
