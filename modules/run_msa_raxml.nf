nextflow.enable.dsl=2

process RUN_MSA_RAXML {
    tag "D=${dist}_L=${len}_chk=${chunk_id}[${rep_start}..${rep_end}]"

    input:
    tuple val(dist), val(len), val(chunk_id), val(rep_start), val(rep_end), path(true_trees), path(msas)

    output:
    path("chunk_msa_raxml_D${dist}_L${len}_chk${chunk_id}.csv"), emit: csv
    path("msa_raxml_*.nwk")
    path("msa_raxml_meta_*.json")

    script:
    def raxml_model = params.containsKey('raxml_model') ? params.raxml_model : 'PROTGAMMALGX'
    """
    for rep in \$(seq ${rep_start} ${rep_end}); do
        python3 ${moduleDir}/../bin/run_msa_raxml.py \\
            --msa msa_\${rep}.fasta \\
            --outtree msa_raxml_\${rep}.nwk \\
            --outjson msa_raxml_meta_\${rep}.json \\
            --model "${raxml_model}" \\
            --seed \${rep} \\
            --threads ${task.cpus}

        python3 ${moduleDir}/../bin/evaluate_trees.py \\
            --truetree true_tree_\${rep}.nwk \\
            --esttree msa_raxml_\${rep}.nwk \\
            --pipeline MSA+RAXML \\
            --distance ${dist} \\
            --length ${len} \\
            --alpha ${params.alpha} \\
            --replicate \${rep} \\
            --json msa_raxml_meta_\${rep}.json \\
            --outcsv chunk_msa_raxml_D${dist}_L${len}_chk${chunk_id}.csv
    done
    """
}
