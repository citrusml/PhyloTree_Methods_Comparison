nextflow.enable.dsl=2

process RUN_MSA_ML {
    tag "D=${dist}_L=${len}_chk=${chunk_id}[${rep_start}..${rep_end}]"

    input:
    tuple val(dist), val(len), val(chunk_id), val(rep_start), val(rep_end), path(true_trees), path(msas)

    output:
    path("chunk_msa_ml_D${dist}_L${len}_chk${chunk_id}.csv"), emit: csv
    path("msa_ml_*.nwk")
    path("msa_ml_meta_*.json")

    script:
    """
    for rep in \$(seq ${rep_start} ${rep_end}); do
        python3 ${moduleDir}/../bin/run_msa_ml.py \\
            --msa msa_\${rep}.fasta \\
            --outtree msa_ml_\${rep}.nwk \\
            --outjson msa_ml_meta_\${rep}.json \\
            --threads ${task.cpus}

        python3 ${moduleDir}/../bin/evaluate_trees.py \\
            --truetree true_tree_\${rep}.nwk \\
            --esttree msa_ml_\${rep}.nwk \\
            --pipeline MSA+ML \\
            --distance ${dist} \\
            --length ${len} \\
            --alpha ${params.alpha} \\
            --replicate \${rep} \\
            --json msa_ml_meta_\${rep}.json \\
            --outcsv chunk_msa_ml_D${dist}_L${len}_chk${chunk_id}.csv
    done
    """
}
