nextflow.enable.dsl=2

process RUN_MSA_BI {
    tag "D=${dist}_L=${len}_chk=${chunk_id}[${rep_start}..${rep_end}]"

    input:
    tuple val(dist), val(len), val(chunk_id), val(rep_start), val(rep_end), path(true_trees), path(msas)

    output:
    path("chunk_msa_bi_D${dist}_L${len}_chk${chunk_id}.csv"), emit: csv
    path("msa_bi_*.nwk")

    script:
    def ngen_val = params.bi_ngen ?: 100000
    def samplefreq_val = params.bi_samplefreq ?: 1000
    def burnin_val = params.bi_burnin ?: 20
    def nchains_val = params.bi_nchains ?: 4
    def model_val = params.bi_model ?: (params.model ? params.model.split('\\+')[0].toLowerCase() : "lg")
    """
    for rep in \$(seq ${rep_start} ${rep_end}); do
        python3 ${moduleDir}/../bin/run_msa_bi.py \\
            --msa msa_\${rep}.fasta \\
            --outtree msa_bi_\${rep}.nwk \\
            --ngen ${ngen_val} \\
            --samplefreq ${samplefreq_val} \\
            --burnin ${burnin_val} \\
            --nchains ${nchains_val} \\
            --model ${model_val} \\
            --threads ${task.cpus}

        python3 ${moduleDir}/../bin/evaluate_trees.py \\
            --truetree true_tree_\${rep}.nwk \\
            --esttree msa_bi_\${rep}.nwk \\
            --pipeline MSA+BI \\
            --distance ${dist} \\
            --length ${len} \\
            --alpha ${params.alpha} \\
            --replicate \${rep} \\
            --outcsv chunk_msa_bi_D${dist}_L${len}_chk${chunk_id}.csv
    done
    """
}
