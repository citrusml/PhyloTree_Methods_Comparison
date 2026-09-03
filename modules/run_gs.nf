nextflow.enable.dsl=2

process RUN_GS {
    tag "D=${dist}_L=${len}_chk=${chunk_id}[${rep_start}..${rep_end}]"

    input:
    tuple val(dist), val(len), val(chunk_id), val(rep_start), val(rep_end), path(true_trees), path(fastas)

    output:
    path("chunk_gs_D${dist}_L${len}_chk${chunk_id}.csv"), emit: csv
    path("gs_*.nwk")

    script:
    def sens_val = params.gs_sensitivity ?: 7.5
    """
    GS_BIN=\$(which gs2 2>/dev/null || echo "${moduleDir}/../bin/gs2")
    if [ ! -x "\$GS_BIN" ]; then
        echo "Error: gs2 binary not found in PATH or ${moduleDir}/../bin/gs2" >&2
        exit 1
    fi

    for rep in \$(seq ${rep_start} ${rep_end}); do
        \${GS_BIN} -s -l -t ${task.cpus} -m ${sens_val} seqs_\${rep}.fasta | tr -d '"' > gs_\${rep}.nwk

        python3 ${moduleDir}/../bin/evaluate_trees.py \\
            --truetree true_tree_\${rep}.nwk \\
            --esttree gs_\${rep}.nwk \\
            --pipeline GS \\
            --distance ${dist} \\
            --length ${len} \\
            --alpha ${params.alpha} \\
            --replicate \${rep} \\
            --outcsv chunk_gs_D${dist}_L${len}_chk${chunk_id}.csv
    done
    """
}
