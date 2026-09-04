nextflow.enable.dsl=2

process CALCULATE_SSS {
    tag "D=${dist}_L=${len}_chk=${chunk_id}[${rep_start}..${rep_end}]"

    input:
    tuple val(dist), val(len), val(chunk_id), val(rep_start), val(rep_end), path(true_trees), path(fastas)

    output:
    path("chunk_sss_D${dist}_L${len}_chk${chunk_id}.csv"), emit: csv

    script:
    """
    python3 ${moduleDir}/../bin/calculate_sss.py \\
        --rep_start ${rep_start} \\
        --rep_end ${rep_end} \\
        --distance ${dist} \\
        --length ${len} \\
        --outcsv chunk_sss_D${dist}_L${len}_chk${chunk_id}.csv
    """
}
