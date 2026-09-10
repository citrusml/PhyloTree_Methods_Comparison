nextflow.enable.dsl=2

process RUN_PWA_NJ {
    tag "${key}"

    input:
    tuple val(key), path(seqs)

    output:
    tuple val(key), path("tree_pwa_nj.nwk"), emit: tree

    script:
    """
    python3 ${projectDir}/bin/run_pwa_nj.py \\
        --fasta ${seqs} \\
        --gap_open ${params.gap_open} \\
        --gap_extend ${params.gap_extend} \\
        --dist_model ${params.dist_model} \\
        --tool ${params.nj_tool} \\
        --outtree tree_pwa_nj.nwk
    """
}

process RUN_MAFFT {
    tag "${key}"

    input:
    tuple val(key), path(seqs)

    output:
    tuple val(key), path("mafft_aln.fasta"), emit: msa

    script:
    """
    MAFFT_BIN=\$(which mafft 2>/dev/null || echo "/opt/homebrew/bin/mafft")
    \${MAFFT_BIN} --auto --quiet ${seqs} > mafft_aln.fasta
    """
}

process RUN_MSA_NJ {
    tag "${key}"

    input:
    tuple val(key), path(mafft_msa)

    output:
    tuple val(key), path("tree_msa_nj.nwk"), emit: tree

    script:
    """
    python3 ${projectDir}/bin/run_msa_nj.py \\
        --msa ${mafft_msa} \\
        --dist_model ${params.dist_model} \\
        --tool ${params.nj_tool} \\
        --outtree tree_msa_nj.nwk
    """
}

process RUN_MSA_ML {
    tag "${key}"

    input:
    tuple val(key), path(mafft_msa)

    output:
    tuple val(key), path("tree_msa_ml.nwk"), emit: tree

    script:
    """
    IQTREE_BIN=\$(which iqtree2 2>/dev/null || which iqtree 2>/dev/null || echo "${projectDir}/bin/iqtree2")
    \${IQTREE_BIN} -s ${mafft_msa} -m ${params.model} --prefix iqtree_msa_ml -T 2 --quiet -redo
    mv iqtree_msa_ml.treefile tree_msa_ml.nwk
    """
}

process RUN_TRUE_MSA_NJ {
    tag "${key}"

    input:
    tuple val(key), path(true_msa)

    output:
    tuple val(key), path("tree_true_msa_nj.nwk"), emit: tree

    script:
    """
    python3 ${projectDir}/bin/run_msa_nj.py \\
        --msa ${true_msa} \\
        --dist_model ${params.dist_model} \\
        --tool ${params.nj_tool} \\
        --outtree tree_true_msa_nj.nwk
    """
}

process RUN_TRUE_MSA_ML {
    tag "${key}"

    input:
    tuple val(key), path(true_msa)

    output:
    tuple val(key), path("tree_true_msa_ml.nwk"), emit: tree

    script:
    """
    IQTREE_BIN=\$(which iqtree2 2>/dev/null || which iqtree 2>/dev/null || echo "${projectDir}/bin/iqtree2")
    \${IQTREE_BIN} -s ${true_msa} -m ${params.model} --prefix iqtree_true_ml -T 2 --quiet -redo
    mv iqtree_true_ml.treefile tree_true_msa_ml.nwk
    """
}

process RUN_SSS {
    tag "${key}"

    input:
    tuple val(key), val(dist), val(len), val(rep), path(seqs)

    output:
    tuple val(key), path("sss.csv"), emit: csv

    script:
    """
    python3 ${projectDir}/bin/calculate_sss.py \\
        --fasta ${seqs} \\
        --distance ${dist} \\
        --length ${len} \\
        --rep_start ${rep} \\
        --rep_end ${rep} \\
        --outcsv sss.csv
    """
}
