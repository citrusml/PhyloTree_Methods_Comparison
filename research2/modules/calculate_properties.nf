nextflow.enable.dsl=2

process CALCULATE_PROPERTIES {
    tag "${key}"

    input:
    tuple val(key),
          val(cond_name), val(simulator), val(ins_rate), val(del_rate), val(dist), val(len), val(rep),
          path(true_tree), path(seqs), path(true_msa),
          path(mafft_msa),
          path(tree_pwa_nj),
          path(tree_msa_nj),
          path(tree_msa_ml),
          path(tree_true_nj),
          path(tree_true_ml),
          path(sss_csv)

    output:
    path("props_${key}.csv"), emit: csv

    script:
    """
    python3 ${moduleDir}/../scripts/calculate_msa_properties.py \\
        --cond_name ${cond_name} \\
        --simulator ${simulator} \\
        --ins_rate ${ins_rate} \\
        --del_rate ${del_rate} \\
        --dist ${dist} \\
        --length ${len} \\
        --rep ${rep} \\
        --true_tree ${true_tree} \\
        --true_msa ${true_msa} \\
        --unaligned ${seqs} \\
        --mafft_msa ${mafft_msa} \\
        --tree_pwa_nj ${tree_pwa_nj} \\
        --tree_msa_nj ${tree_msa_nj} \\
        --tree_msa_ml ${tree_msa_ml} \\
        --tree_true_msa_nj ${tree_true_nj} \\
        --tree_true_msa_ml ${tree_true_ml} \\
        --sss_csv ${sss_csv} \\
        --outfile props_${key}.csv
    """
}
