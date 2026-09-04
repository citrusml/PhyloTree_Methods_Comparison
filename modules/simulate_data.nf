nextflow.enable.dsl=2

process SIMULATE_DATA {
    tag "D=${dist}_L=${len}_chk=${chunk_id}[${rep_start}..${rep_end}]"

    input:
    tuple val(dist), val(len), val(chunk_id), val(rep_start), val(rep_end)

    output:
    tuple val(dist), val(len), val(chunk_id), val(rep_start), val(rep_end), path("true_tree_*.nwk"), path("seqs_*.fasta"), emit: sim_data
    tuple val(dist), val(len), val(chunk_id), val(rep_start), val(rep_end), path("true_tree_*.nwk"), path("true_msa_*.fasta"), emit: true_msa_data

    script:
    def model_str = params.alpha ? "${params.model}{${params.alpha}}" : "${params.model}"
    def indel_size_arg = params.indel_size ? "--indel-size \"${params.indel_size}\"" : ""
    def tree_model = params.containsKey('tree_model') ? params.tree_model.toString().toLowerCase() : 'birth_death'
    def rate_sd = params.containsKey('rate_sd') ? params.rate_sd : 0.0
    def lba_ratio = params.containsKey('lba_ratio') ? params.lba_ratio : 1.0
    def is_paper_yule = (tree_model == 'paper_yule' || tree_model == 'yule')
    """
    IQTREE_BIN=\$(which iqtree2 2>/dev/null || which iqtree3 2>/dev/null || which iqtree 2>/dev/null)
    if [ -z "\$IQTREE_BIN" ]; then
        echo "Error: IQ-TREE / AliSim not found in PATH" >&2
        exit 1
    fi
    for rep in \$(seq ${rep_start} ${rep_end}); do
        if [ "${is_paper_yule}" = "true" ]; then
            python3 ${moduleDir}/../bin/generate_tree.py \\
                --taxa ${params.taxa} \\
                --scale ${dist} \\
                --seed \${rep} \\
                --model paper_yule \\
                --rate_sd ${rate_sd} \\
                --lba_ratio ${lba_ratio} \\
                --outtree true_tree_\${rep}.nwk

            \${IQTREE_BIN} --alisim sim_\${rep} \\
                -m "${model_str}" \\
                --length ${len} \\
                -t true_tree_\${rep}.nwk \\
                --indel ${params.insert_rate},${params.delete_rate} \\
                ${indel_size_arg} \\
                -af fasta \\
                -seed \${rep} \\
                --redo
        else
            \${IQTREE_BIN} --alisim sim_\${rep} \\
                -m "${model_str}" \\
                --length ${len} \\
                -t "RANDOM{bd{${params.birth_rate}/${params.death_rate}}/${params.taxa}}" \\
                --indel ${params.insert_rate},${params.delete_rate} \\
                ${indel_size_arg} \\
                --branch-scale ${dist} \\
                -af fasta \\
                -seed \${rep} \\
                --redo

            if [ -f sim_\${rep}.treefile ]; then
                mv sim_\${rep}.treefile true_tree_\${rep}.nwk
            else
                mv sim_\${rep}.tree true_tree_\${rep}.nwk
            fi
        fi

        mv sim_\${rep}.unaligned.fa seqs_\${rep}.fasta
        mv sim_\${rep}.fa true_msa_\${rep}.fasta
    done
    """
}

