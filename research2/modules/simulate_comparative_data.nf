nextflow.enable.dsl=2

process SIMULATE_COMPARATIVE_DATA {
    tag "${key}"

    input:
    tuple val(cond_name), val(simulator), val(ins_rate), val(del_rate), val(dist), val(len), val(rep)

    output:
    tuple val(key), val(cond_name), val(simulator), val(ins_rate), val(del_rate), val(dist), val(len), val(rep), path("true_tree_${rep}.nwk"), path("seqs_${rep}.fasta"), path("true_msa_${rep}.fasta"), emit: sim_data

    script:
    key = "${cond_name}_D${dist}_L${len}_rep${rep}"
    def cur_seed = 100000 + rep * 1000 + (dist * 100 as int)
    def indel_size_str = params.containsKey('indel_size') ? params.indel_size.toString() : "POW 1.7 50"
    def indel_arg = "--indel-size \"${indel_size_str}\""

    // INDELible specific parsing
    def raw_model = params.model ? params.model.toString().toUpperCase() : 'LG'
    def submodel_str = raw_model.contains('LG') ? 'LG'
                     : raw_model.contains('WAG') ? 'WAG'
                     : raw_model.contains('JTT') ? 'JTT'
                     : raw_model.contains('DAYHOFF') ? 'DAYHOFF'
                     : raw_model.contains('BLOSUM62') ? 'BLOSUM62'
                     : 'LG'
    def spec = indel_size_str.replaceAll(/[\{\},]/, ' ').replaceAll(/\//, ' ').replaceAll(/\s+/, ' ').trim()
    def parts = spec.split('\\s+')
    def indel_type = parts.length > 0 ? parts[0] : 'POW'
    def indel_p1 = parts.length > 1 ? parts[1] : '1.7'
    def indel_p2 = parts.length > 2 ? parts[2] : '50'
    def indelible_cmd = "${indel_type} ${indel_p1} ${indel_p2}"
    def alisim_indel_arg = "--indel-size \"${indel_type}{${indel_p1}/${indel_p2}},${indel_type}{${indel_p1}/${indel_p2}}\""

    """
    if [ "${simulator}" = "indelible" ]; then
        INDELIBLE_BIN=\$(which indelible 2>/dev/null || echo "${projectDir}/bin/indelible")
        if [ ! -x "\${INDELIBLE_BIN}" ]; then
            if [ -f "${projectDir}/bin/install_indelible.sh" ]; then
                bash "${projectDir}/bin/install_indelible.sh" >&2
            fi
        fi
        if [ ! -x "\${INDELIBLE_BIN}" ]; then
            echo "Error: INDELible binary not found at \${INDELIBLE_BIN}" >&2
            exit 1
        fi
    else
        IQTREE_BIN=\$(which iqtree2 2>/dev/null || which iqtree 2>/dev/null || echo "${projectDir}/bin/iqtree2")
        if [ -z "\$IQTREE_BIN" ] || [ ! -x "\$IQTREE_BIN" ]; then
            echo "Error: IQ-TREE / AliSim not found at \${IQTREE_BIN}" >&2
            exit 1
        fi
    fi

    sim_ok=false
    for attempt in \$(seq 0 4); do
        run_seed=\$(( ${cur_seed} + attempt * 10000 ))

        # 1. 系統樹生成 (paper_yule)
        python3 ${projectDir}/bin/generate_tree.py \\
            --taxa ${params.taxa} \\
            --scale ${dist} \\
            --seed \${run_seed} \\
            --model paper_yule \\
            --rate_sd 0.0 \\
            --lba_ratio 1.0 \\
            --outtree true_tree_${rep}.nwk

        # 2. シミュレーション実行
        if [ "${simulator}" = "indelible" ]; then
            TREE_STR=\$(cat true_tree_${rep}.nwk | tr -d '\\n\\r ')
            if [[ ! "\${TREE_STR}" =~ \\;\$ ]]; then
                TREE_STR="\${TREE_STR};"
            fi

            WORK_IND="work_indelible_${rep}_\${attempt}"
            rm -rf "\${WORK_IND}" && mkdir -p "\${WORK_IND}"

            cat << EOF > "\${WORK_IND}/control.txt"
[TYPE] AMINOACID 1
[SETTINGS]
  [randomseed] \${run_seed}
  [output] FASTA
[MODEL] mymodel
  [submodel] ${submodel_str}
  [rates] 0 ${params.alpha} 4
  [insertmodel] ${indelible_cmd}
  [deletemodel] ${indelible_cmd}
  [insertrate] ${ins_rate}
  [deleterate] ${del_rate}
[TREE] mytree \${TREE_STR}
[PARTITIONS] mypart [mytree mymodel ${len}]
[EVOLVE] mypart 1 sim_${rep}
EOF

            if (cd "\${WORK_IND}" && "\${INDELIBLE_BIN}" control.txt > indelible.log 2>&1); then
                if [ -f "\${WORK_IND}/sim_${rep}.fas" ] && [ -f "\${WORK_IND}/sim_${rep}_TRUE.fas" ]; then
                    if python3 -c '
import sys
fname = sys.argv[1]
taxa = int(sys.argv[2])
count, valid, seq_len, in_header = 0, True, 0, False
with open(fname) as f:
    for line in f:
        line = line.strip()
        if not line: continue
        if line.startswith(">"):
            if in_header and seq_len == 0:
                valid = False; break
            count += 1; in_header = True; seq_len = 0
        else:
            seq_len += len(line)
    if in_header and seq_len == 0:
        valid = False
sys.exit(0 if (valid and count == taxa) else 1)
' "\${WORK_IND}/sim_${rep}.fas" "${params.taxa}"; then
                        mv "\${WORK_IND}/sim_${rep}.fas" seqs_${rep}.fasta
                        mv "\${WORK_IND}/sim_${rep}_TRUE.fas" true_msa_${rep}.fasta
                        rm -rf "\${WORK_IND}"
                        sim_ok=true
                        break
                    fi
                fi
            fi
            rm -rf "\${WORK_IND}"
        else
            # AliSim 実行
            if \${IQTREE_BIN} --alisim sim_${rep} \\
                -m "${params.model}" \\
                --length ${len} \\
                -t true_tree_${rep}.nwk \\
                --indel ${ins_rate},${del_rate} \\
                ${alisim_indel_arg} \\
                -af fasta \\
                -seed \${run_seed} \\
                --redo > alisim_${rep}.log 2>&1; then
                if [ -f "sim_${rep}.unaligned.fa" ] && [ -f "sim_${rep}.fa" ]; then
                    mv sim_${rep}.unaligned.fa seqs_${rep}.fasta
                    mv sim_${rep}.fa true_msa_${rep}.fasta
                    rm -f alisim_${rep}.log
                    sim_ok=true
                    break
                fi
            fi
        fi
    done

    if [ "\${sim_ok}" != "true" ]; then
        echo "Error: Simulation failed for ${cond_name} rep ${rep}" >&2
        exit 2
    fi
    """
}
