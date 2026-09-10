nextflow.enable.dsl=2

process SIMULATE_DATA {
    tag "D=${dist}_L=${len}_chk=${chunk_id}[${rep_start}..${rep_end}]"

    input:
    tuple val(dist), val(len), val(chunk_id), val(rep_start), val(rep_end)

    output:
    tuple val(dist), val(len), val(chunk_id), val(rep_start), val(rep_end), path("true_tree_*.nwk"), path("seqs_*.fasta"), emit: sim_data
    tuple val(dist), val(len), val(chunk_id), val(rep_start), val(rep_end), path("true_tree_*.nwk"), path("true_msa_*.fasta"), emit: true_msa_data

    script:
    def simulator = params.containsKey('simulator') ? params.simulator.toString().toLowerCase() : 'alisim'
    def model_str = params.alpha ? "${params.model}{${params.alpha}}" : "${params.model}"
    def indel_size_arg = params.indel_size ? "--indel-size \"${params.indel_size}\"" : ""
    def tree_model = params.containsKey('tree_model') ? params.tree_model.toString().toLowerCase() : 'birth_death'
    def rate_sd = params.containsKey('rate_sd') ? params.rate_sd : 0.0
    def lba_ratio = params.containsKey('lba_ratio') ? params.lba_ratio : 1.0
    def is_paper_yule = (tree_model == 'paper_yule' || tree_model == 'yule')

    // INDELible specific parameter parsing
    def raw_model = params.model ? params.model.toString().toUpperCase() : 'WAG'
    def submodel_str = raw_model.contains('WAG') ? 'WAG'
                     : raw_model.contains('JTT') ? 'JTT'
                     : raw_model.contains('DAYHOFF') ? 'DAYHOFF'
                     : raw_model.contains('BLOSUM62') ? 'BLOSUM62'
                     : 'WAG'
    def indel_str = params.containsKey('indel_size') ? params.indel_size.toString() : "POW 1.7 50"
    def indelible_spec = indel_str.replaceAll(/[\{\},]/, ' ').replaceAll(/\//, ' ').replaceAll(/\s+/, ' ').trim()
    def m_parts = indelible_spec.split('\\s+')
    def indel_type = m_parts.length > 0 ? m_parts[0] : 'POW'
    def indel_p1 = m_parts.length > 1 ? m_parts[1] : '1.7'
    def indel_p2 = m_parts.length > 2 ? m_parts[2] : '50'
    def indelible_indel_cmd = "${indel_type} ${indel_p1} ${indel_p2}"

    """
    if [ "${simulator}" = "indelible" ]; then
        INDELIBLE_BIN=\$(which indelible 2>/dev/null || echo "${moduleDir}/../bin/indelible")
        if [ ! -x "\${INDELIBLE_BIN}" ]; then
            if [ -f "${moduleDir}/../bin/install_indelible.sh" ]; then
                echo "Compiling INDELible binary via install_indelible.sh..." >&2
                bash "${moduleDir}/../bin/install_indelible.sh" >&2
            fi
        fi
        if [ ! -x "\${INDELIBLE_BIN}" ]; then
            echo "Error: INDELible binary not found or not executable at \${INDELIBLE_BIN}" >&2
            exit 1
        fi
    else
        IQTREE_BIN=\$(which iqtree2 2>/dev/null || which iqtree3 2>/dev/null || which iqtree 2>/dev/null || echo "${moduleDir}/../bin/iqtree2")
        if [ -z "\$IQTREE_BIN" ] || [ ! -x "\$IQTREE_BIN" ]; then
            echo "Error: IQ-TREE / AliSim not found in PATH or bin/" >&2
            exit 1
        fi
    fi

    for rep in \$(seq ${rep_start} ${rep_end}); do
        sim_ok=false
        for attempt in \$(seq 0 4); do
            cur_seed=\$((rep + attempt * 100000))

            # 1. 系統樹生成
            if [ "${is_paper_yule}" = "true" ]; then
                python3 ${moduleDir}/../bin/generate_tree.py \\
                    --taxa ${params.taxa} \\
                    --scale ${dist} \\
                    --seed \${cur_seed} \\
                    --model paper_yule \\
                    --rate_sd ${rate_sd} \\
                    --lba_ratio ${lba_ratio} \\
                    --outtree true_tree_\${rep}.nwk
            else
                python3 ${moduleDir}/../bin/generate_tree.py \\
                    --taxa ${params.taxa} \\
                    --scale ${dist} \\
                    --seed \${cur_seed} \\
                    --model birth_death \\
                    --birth_rate ${params.birth_rate} \\
                    --death_rate ${params.death_rate} \\
                    --outtree true_tree_\${rep}.nwk 2>/dev/null || \\
                python3 ${moduleDir}/../bin/generate_tree.py \\
                    --taxa ${params.taxa} \\
                    --scale ${dist} \\
                    --seed \${cur_seed} \\
                    --model paper_yule \\
                    --outtree true_tree_\${rep}.nwk
            fi

            # 2. シミュレーション実行 (INDELible vs AliSim)
            if [ "${simulator}" = "indelible" ]; then
                TREE_STR=\$(cat true_tree_\${rep}.nwk | tr -d '\\n\\r ')
                if [[ ! "\${TREE_STR}" =~ \\;\$ ]]; then
                    TREE_STR="\${TREE_STR};"
                fi

                WORK_IND="work_indelible_\${rep}_\${attempt}"
                rm -rf "\${WORK_IND}" && mkdir -p "\${WORK_IND}"

                cat << EOF > "\${WORK_IND}/control.txt"
[TYPE] AMINOACID 1
[SETTINGS]
  [randomseed] \${cur_seed}
  [output] FASTA
[MODEL] mymodel
  [submodel] ${submodel_str}
  [rates] 0 ${params.alpha} 4
  [insertmodel] ${indelible_indel_cmd}
  [deletemodel] ${indelible_indel_cmd}
  [insertrate] ${params.insert_rate}
  [deleterate] ${params.delete_rate}
[TREE] mytree \${TREE_STR}
[PARTITIONS] mypart [mytree mymodel ${len}]
[EVOLVE] mypart 1 sim_\${rep}
EOF

                if (cd "\${WORK_IND}" && "\${INDELIBLE_BIN}" control.txt > indelible.log 2>&1); then
                    if [ -f "\${WORK_IND}/sim_\${rep}.fas" ] && [ -f "\${WORK_IND}/sim_\${rep}_TRUE.fas" ]; then
                        if python3 -c '
import sys
fname = sys.argv[1]
taxa = int(sys.argv[2])
count, valid, seq_len, in_header = 0, True, 0, False
with open(fname) as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if in_header and seq_len == 0:
                valid = False; break
            count += 1; in_header = True; seq_len = 0
        else:
            seq_len += len(line)
    if in_header and seq_len == 0:
        valid = False
sys.exit(0 if (valid and count == taxa) else 1)
' "\${WORK_IND}/sim_\${rep}.fas" "${params.taxa}"; then
                            mv "\${WORK_IND}/sim_\${rep}.fas" seqs_\${rep}.fasta
                            mv "\${WORK_IND}/sim_\${rep}_TRUE.fas" true_msa_\${rep}.fasta
                            rm -rf "\${WORK_IND}"
                            sim_ok=true
                            break
                        fi
                    fi
                fi
                if [ -f "\${WORK_IND}/indelible.log" ]; then
                    echo "--- INDELible Log (attempt \${attempt}) ---" >&2
                    cat "\${WORK_IND}/indelible.log" >&2
                fi
                rm -rf "\${WORK_IND}"
            else
                # AliSim 実行ブロック
                if [ "${is_paper_yule}" = "true" ]; then
                    if \${IQTREE_BIN} --alisim sim_\${rep} \\
                        -m "${model_str}" \\
                        --length ${len} \\
                        -t true_tree_\${rep}.nwk \\
                        --indel ${params.insert_rate},${params.delete_rate} \\
                        ${indel_size_arg} \\
                        -af fasta \\
                        -seed \${cur_seed} \\
                        --redo > alisim_\${rep}.log 2>&1; then
                        if [ -f "sim_\${rep}.unaligned.fa" ] && [ -f "sim_\${rep}.fa" ]; then
                            if python3 -c '
import sys
fname = sys.argv[1]
taxa = int(sys.argv[2])
count, valid, seq_len, in_header = 0, True, 0, False
with open(fname) as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if in_header and seq_len == 0:
                valid = False; break
            count += 1; in_header = True; seq_len = 0
        else:
            seq_len += len(line)
    if in_header and seq_len == 0:
        valid = False
sys.exit(0 if (valid and count == taxa) else 1)
' "sim_\${rep}.unaligned.fa" "${params.taxa}"; then
                                mv sim_\${rep}.unaligned.fa seqs_\${rep}.fasta
                                mv sim_\${rep}.fa true_msa_\${rep}.fasta
                                rm -f alisim_\${rep}.log
                                sim_ok=true
                                break
                            fi
                        fi
                    fi
                else
                    if \${IQTREE_BIN} --alisim sim_\${rep} \\
                        -m "${model_str}" \\
                        --length ${len} \\
                        -t "RANDOM{bd{${params.birth_rate}/${params.death_rate}}/${params.taxa}}" \\
                        --indel ${params.insert_rate},${params.delete_rate} \\
                        ${indel_size_arg} \\
                        --branch-scale ${dist} \\
                        -af fasta \\
                        -seed \${cur_seed} \\
                        --redo > alisim_\${rep}.log 2>&1; then
                        if [ -f "sim_\${rep}.unaligned.fa" ] && [ -f "sim_\${rep}.fa" ]; then
                            if python3 -c '
import sys
fname = sys.argv[1]
taxa = int(sys.argv[2])
count, valid, seq_len, in_header = 0, True, 0, False
with open(fname) as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if in_header and seq_len == 0:
                valid = False; break
            count += 1; in_header = True; seq_len = 0
        else:
            seq_len += len(line)
    if in_header and seq_len == 0:
        valid = False
sys.exit(0 if (valid and count == taxa) else 1)
' "sim_\${rep}.unaligned.fa" "${params.taxa}"; then
                                if [ -f sim_\${rep}.treefile ]; then
                                    mv sim_\${rep}.treefile true_tree_\${rep}.nwk
                                elif [ -f sim_\${rep}.tree ]; then
                                    mv sim_\${rep}.tree true_tree_\${rep}.nwk
                                fi
                                mv sim_\${rep}.unaligned.fa seqs_\${rep}.fasta
                                mv sim_\${rep}.fa true_msa_\${rep}.fasta
                                rm -f alisim_\${rep}.log
                                sim_ok=true
                                break
                            fi
                        fi
                    fi
                fi
            fi
        done

        if [ "\${sim_ok}" != "true" ]; then
            echo "Error: Simulation (${simulator}) failed for replicate \${rep} after 5 attempts." >&2
            exit 2
        fi
    done
    """
}
