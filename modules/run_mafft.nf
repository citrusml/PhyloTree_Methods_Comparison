nextflow.enable.dsl=2

process RUN_MAFFT {
    tag "D=${dist}_L=${len}_chk=${chunk_id}[${rep_start}..${rep_end}]"

    publishDir "${params.outdir}/replications", mode: 'copy',
        enabled: (params.containsKey('save_replications') ? params.save_replications : true),
        saveAs: { filename ->
            def m_msa = filename =~ /^msa_(\d+)\.fasta$/
            if (m_msa) return "D${dist}_L${len}_rep${m_msa[0][1]}/msa.fasta"
            return null
        }

    input:
    tuple val(dist), val(len), val(chunk_id), val(rep_start), val(rep_end), path(true_trees), path(fastas)

    output:
    tuple val(dist), val(len), val(chunk_id), val(rep_start), val(rep_end), path(true_trees), path("msa_*.fasta"), emit: msa_data

    script:
    """
    for rep in \$(seq ${rep_start} ${rep_end}); do
        mafft --threadit 0 --auto --thread ${task.cpus} --quiet seqs_\${rep}.fasta > msa_\${rep}.fasta
    done
    """
}
