nextflow.enable.dsl=2

process RUN_MAFFT_LINSI {
    tag "D=${dist}_L=${len}_chk=${chunk_id}[${rep_start}..${rep_end}]"

    publishDir "${params.outdir}/replications", mode: 'copy',
        enabled: (params.containsKey('save_replications') ? params.save_replications : true),
        saveAs: { filename ->
            def m_msa = filename =~ /^msa_linsi_(\d+)\.fasta$/
            if (m_msa) return "D${dist}_L${len}_rep${m_msa[0][1]}/msa_linsi.fasta"
            return null
        }

    input:
    tuple val(dist), val(len), val(chunk_id), val(rep_start), val(rep_end), path(true_trees), path(fastas)

    output:
    tuple val(dist), val(len), val(chunk_id), val(rep_start), val(rep_end), path(true_trees), path("msa_linsi_*.fasta"), emit: msa_data

    script:
    """
    for rep in \$(seq ${rep_start} ${rep_end}); do
        mafft-linsi --thread ${task.cpus} --quiet seqs_\${rep}.fasta > msa_linsi_\${rep}.fasta
    done
    """
}
