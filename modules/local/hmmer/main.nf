process HMMSEARCH {
    tag "$meta.id"
    //label "process_multi"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmmer:3.4--h503566f_3' :
        'quay.io/biocontainers/hmmer:3.4--h503566f_3' }"

    publishDir "${params.outdir}/${meta.id}/hmmsearch", mode: 'copy', overwrite: true//, pattern: 'hmmsearch.out *alignment.sto *domtblout'


    input:
    tuple val(meta), path(hmmfile), path(seqdb)

    output:
    path 'hmmsearch.out',  emit: hmmsearch_out
    path 'alignment.sto',    emit: alignment
    path 'domtblout',        emit: domtblout
    path 'versions.yml',     emit: versions

    script:
    // <https://www.nextflow.io/docs/latest/script.html#conditional-execution>
    def seq_eval     = params.seq_eval == null? "": "--incE ${params.seq_eval}"
    def seq_bitscore = params.seq_bitscore == null? "": "--incT ${params.seq_bitscore}"
    def dom_eval     = params.dom_eval == null? "": "--incdomE ${params.dom_eval}"
    def dom_bitscore = params.dom_bitscore == null? "": "--incdomT ${params.dom_bitscore}"
    def z_seq_evalue = params.z_seq_evalue == null? "": "-Z ${params.z_seq_evalue}"
    def z_dom_evalue = params.z_dom_evalue == null? "": "--domZ ${params.z_dom_evalue}"

    """
    hmmsearch \
        --domtblout domtblout \
        -A alignment.sto \
        ${seq_eval} \
        ${seq_bitscore} \
        ${dom_eval} \
        ${dom_bitscore} \
        ${z_seq_evalue} \
        ${z_dom_evalue} \
        $hmmfile \
        $seqdb > hmmsearch.out

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        hmmsearch: \$(hmmsearch -h | head -n 2 | tail -n 1 | awk '{print \$3}')
    END_VERSIONS
    """
}
