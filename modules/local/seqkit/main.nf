process SEQKIT_TRANSLATE {
    tag "$meta.id"
    //label "process_multi"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqkit:2.9.0--h9ee0642_0' :
        'quay.io/biocontainers/seqkit:2.9.0--h9ee0642_0' }"

    publishDir "${params.outdir}/${meta.id}/seqkit", mode: 'copy', overwrite: true, pattern: '*translated.fasta'


    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path('*.translated.fasta'),  emit: fasta_translated
    path 'versions.yml',                          emit: versions
    """
    seqkit translate \
        --append-frame \
        --out-subseqs \
        --frame 6 \
        --min-len ${params.translate_min_len} \
        --seq-type 'dna' \
        -o ${fasta.baseName}.translated.fasta \
        ${fasta}

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        seqkit: \$(seqkit version | sed -e "s/seqkit v//g")
    END_VERSIONS
    """
}
