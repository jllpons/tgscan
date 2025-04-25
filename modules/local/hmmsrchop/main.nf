process HMMSRCHOP_PARSE {
    tag "$meta.id"
    //label "process_multi"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.79' :
        'quay.io/biocontainers/biopython:1.79' }"

    publishDir "${params.outdir}/${meta.id}/hmmsrchop", mode: 'copy', overwrite: true, pattern: 'hmmsearch.out.serialized.json'


    input:
    tuple val(meta), path(hmmsearch_out)

    output:
    tuple val(meta), path('hmmsearch.out.serialized.json'), emit: hmmsearch_out_json
    path 'versions.yml',                                    emit: versions

    script:
    """
    hmmsrchop.py parse ${hmmsearch_out} > hmmsearch.out.serialized.json

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        Python: \$(python -V | awk '{print \$2}')
        Biopython: \$(python -c "import Bio; print(Bio.__version__)")
    END_VERSIONS
    """
}

process HMMSRCHOP_TOGFF {
    tag "$meta.id"
    //label "process_multi"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.79' :
        'quay.io/biocontainers/biopython:1.79' }"

    publishDir "${params.outdir}/${meta.id}/hmmsrchop", mode: 'copy', overwrite: true, pattern: 'hmmsearch.out.profile_alignments.gff3'


    input:
    tuple val(meta), path(hmmsearch_out_json)

    output:
    tuple val(meta), path('hmmsearch.out.profile_alignments.gff3'), emit: hmmsearch_out_togff
    path 'versions.yml',                                            emit: versions

    script:
    """
    hmmsrchop.py togff \
        --mode ${params.togff_mode} \
        --metric ${params.togff_metric} \
        ${hmmsearch_out_json} > hmmsearch.out.profile_alignments.gff3

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        Python: \$(python -V | awk '{print \$2}')
        Biopython: \$(python -c "import Bio; print(Bio.__version__)")
    END_VERSIONS
    """
}
