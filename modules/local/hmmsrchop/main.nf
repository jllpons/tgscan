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

process HMMSRCHOP_PLOT {
    tag "$meta.id"
    //label "process_multi"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.79' :
        'quay.io/biocontainers/biopython:1.79' }"

    publishDir "${params.outdir}/${meta.id}/hmmsrchop/plots", mode: 'copy', overwrite: true, pattern: '*.jpeg'


    input:
    tuple val(meta), path(hmmsearch_out_json)

    output:
    tuple val(meta), path('*.jpeg'),    emit: plots, optional: true
    path 'versions.yml',                emit: versions

    script:
    def seq_eval_threshold_arg = params.seq_eval == null ? "" : "--seq-eval-threshold ${params.seq_eval}"
    """
    hmmsrchop.py plot \
        ${seq_eval_threshold_arg} \
        ${hmmsearch_out_json}

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        Python: \$(python -V | awk '{print \$2}')
        Biopython: \$(python -c "import Bio; print(Bio.__version__)")
    END_VERSIONS
    """
}
