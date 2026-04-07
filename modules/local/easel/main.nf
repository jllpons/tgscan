process TRANSLATE {
    tag "$sample_id.id"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/easel:0.49--h7b50bb2_1' :
        'quay.io/biocontainers/easel:0.49--hb6cb901_3' }"

    publishDir "${params.outdir}/${sample_id.id}/translate", mode: 'copy', overwrite: true, pattern: "*.6frameORFs.fasta.gz"

    input:
    tuple val(sample_id), path(genome)

    output:
    tuple val(sample_id), path("*.6frameORFs.fasta.gz"),    emit: orfs
    path 'versions.yml',                                    emit: versions

    script:
    """
    esl-translate "${genome}" | gzip --to-stdout > "${genome.baseName}.6frameORFs.fasta.gz"


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        easl-translate: \$(esl-tranlsate -h | grep -o '^# Easel [0-9.]*' | sed 's/^# Easel *//')
    END_VERSIONS
    """
}


