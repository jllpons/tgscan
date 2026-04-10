process HMMSEARCH_DOMTBLOUT_TO_GFF {
    tag "$sample_id.id"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:5.3.0' :
        'quay.io/biocontainers/gawk:5.3.0' }"

    publishDir "${params.outdir}/${sample_id.id}/hmmsearch", mode: 'copy', overwrite: true, pattern: "*.gff3.gz"

    input:
    tuple val(sample_id), path(hmmsearch_domtbl)

    output:
    tuple val(sample_id), path("*.gff3.gz"),    emit: hmmsearch_gff
    path 'versions.yml',                        emit: versions

    script:
    """
    tmp_domtbl_uncompressed="_${hmmsearch_domtbl.baseName}.domtbl"
    tmp_gff_uncompressed="_${hmmsearch_domtbl.baseName}.gff3"

    gzip -dc ${hmmsearch_domtbl} > \$tmp_domtbl_uncompressed

    hmmsearch_tblout2gff.awk \$tmp_domtbl_uncompressed > \$tmp_gff_uncompressed

    gzip -c \$tmp_gff_uncompressed > ${hmmsearch_domtbl.baseName}.gff3.gz

    rm -f \$tmp_domtbl_uncompressed \$tmp_gff_uncompressed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gawk: \$(gawk --version | gawk 'NR==1{print \$3}')
    END_VERSIONS
    """
}
