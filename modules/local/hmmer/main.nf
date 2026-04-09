process HMMSEARCH {
    tag "$sample_id.id"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmmer:3.4--h503566f_3' :
        'quay.io/biocontainers/hmmer:3.4--hb6cb901_4' }"

    publishDir "${params.outdir}/${sample_id.id}/hmmsearch", mode: 'copy', overwrite: true, pattern: "hmmsearch.*.{txt,sto,tbl,domtbl}.gz"

    input:
    tuple val(sample_id), path(hmm)
    tuple val(sample_id_b), path(orfs)

    output:
    tuple val(sample_id), path("hmmsearch.${sample_id.id}.txt.gz"),     emit: output
    tuple val(sample_id), path("hmmsearch.${sample_id.id}.sto.gz"),     emit: stockholm
    tuple val(sample_id), path("hmmsearch.${sample_id.id}.tbl.gz"),     emit: tbloutput
    tuple val(sample_id), path("hmmsearch.${sample_id.id}.domtbl.gz"),  emit: domtbloutput
    path 'versions.yml',                                                emit: versions

    script:
    arg_T       = params.hmmsearch_T    ? "-T ${params.hmmsearch_T}" : ''
    arg_domT    = params.hmmsearch_domT ? "--domT ${params.hmmsearch_domT}" : ''
    """
    TMP_ORFS=_tmp.${orfs.baseName}
    gunzip -c ${orfs} > \$TMP_ORFS

    hmmsearch \\
        $arg_T \\
        $arg_domT \\
        -o hmmsearch.${sample_id.id}.txt \\
        -A hmmsearch.${sample_id.id}.sto \\
        --tblout hmmsearch.${sample_id.id}.tbl \\
        --domtblout hmmsearch.${sample_id.id}.domtbl \\
        --notextw \\
        $hmm \\
        \$TMP_ORFS

    rm \$TMP_ORFS

    gzip --no-name *.txt *.sto *.tbl *.domtbl

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hmmer: \$(hmmsearch -h | grep -o '^# HMMER [0-9.]*' | sed 's/^# HMMER *//')
    END_VERSIONS
    """
}


