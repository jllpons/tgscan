include { SEQKIT_TRANSLATE } from '../modules/local/seqkit/main'
include { HMMSEARCH        } from '../modules/local/hmmer/main'

workflow TGSCAN {

    take:
    samplesheet // channel: samplesheet read in from --input

    main:

    ch_versions = Channel.empty()

    ch_fasta = samplesheet
        .map { meta, fasta, _hmmfile ->
            [meta, fasta]
        }
    ch_hmmfile = samplesheet
        .map { meta, _fasta, hmmfile ->
            [meta, hmmfile]
        }

    SEQKIT_TRANSLATE(
        ch_fasta,
    )

    ch_hmmfile_fastaTranslated = samplesheet
        .join(SEQKIT_TRANSLATE.out.fasta_translated, by: 0)
        .map { meta, _fasta,  hmmfile, fasta_translated -> [meta, hmmfile, fasta_translated] }

    HMMSEARCH(
        ch_hmmfile_fastaTranslated,
    )

    ch_versions = ch_versions.mix(
        SEQKIT_TRANSLATE.out.versions,
        HMMSEARCH.out.versions,
    )


    emit:
    versions = ch_versions

}
