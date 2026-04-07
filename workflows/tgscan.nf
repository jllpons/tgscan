include { TRANSLATE                     } from '../modules/local/easel/main.nf'
include { HMMSEARCH                     } from '../modules/local/hmmer/main.nf'
include { HMMSEARCH_DOMTBLOUT_TO_GFF    } from '../modules/local/awk/main.nf'


workflow TGSCAN {


    take:
    samplesheet // channel: samplesheet read in from --input

    main:

    ch_versions = Channel.empty()

    ch_genome = samplesheet
        .map { meta, genome, _hmm ->
            [meta, genome]
        }

    ch_hmm  = samplesheet
        .map { meta, _genome, hmm ->
            [meta, hmm]
        }

    // Step 1: 6-frame translation of the input genome to individual ORFs
    TRANSLATE(
        ch_genome,
    )
    ch_versions = ch_versions.mix(TRANSLATE.out.versions)

    ch_hmmserarch_input = TRANSLATE.out.orfs
        .combine(ch_hmm)

    // Step 2: HMMER search to query the provided HMMs against the translated ORFs
    HMMSEARCH(
        ch_hmm,
        TRANSLATE.out.orfs
    )
    ch_versions = ch_versions.mix(HMMSEARCH.out.versions)

    // Step 3: Convert HMMER domtblout output to GFF format for downstream analysis
    HMMSEARCH_DOMTBLOUT_TO_GFF(
        HMMSEARCH.out.domtbloutput,
    )

    emit:
    versions = ch_versions
}
