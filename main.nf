import groovy.yaml.YamlBuilder

include { samplesheetToList } from 'plugin/nf-schema'

include { TGSCAN } from './workflows/tgscan'

help_message = """
tgscan.nf
======================================================================
Scan a six-frame translated genome for sequences of interest using HMM

Usage:
    nextflow run exscan.nf --input <samplesheet.csv> [options]
    nextflow run exscan.nf -params-file <params.yaml>

Required Arguments:
  --input               : Path to comma-separated file containing each sample to be processed.
  --outdir              : Path to the output directory where results will be saved

Alternative to Required Arguments:
  -params-file <yaml>   : YAML/JSON file with the parameters

Pipeline Parameters:

  --translate_min_len   : Minimum length of a translated sequence to be retained for
                          downstream analysis. [default: ${params.translate_min_len}]

  --seq_eval            : Cosider targets with a per sequence E-value <= of this value as true positives.
                          Mutually exclusive with --seq_bitscore. [default: ${params.seq_eval}]

  --seq_bitscore        : Consider targets with a per sequence bitscore >= of this value as true positives.
                          Mutually exclusive with --seq_eval. [default: ${params.seq_bitscore}]

  --dom_eval            : For target sequences that sequences that have already satisfied the
                          --seq_eval or --seq_bitscore thresholds, consider domains with a
                          E-value <= of this value as true positives.
                          Mutually exclusive with --dom_bitscore. [default: ${params.dom_eval}]

  --dom_bitscore        : For target sequences that sequences that have already satisfied the
                          --seq_eval or --seq_bitscore thresholds, consider domains with a
                          E-value <= of this value as true positives.
                          per domain bitscore >= of this value.
                          Given value will be passed as '--domT' argument to hmmsearch.
                          Mutually exclusive with --dom_eval. [default: ${params.dom_bitscore}]

  --z_seq_evalue        : Set the #  of comparasions made for sequence E-value calculation.
                          Given value will be passed as '-Z' argument to hmmsearch.
                          [default: ${params.z_seq_evalue}]

  --z_dom_evalue        : Set the #  of comparasions made for domain E-value calculation.
                          Given value will be passed as '--domZ' argument to hmmsearch.
                          [default: ${params.z_dom_evalue}]

  --togff_mode          : Mode of operation when writing profile alignments coordinates as GFF3.
                          Options are:
                            all-alignments   : All alignments from all profiles
                            best-profile     : All alignments from the best performing profile
                            best-alignment   : Single best alignment from the best performing profile
                            best-per-profile : Single best alignment from each profile
                          [default: ${params.togff_mode}]

  --togff_metric        : Metric to use for selecting the best scoring profile or alignment.
                          Options are:
                            bitscore         : Bitscore
                            evalue           : E-value for profile, conditional E-value for alignment
                          [default: ${params.togff_metric}]

  --help                : Print help message and exit

  --version             : Print version and exit
"""


init_summary = """
T G S C A N   P I P E L I N E   v${params.manifest.version}
======================================
input                   : ${params.input}
outdir                  : ${params.outdir}
translate_min_len       : ${params.translate_min_len}
seq_eval                : ${params.seq_eval}
seq_bitscore            : ${params.seq_bitscore}
dom_eval                : ${params.dom_eval}
dom_bitscore            : ${params.dom_bitscore}
z_seq_evalue            : ${params.z_seq_evalue}
z_dom_evalue            : ${params.z_dom_evalue}
togff_mode              : ${params.togff_mode}
togff_metric            : ${params.togff_metric}

--

Run as                  : ${workflow.commandLine}
Started at              : ${workflow.start}
Config files            : ${workflow.configFiles}

--
"""
// container images : ${workflow.containerEngine}:${workflow.container}


// DESC: Validate input arguments and initialize pipeline, printing a small summary
// ARGS: None, uses variables defined at the beginning of the script
// OUTS: None
// RETS: None
def validateParams() {

    // `--help` and `--version` flags
    if (params.help) {
        println help_message
        System.exit(0)
    }
    if (params.version) {
        println "${params.manifest.name} v${params.manifest.version}"
        System.exit(0)
    }

    // Check required arguments
    if (params.input == null) {
        println help_message
        log.error "Missing required argument: --input"
        System.exit(1)
    }
    if (!file(params.input).exists()) {
        log.error "File not found: ${params.input}"
        System.exit(1)
    }
    if (params.outdir == null) {
        println help_message
        log.error "Missing required argument: --outdir"
        System.exit(1)
    }

    // Check for mutually exclusive arguments
    if (params.seq_eval && params.seq_bitscore) {
        log.error "Arguments --seq_eval and --seq_bitscore are mutually exclusive"
        System.exit(1)
    }
    if (params.dom_eval && params.dom_bitscore) {
        log.error "Arguments --dom_eval and --dom_bitscore are mutually exclusive"
        System.exit(1)
    }
    if (params.seq_eval == null && params.seq_bitscore == null) {
        log.error "At least one of --seq_eval or --seq_bitscore must be provided"
        System.exit(1)
    }
    if (params.dom_eval == null && params.dom_bitscore == null) {
        log.error "At least one of --dom_eval or --dom_bitscore must be provided"
        System.exit(1)
    }

    // Check for valid togff_mode and togff_metric values
    def valid_togff_modes = ['all-alignments', 'best-profile', 'best-alignment', 'best-per-profile']
    if (!valid_togff_modes.contains(params.togff_mode)) {
        log.error "Invalid value for --togff_mode: ${params.togff_mode}. Valid options are: {${valid_togff_modes.join(', ')}}"
        System.exit(1)
    }

    def valid_togff_metrics = ['bitscore', 'evalue']
    if (!valid_togff_metrics.contains(params.togff_metric)) {
        log.error "Invalid value for --togff_metric: ${params.togff_metric}. Valid options are: {${valid_togff_metrics.join(', ')}}"
        System.exit(1)
    }

    // Check for valid z_seq_evalue and z_dom_evalue values
    if (params.z_seq_evalue < 0) {
        log.error "Invalid value for --z_seq_evalue: ${params.z_seq_evalue}. Must be >= 0"
        System.exit(1)
    }

    if (params.z_dom_evalue < 0) {
        log.error "Invalid value for --z_dom_evalue: ${params.z_dom_evalue}. Must be >= 0"
        System.exit(1)
    }

    // Check for valid translate_min_len value
    if (params.translate_min_len < 0) {
        log.error "Invalid value for --translate_min_len: ${params.translate_min_len}. Must be >= 0"
        System.exit(1)
    }
}


// DESC: Dump parameters to a YAML file
// ARGS: None, uses variables defined at the beginning of the script
// OUTS: None
// RETS: Channel with the samplesheet converted to a list
def dumpParametersToYaml() {
    def paramMap = [
        input: params.input,
        outdir: params.outdir,
        translate_min_len: params.translate_min_len,
        seq_eval: params.seq_eval,
        seq_bitscore: params.seq_bitscore,
        dom_eval: params.dom_eval,
        dom_bitscore: params.dom_bitscore,
    ]

    def yaml = new YamlBuilder()
    yaml(paramMap)

    def outputDir = new File("${params.outdir}/pipeline_info")
    outputDir.mkdirs()
    def outputFile = new File(outputDir, "params_used.yaml")
    outputFile.text = yaml.toString()

}


// DESC: Parse the samplesheet and convert it to a list
// ARGS: None, uses variables defined at the beginning of the script
// OUTS: Channel with the samplesheet converted to a list
// RETS: Channel with the samplesheet converted to a list
def parseSamplesheet() {
    def ch_samplesheet = Channel
        .fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))

    return ch_samplesheet
}

// DESC: Display completion message based on workflow status
// ARGS: None, uses variables defined at the beginning of the script
// OUTS: Completion message at `INFO` or `ERROR` level
// RETS: None
def completionMsg() {

    if (workflow.success) {
        if (workflow.stats.ignoredCount == 0) {
            log.info "Pipeline completed successfully!"
        }
        else {
            log.info "Pipeline completed successully, but with errored processes"
        }
    }
    else {
        log.error "Pipeline completed with errors"
    }

}

// Main workflow
workflow {

    main:

    // Validate input parameters
    validateParams()
    // Initialization Summary - Everything looks good so far
    log.info init_summary

    // Parse and validate samplesheet.csv from `--input`
    ch_samplesheet = parseSamplesheet()

    // Dump parameters to YAML file
    dumpParametersToYaml()


    TGSCAN(
        ch_samplesheet,
    )


    ch_versions = TGSCAN.out.versions
    ch_versions.collectFile(
        storeDir: "${params.outdir}/pipeline_info/",
        name: 'versions.yml',
        sort: true,
        newLine: true
    )

    // Display any error encountered during the workflow
    workflow.onComplete {
        completionMsg()
    }

}


