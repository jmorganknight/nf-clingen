#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

import groovy.json.JsonSlurper

include { SF_UPSTREAM } from './subworkflows/local/upstream'
include { SF_ALIGNMENT } from './subworkflows/local/alignment'
include { SF_PREPROCESS } from './subworkflows/local/preprocess'
include { SF_CALLING } from './subworkflows/local/calling'
include { SF_ENDPOINTS } from './subworkflows/local/endpoints'

def schema() {
    new JsonSlurper().parse(file("${projectDir}/nextflow_schema.json"))
}

def asBoolean(value) {
    if( value instanceof Boolean ) {
        return value
    }
    return value?.toString()?.trim()?.toLowerCase() in ['1', 'true', 'yes', 'y']
}

def normalizeParams() {
    ['aligner', 'preprocess', 'caller', 'workflow'].each { key ->
        params[key] = params[key].toString().trim().toLowerCase()
    }

    ['skip_qc', 'help', 'validate_only'].each { key ->
        params[key] = asBoolean(params[key])
    }

    params.max_report_variants = params.max_report_variants as Integer
    if( params.clinvar_vcf?.toString()?.trim() in ['', 'null', 'none'] ) {
        params.clinvar_vcf = null
    }
}

def validateParams() {
    def spec = schema()

    spec.required.each { key ->
        def value = params[key]
        if( value == null || value.toString().trim() == '' ) {
            throw new IllegalArgumentException("Missing required parameter --${key}")
        }
    }

    spec.properties.each { key, property ->
        def value = params[key]
        if( value == null ) {
            return
        }

        switch( property.type ) {
            case 'boolean':
                if( !(value instanceof Boolean) ) {
                    throw new IllegalArgumentException("Parameter --${key} must be a boolean")
                }
                break
            case 'integer':
                if( !(value instanceof Integer) ) {
                    throw new IllegalArgumentException("Parameter --${key} must be an integer")
                }
                if( property.minimum != null && value < (property.minimum as Integer) ) {
                    throw new IllegalArgumentException("Parameter --${key} must be >= ${property.minimum}")
                }
                break
            case 'string':
                if( !(value instanceof CharSequence) ) {
                    throw new IllegalArgumentException("Parameter --${key} must be a string")
                }
                break
        }

        if( property.enum && !property.enum.contains(value) ) {
            throw new IllegalArgumentException("Unsupported --${key} '${value}'. Supported values: ${property.enum.join(', ')}")
        }
    }
}

def renderHelp() {
    def spec = schema()
    def lines = []
    lines << 'nf-prism'
    lines << ''
    lines << spec.description.toString()
    lines << ''
    lines << 'Parameters:'

    spec.properties.each { key, property ->
        def defaultValue = property.containsKey('default') ? property.default : '(none)'
        def enumText = property.enum ? " Allowed: ${property.enum.join(', ')}." : ''
        lines << String.format('  --%-20s %s Default: %s.%s', key, property.description.toString(), defaultValue, enumText)
    }

    lines << ''
    lines << 'Example:'
    lines << "  nextflow run . -profile docker --workflow clinical --caller deepvariant --clinvar_vcf ${projectDir}/assets/clinvar.vcf.gz"
    lines.join('\n')
}

def runtimeBanner() {
    """
    ==============================================
                     nf-prism runtime
    ==============================================
    Reads pattern         : ${params.reads}
    Reference FASTA       : ${params.reference}
    Output directory      : ${params.outdir}
    Skip QC               : ${params.skip_qc}
    Alignment engine      : ${params.aligner}
    Preprocess engine     : ${params.preprocess}
    Variant caller        : ${params.caller}
    Endpoint workflow     : ${params.workflow}
    ClinVar resource      : ${params.clinvar_vcf ?: 'disabled'}
    Report variant limit  : ${params.max_report_variants}
    ==============================================
    Reference indexes are built on demand inside the workflow.
    Clinical reporting uses bcftools-based local triage with optional ClinVar overlay and WeasyPrint PDF generation.
    """.stripIndent()
}

workflow {
    normalizeParams()

    if( params.help ) {
        log.info(renderHelp())
        return
    }

    validateParams()
    log.info(runtimeBanner())

    if( params.validate_only ) {
        log.info 'Schema validation succeeded; exiting because --validate_only was set.'
        return
    }

    def samples_ch = Channel
        .fromFilePairs(params.reads, checkIfExists: true)

    def reference_source_ch = Channel.of(file(params.reference, checkIfExists: true))

    def upstream = SF_UPSTREAM(samples_ch, reference_source_ch)
    def aligned = SF_ALIGNMENT(upstream.analysis_reads_ch, upstream.fasta_index_ch, upstream.bwa_reference_ch)
    def preprocessed = SF_PREPROCESS(aligned.aligned_ch, upstream.fasta_index_ch)
    def called = SF_CALLING(preprocessed.processed_bam_ch, reference_source_ch, upstream.fasta_index_ch, upstream.fasta_index_ch)
    SF_ENDPOINTS(called.called_vcf_ch)
}

workflow.onComplete {
    log.info "nf-prism finished with status: ${workflow.success ? 'SUCCESS' : 'FAILED'}"
    log.info "Published results: ${params.outdir}"
}