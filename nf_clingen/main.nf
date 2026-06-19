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

def inferBuildTag(value) {
    def text = value?.toString()?.toLowerCase()
    if( !text ) return null

    if( text =~ /grch38|hg38|b38/ ) {
        return 'grch38'
    }
    if( text =~ /grch37|hg37|b37/ ) {
        return 'grch37'
    }
    return null
}

def inferPanelBuildTag(panelPath) {
    if( !panelPath ) return null

    def path = panelPath.toString()
    def candidates = []

    if( new File(path).exists() && new File(path).isDirectory() ) {
        candidates.addAll(new File(path).listFiles()?.collect { it.name } ?: [])
    }
    else {
        candidates << new File(path).name
    }

    def tags = candidates.collect { inferBuildTag(it) }.findAll { it }
    if( tags.contains('grch38') ) return 'grch38'
    if( tags.contains('grch37') ) return 'grch37'
    return null
}

def normalizeParams() {
    ['aligner', 'preprocess', 'caller', 'workflow', 'execution_mode'].each { key ->
        params[key] = params[key].toString().trim().toLowerCase()
    }

    ['skip_qc', 'help', 'validate_only'].each { key ->
        params[key] = asBoolean(params[key])
    }

    params.max_report_variants = params.max_report_variants as Integer
    params.deepvariant_memory = params.deepvariant_memory?.toString()
    params.deepvariant_time = params.deepvariant_time?.toString()
    if( params.clinvar_vcf?.toString()?.trim() in ['', 'null', 'none'] ) {
        params.clinvar_vcf = null
    }
    if( params.gnomad_vcf?.toString()?.trim() in ['', 'null', 'none'] ) {
        params.gnomad_vcf = null
    }
    if( params.samplesheet?.toString()?.trim() in ['', 'null', 'none'] ) {
        params.samplesheet = null
    }
    if( params.input_vcf?.toString()?.trim() in ['', 'null', 'none'] ) {
        params.input_vcf = null
    }
    if( params.input_vcf_tbi?.toString()?.trim() in ['', 'null', 'none'] ) {
        params.input_vcf_tbi = null
    }
    if( params.input_sample_id?.toString()?.trim() in ['', 'null', 'none'] ) {
        params.input_sample_id = null
    }
}

def validateParams() {
    def spec = schema()

    if( !params.input_vcf && !params.samplesheet && (!params.reads || params.reads.toString().trim() == '') ) {
        throw new IllegalArgumentException("Provide either --input_vcf, --reads or --samplesheet")
    }

    if( params.execution_mode !in ['local', 'aws'] ) {
        throw new IllegalArgumentException("Unsupported --execution_mode '${params.execution_mode}'. Supported values: local, aws")
    }

    if( params.execution_mode == 'aws' ) {
        if( !params.aws_workdir || !params.aws_workdir.toString().startsWith('s3://') ) {
            throw new IllegalArgumentException("AWS mode requires --aws_workdir with an s3:// path")
        }
        if( params.scratch_dir ) {
            log.warn "Ignoring --scratch_dir in AWS mode"
        }
    }

    if( params.workflow == 'genealogy' ) {
        def referenceBuild = inferBuildTag(params.reference)
        def eagleMapBuild = inferPanelBuildTag(params.eagle2_genetic_map)
        def panelBuild = inferPanelBuildTag(params.beagle_ref_panel)

        def knownBuilds = [referenceBuild, eagleMapBuild, panelBuild].findAll { it }
        if( knownBuilds && knownBuilds.unique().size() > 1 ) {
            throw new IllegalArgumentException(
                "Genealogy build mismatch detected: reference=${referenceBuild ?: 'unknown'}, " +
                "eagle2_map=${eagleMapBuild ?: 'unknown'}, beagle_panel=${panelBuild ?: 'unknown'}. " +
                "Use a consistent genome build across reference, genetic map, and reference panel files."
            )
        }

        if( referenceBuild && panelBuild && referenceBuild != panelBuild ) {
            throw new IllegalArgumentException(
                "Genealogy reference build (${referenceBuild}) does not match Beagle panel build (${panelBuild}). " +
                "Use a matching panel set or a matching reference FASTA."
            )
        }
    }

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
                    // Some config defaults (for example MemoryUnit/Duration) are typed
                    // objects at parse time; normalize them to plain strings for schema
                    // validation and consistent downstream behavior.
                    params[key] = value.toString()
                    value = params[key]
                }
                break
        }

        if( property.enum && !property.enum.contains(value) ) {
            throw new IllegalArgumentException("Unsupported --${key} '${value}'. Supported values: ${property.enum.join(', ')}")
        }
    }
}

def samplesChannelFromSheet(sheetPath) {
    Channel
        .fromPath(sheetPath, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            def sampleId = row.sample_id?.toString()?.trim()
            def r1 = row.r1?.toString()?.trim()
            def r2 = row.r2?.toString()?.trim()

            if( !sampleId || !r1 || !r2 ) {
                throw new IllegalArgumentException("Samplesheet requires columns sample_id,r1,r2 with non-empty values")
            }

            tuple(sampleId, [file(r1, checkIfExists: true), file(r2, checkIfExists: true)])
        }
}

def renderHelp() {
    def spec = schema()
    def lines = []
    lines << 'nf-clingen'
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
                     nf-clingen runtime
    ==============================================
    Reads pattern         : ${params.reads}
    Input VCF             : ${params.input_vcf ?: 'disabled'}
    Reference FASTA       : ${params.reference}
    Output directory      : ${params.outdir}
    Skip QC               : ${params.skip_qc}
    Alignment engine      : ${params.aligner}
    Preprocess engine     : ${params.preprocess}
    Variant caller        : ${params.caller}
    Endpoint workflow     : ${params.workflow}
    ClinVar resource      : ${params.clinvar_vcf ?: 'disabled'}
    gnomAD resource       : ${params.gnomad_vcf ?: 'disabled'}
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

    def reference_source_ch = Channel.of(file(params.reference, checkIfExists: true))
    def called_vcf_ch

    if( params.input_vcf ) {
        def input_vcf = file(params.input_vcf, checkIfExists: true)
        def input_vcf_index = file(params.input_vcf_tbi ?: "${params.input_vcf}.tbi", checkIfExists: true)
        def sample_id = params.input_sample_id ?: input_vcf.name.replaceFirst(/\.vcf(\.gz)?$|\.bcf$/, '')
        log.info "VCF-only routing: skipping QC, alignment, preprocessing, and variant calling. Reusing ${input_vcf}."
        called_vcf_ch = Channel.of(tuple(sample_id, input_vcf, input_vcf_index))
    }
    else {
        def samples_ch = params.samplesheet
            ? samplesChannelFromSheet(params.samplesheet)
            : Channel.fromFilePairs(params.reads, checkIfExists: true)

        def upstream = SF_UPSTREAM(samples_ch, reference_source_ch)
        def aligned = SF_ALIGNMENT(upstream.analysis_reads_ch, upstream.fasta_index_ch, upstream.bwa_reference_ch)
        def preprocessed = SF_PREPROCESS(aligned.aligned_ch, upstream.fasta_index_ch)
        def called = SF_CALLING(preprocessed.processed_bam_ch, reference_source_ch, upstream.fasta_index_ch, upstream.fasta_index_ch)
        called_vcf_ch = called.called_vcf_ch
    }

    SF_ENDPOINTS(called_vcf_ch)
}

workflow.onComplete {
    log.info "nf-clingen finished with status: ${workflow.success ? 'SUCCESS' : 'FAILED'}"
    log.info "Published results: ${params.outdir}"
}