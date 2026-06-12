include { INDEX_CLINVAR_RESOURCE; INDEX_GNOMAD_RESOURCE; RUN_GNOMAD_ANNOTATION; RUN_CLINICAL_ANNOTATION; RUN_CLINVAR_ANNOTATION; BUILD_CLINICAL_REPORT } from '../../modules/local/clinical'
include { FILTER_BY_PHENOTYPE } from '../../modules/local/phenotype'
include { PHASE_EAGLE2; IMPUTE_BEAGLE; CONCAT_IMPUTED; BUILD_GENEALOGY_MANIFEST } from '../../modules/local/genealogy'

workflow SF_ENDPOINTS {
    take:
    called_vcf_ch

    main:
    switch( params.workflow ) {
        case 'clinical':
            log.info 'Endpoint routing: clinical annotation and PDF compilation selected.'
            
            // Validate and resolve report_mode
            def report_mode = params.report_mode
            if( report_mode == 'auto' ) {
                report_mode = params.patient_phenotype ? 'phenotype' : 'clinical'
            }
            log.info "Report mode: ${report_mode} (patient_phenotype: ${params.patient_phenotype ?: 'none'})"
            
            if( report_mode in ['phenotype', 'combined'] && !params.patient_phenotype ) {
                error "report_mode '${report_mode}' requires --patient_phenotype to be set"
            }
            
            // Optional gnomAD pre-annotation for AF-based clinical triage.
            def annotation_input_ch = called_vcf_ch
            if( params.gnomad_vcf ) {
                log.info "Clinical annotations: gnomAD AF overlay enabled from ${params.gnomad_vcf}."
                def gnomad_resource_ch = INDEX_GNOMAD_RESOURCE(Channel.value(file(params.gnomad_vcf, checkIfExists: true))).resource
                def gnomad_input_ch = called_vcf_ch
                    .combine(gnomad_resource_ch)
                    .map { row -> tuple(row[0], row[1], row[2], row[3], row[4]) }
                annotation_input_ch = RUN_GNOMAD_ANNOTATION(gnomad_input_ch).annotated
            }

            // Always run clinical annotation (optionally after gnomAD AF annotation)
            def clinical_annotation_ch
            if( params.clinvar_vcf ) {
                log.info "Clinical annotations: ClinVar overlay enabled from ${params.clinvar_vcf}."
                def clinvar_resource_ch = INDEX_CLINVAR_RESOURCE(Channel.value(file(params.clinvar_vcf, checkIfExists: true))).resource
                def clinvar_input_ch = annotation_input_ch
                    .combine(clinvar_resource_ch)
                    .map { row -> tuple(row[0], row[1], row[2], row[3], row[4]) }
                clinical_annotation_ch = RUN_CLINVAR_ANNOTATION(clinvar_input_ch).annotation
            }
            else {
                log.info 'Clinical annotations: using intrinsic VCF annotations only.'
                clinical_annotation_ch = RUN_CLINICAL_ANNOTATION(annotation_input_ch).annotation
            }
            
            // Apply phenotype filtering if report_mode requests it
            def report_input_ch
            if( report_mode == 'clinical' ) {
                // Use full clinical annotation for report
                report_input_ch = clinical_annotation_ch
            }
            else if( report_mode in ['phenotype', 'combined'] ) {
                // Run phenotype filtering and use filtered variants for report
                log.info "Phenotype filtering enabled: ${params.patient_phenotype}"
                def phenotype_map_file = file(params.phenotype_gene_map, checkIfExists: true)
                def phenotype_filter_input = clinical_annotation_ch
                    .map { row -> tuple(row[0], row[3], params.patient_phenotype, phenotype_map_file) }
                FILTER_BY_PHENOTYPE(phenotype_filter_input)
                
                // Use phenotype-filtered TSV for report
                report_input_ch = clinical_annotation_ch
                    .join(FILTER_BY_PHENOTYPE.out.filtered)
                    .map { row -> tuple(row[0], row[1], row[2], row[5], row[4]) }
            }

            def report_builder_script = file("${projectDir}/scripts/build_clinical_report.py", checkIfExists: true)
            report_input_ch = report_input_ch.map { row -> tuple(row[0], row[1], row[2], row[3], row[4], report_builder_script) }
            
            BUILD_CLINICAL_REPORT(report_input_ch)
            break
        case 'genealogy':
            log.info 'Endpoint routing: Eagle2 phasing → Beagle imputation selected.'
            if( !params.eagle2_genetic_map ) { error "--eagle2_genetic_map is required for workflow=genealogy" }
            if( !params.beagle_ref_panel )   { error "--beagle_ref_panel is required for workflow=genealogy" }
            def genetic_map_ch  = Channel.value(file(params.eagle2_genetic_map, checkIfExists: true))
            def ref_panel_ch    = Channel.value(file(params.beagle_ref_panel, checkIfExists: true))
            def beagle_map_ch   = params.beagle_genetic_map
                ? Channel.value(file(params.beagle_genetic_map, checkIfExists: true))
                : Channel.value(file("${projectDir}/data"))
            def phased_ch      = PHASE_EAGLE2(called_vcf_ch, genetic_map_ch, ref_panel_ch).phased
            def imputed_raw_ch = IMPUTE_BEAGLE(phased_ch, ref_panel_ch, beagle_map_ch).chrom_vcfs
            def imputed_ch     = CONCAT_IMPUTED(imputed_raw_ch).imputed
            def merged_ch      = phased_ch
                .join(imputed_ch)
                .map { row -> tuple(row[0], row[1], row[2], row[3], row[4]) }
            BUILD_GENEALOGY_MANIFEST(merged_ch)
            break
    }
}
