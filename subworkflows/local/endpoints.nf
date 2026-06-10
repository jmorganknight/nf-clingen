include { INDEX_CLINVAR_RESOURCE; RUN_CLINICAL_ANNOTATION; RUN_CLINVAR_ANNOTATION; BUILD_CLINICAL_REPORT } from '../../modules/local/clinical'
include { FILTER_BY_PHENOTYPE } from '../../modules/local/phenotype'
include { PHASE_EAGLE2; IMPUTE_BEAGLE; BUILD_GENEALOGY_MANIFEST } from '../../modules/local/genealogy'

workflow SF_ENDPOINTS {
    take:
    called_vcf_ch

    main:
    switch( params.workflow ) {
        case 'clinical':
            log.info 'Endpoint routing: clinical annotation and PDF compilation selected.'
            def clinical_annotation_ch
            if( params.clinvar_vcf ) {
                log.info "Clinical annotations: ClinVar overlay enabled from ${params.clinvar_vcf}."
                def clinvar_resource_ch = INDEX_CLINVAR_RESOURCE(Channel.value(file(params.clinvar_vcf, checkIfExists: true))).resource
                def clinvar_input_ch = called_vcf_ch
                    .combine(clinvar_resource_ch)
                    .map { row -> tuple(row[0], row[1], row[2], row[3], row[4]) }
                clinical_annotation_ch = RUN_CLINVAR_ANNOTATION(clinvar_input_ch).annotation
            }
            else {
                log.info 'Clinical annotations: using intrinsic VCF annotations only.'
                clinical_annotation_ch = RUN_CLINICAL_ANNOTATION(called_vcf_ch).annotation
            }
            
            // Apply phenotype filtering if patient phenotype is provided
            def report_input_ch
            if( params.patient_phenotype ) {
                log.info "Phenotype filtering enabled: ${params.patient_phenotype}"
                def phenotype_map_file = file(params.phenotype_gene_map, checkIfExists: true)
                def phenotype_filter_input = clinical_annotation_ch
                    .map { row -> tuple(row[0], row[3], params.patient_phenotype, phenotype_map_file) }
                FILTER_BY_PHENOTYPE(phenotype_filter_input)
                // Build report input: join and use phenotype-filtered TSV (row[5]) in place of raw annotation
                // Join result: [sample_id, vcf, vcf_tbi, annotation_tsv, summary_json, phenotype_filtered_tsv]
                report_input_ch = clinical_annotation_ch
                    .join(FILTER_BY_PHENOTYPE.out.filtered)
                    .map { row -> tuple(row[0], row[1], row[2], row[5], row[4]) }
            }
            else {
                report_input_ch = clinical_annotation_ch
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
            def phased_ch   = PHASE_EAGLE2(called_vcf_ch, genetic_map_ch).phased
            def imputed_ch  = IMPUTE_BEAGLE(phased_ch, ref_panel_ch, beagle_map_ch).imputed
            def merged_ch   = phased_ch
                .join(imputed_ch)
                .map { row -> tuple(row[0], row[1], row[2], row[3], row[4]) }
            BUILD_GENEALOGY_MANIFEST(merged_ch)
            break
    }
}
