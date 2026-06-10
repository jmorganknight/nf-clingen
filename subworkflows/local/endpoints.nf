include { INDEX_CLINVAR_RESOURCE; RUN_CLINICAL_ANNOTATION; RUN_CLINVAR_ANNOTATION; BUILD_CLINICAL_REPORT } from '../../modules/local/clinical'
include { FILTER_BY_PHENOTYPE } from '../../modules/local/phenotype'
include { RUN_GENEALOGY_STACK } from '../../modules/local/genealogy'

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
                // Build report input from phenotype-filtered results
                report_input_ch = clinical_annotation_ch
                    .join(FILTER_BY_PHENOTYPE.out.filtered)
                    .map { row -> tuple(row[0], row[1], row[2], row[4], row[3]) }
            }
            else {
                report_input_ch = clinical_annotation_ch
            }

            def report_builder_script = file("${projectDir}/scripts/build_clinical_report.py", checkIfExists: true)
            report_input_ch = report_input_ch.map { row -> tuple(row[0], row[1], row[2], row[3], row[4], report_builder_script) }
            
            BUILD_CLINICAL_REPORT(report_input_ch)
            break
        case 'genealogy':
            log.info 'Endpoint routing: genealogy phasing/imputation stub selected.'
            RUN_GENEALOGY_STACK(called_vcf_ch)
            break
    }
}
