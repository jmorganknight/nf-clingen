include { INDEX_CLINVAR_RESOURCE; RUN_CLINICAL_ANNOTATION; RUN_CLINVAR_ANNOTATION; BUILD_CLINICAL_REPORT } from '../../modules/local/clinical'
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
            BUILD_CLINICAL_REPORT(clinical_annotation_ch)
            break
        case 'genealogy':
            log.info 'Endpoint routing: genealogy phasing/imputation stub selected.'
            RUN_GENEALOGY_STACK(called_vcf_ch)
            break
    }
}
