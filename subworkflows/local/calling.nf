include { BUILD_SEQUENCE_DICTIONARY } from '../../modules/local/reference'
include { CALL_GATK; FILTER_GATK; CALL_DEEPVARIANT; VQSR_GATK } from '../../modules/local/variant_calling'

workflow SF_CALLING {
    take:
    processed_bam_ch
    reference_source_ch
    fasta_for_deepvariant_ch
    fasta_for_gatk_ch

    main:
    called_vcf_ch = null
    if( params.caller == 'haplotypecaller' ) {
        log.info 'Variant routing: GATK HaplotypeCaller selected.'
        def dict_ch = BUILD_SEQUENCE_DICTIONARY(reference_source_ch).dict
        def gatk_ref_bundle_ch = fasta_for_gatk_ch
            .combine(dict_ch)
            .map { row -> tuple(row[0], row[1], row[2]) }
        def gatk_input_ch = processed_bam_ch
            .combine(gatk_ref_bundle_ch)
            .map { row -> tuple(row[0], row[1], row[2], row[3], row[4], row[5]) }
        def raw_vcf_ch = CALL_GATK(gatk_input_ch).raw_vcf
        def dbsnp_resource = params.vqsr_dbsnp_resource ?: params.vqsr_snp_resource
        def has_all_vqsr_resources = [
            params.vqsr_snp_resource,
            params.vqsr_omni_resource,
            params.vqsr_indel_resource,
            params.vqsr_1kg_resource,
            dbsnp_resource
        ].every { it != null && it.toString().trim() }

        if( has_all_vqsr_resources ) {
            log.info 'VQSR routing: all required resources detected; applying VQSR recalibration.'
            def vqsr_resources_ch = channel
                .value([params.vqsr_snp_resource, params.vqsr_omni_resource, params.vqsr_indel_resource, params.vqsr_1kg_resource, dbsnp_resource])
            def vqsr_input_ch = raw_vcf_ch
                .combine(gatk_ref_bundle_ch)
                .combine(vqsr_resources_ch)
                .map { row -> tuple(row[0], row[1], row[2], row[3], row[4], row[5], row[6], row[7], row[8], row[9], row[10]) }
            called_vcf_ch = VQSR_GATK(vqsr_input_ch).vcf
        }
        else {
            log.warn 'VQSR resources incomplete or missing; falling back to GATK hard filtering.'
            def hard_filter_input_ch = raw_vcf_ch
                .combine(gatk_ref_bundle_ch)
                .map { row -> tuple(row[0], row[1], row[2], row[3], row[4], row[5]) }
            called_vcf_ch = FILTER_GATK(hard_filter_input_ch).vcf
        }
    }
    else if( params.caller == 'deepvariant' ) {
        log.info 'Variant routing: DeepVariant selected.'
        def deepvariant_input_ch = processed_bam_ch
            .combine(fasta_for_deepvariant_ch)
            .map { row -> tuple(row[0], row[1], row[2], row[3], row[4]) }
        called_vcf_ch = CALL_DEEPVARIANT(deepvariant_input_ch).vcf
    }

    emit:
    called_vcf_ch
}
