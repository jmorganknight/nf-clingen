include { BUILD_SEQUENCE_DICTIONARY } from '../../modules/local/reference'
include { CALL_GATK; CALL_DEEPVARIANT } from '../../modules/local/variant_calling'

workflow SF_CALLING {
    take:
    processed_bam_ch
    reference_source_ch
    fasta_for_deepvariant_ch
    fasta_for_gatk_ch

    main:
    called_vcf_ch = null
    switch( params.caller ) {
        case 'haplotypecaller':
            log.info 'Variant routing: GATK HaplotypeCaller selected.'
            def dict_ch = BUILD_SEQUENCE_DICTIONARY(reference_source_ch).dict
            def gatk_ref_bundle_ch = fasta_for_gatk_ch
                .combine(dict_ch)
                .map { row -> tuple(row[0], row[1], row[2]) }
            def gatk_input_ch = processed_bam_ch
                .combine(gatk_ref_bundle_ch)
                .map { row -> tuple(row[0], row[1], row[2], row[3], row[4], row[5]) }
            called_vcf_ch = CALL_GATK(gatk_input_ch).vcf
            break
        case 'deepvariant':
            log.info 'Variant routing: DeepVariant selected.'
            def deepvariant_input_ch = processed_bam_ch
                .combine(fasta_for_deepvariant_ch)
                .map { row -> tuple(row[0], row[1], row[2], row[3], row[4]) }
            called_vcf_ch = CALL_DEEPVARIANT(deepvariant_input_ch).vcf
            break
    }

    emit:
    called_vcf_ch
}
