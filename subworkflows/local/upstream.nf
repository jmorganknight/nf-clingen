include { RUN_FASTQC; RUN_FASTP } from '../../modules/local/qc'
include { BUILD_FASTA_INDEX; BUILD_BWAMEM2_INDEX } from '../../modules/local/reference'

workflow SF_UPSTREAM {
    take:
    samples_ch
    reference_source_ch

    main:
    fasta_index_result = BUILD_FASTA_INDEX(reference_source_ch)
    fasta_index_ch = fasta_index_result.reference

    if( params.skip_qc ) {
        log.info 'QC routing: FASTQC and FASTP skipped; raw reads move directly to alignment.'
        analysis_reads_ch = samples_ch
    }
    else {
        log.info 'QC routing: FASTQC and FASTP enabled.'
        RUN_FASTQC(samples_ch)
        analysis_reads_ch = RUN_FASTP(samples_ch).trimmed_reads
    }

    bwa_reference_ch = Channel.empty()
    if( params.aligner == 'bwamem2' ) {
        bwa_reference_ch = BUILD_BWAMEM2_INDEX(reference_source_ch).reference
    }

    emit:
    analysis_reads_ch
    fasta_index_ch
    bwa_reference_ch
}
