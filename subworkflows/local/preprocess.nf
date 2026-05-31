include { RUN_SAMTOOLS; RUN_ELPREP; INDEX_BAM } from '../../modules/local/preprocess'

workflow SF_PREPROCESS {
    take:
    aligned_ch
    fasta_for_elprep_ch

    main:
    processed_bam_ch = null
    switch( params.preprocess ) {
        case 'samtools':
            log.info 'Preprocess routing: samtools sort/index selected.'
            processed_bam_ch = RUN_SAMTOOLS(aligned_ch).bam
            break
        case 'elprep':
            log.info 'Preprocess routing: elPrep selected for reduced intermediate I/O.'
            def elprep_input_ch = aligned_ch
                .combine(fasta_for_elprep_ch)
                .map { row -> tuple(row[0], row[1], row[2]) }
            def elprep_bam_ch = RUN_ELPREP(elprep_input_ch).bam
            processed_bam_ch = INDEX_BAM(elprep_bam_ch).bam
            break
    }

    emit:
    processed_bam_ch
}
