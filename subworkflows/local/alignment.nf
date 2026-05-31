include { ALIGN_MINIMAP2; ALIGN_BWAMEM2 } from '../../modules/local/alignment'

workflow SF_ALIGNMENT {
    take:
    analysis_reads_ch
    fasta_for_alignment_ch
    bwa_reference_ch

    main:
    aligned_ch = null
    switch( params.aligner ) {
        case 'minimap2':
            log.info 'Alignment routing: minimap2 short-read mode selected.'
            def minimap_input_ch = analysis_reads_ch
                .combine(fasta_for_alignment_ch)
                .map { row -> tuple(row[0], row[1], row[2]) }
            aligned_ch = ALIGN_MINIMAP2(minimap_input_ch).aligned_sam
            break
        case 'bwamem2':
            log.info 'Alignment routing: bwamem2 selected with explicit reference indexing.'
            def bwa_input_ch = analysis_reads_ch
                .combine(bwa_reference_ch)
                .map { row -> tuple(row[0], row[1], row[2], row[3]) }
            aligned_ch = ALIGN_BWAMEM2(bwa_input_ch).aligned_sam
            break
    }

    emit:
    aligned_ch
}
