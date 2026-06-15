include { ALIGN_MINIMAP2; ALIGN_BWAMEM2 } from '../../modules/local/alignment'

workflow SF_ALIGNMENT {
    take:
    analysis_reads_ch
    fasta_for_alignment_ch
    bwa_reference_ch

    main:
    // Normalizes read tuples from either:
    // 1) fromFilePairs => (sample_id, [r1, r2])
    // 2) FASTP output  => (sample_id, r1, r2)
    def normalized_reads_ch = analysis_reads_ch.map { row ->
        def reads = null
        if( row.size() > 2 ) {
            reads = [row[1], row[2]]
        }
        else if( row[1] instanceof Collection ) {
            reads = row[1].toList()
        }
        else {
            throw new IllegalArgumentException("Unsupported read tuple shape for sample '${row[0]}'")
        }
        tuple(row[0], reads)
    }

    aligned_ch = null
    switch( params.aligner ) {
        case 'minimap2':
            log.info 'Alignment routing: minimap2 short-read mode selected.'
            def minimap_input_ch = normalized_reads_ch
                .combine(fasta_for_alignment_ch)
                .map { row -> tuple(row[0], row[1], row[2]) }
            aligned_ch = ALIGN_MINIMAP2(minimap_input_ch).aligned_sam
            break
        case 'bwamem2':
            log.info 'Alignment routing: bwamem2 selected with explicit reference indexing.'
            def bwa_input_ch = normalized_reads_ch
                .combine(bwa_reference_ch)
                .map { row -> tuple(row[0], row[1], row[2], row[3]) }
            aligned_ch = ALIGN_BWAMEM2(bwa_input_ch).aligned_sam
            break
    }

    emit:
    aligned_ch
}
