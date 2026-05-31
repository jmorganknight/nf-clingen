process ALIGN_MINIMAP2 {
    tag "${sample_id}"
    publishDir "${params.outdir}/alignment", mode: 'copy', pattern: "${sample_id}.aligned.sam"

    input:
    tuple val(sample_id), path(reads), path(reference)

    output:
    tuple val(sample_id), path("${sample_id}.aligned.sam"), emit: aligned_sam

    script:
    def read1 = reads[0]
    def read2 = reads[1]
    """
    set -euo pipefail
    minimap2 -t ${task.cpus} -ax sr ${reference} ${read1} ${read2} > ${sample_id}.aligned.sam
    """

    stub:
    """
    touch ${sample_id}.aligned.sam
    """
}

process ALIGN_BWAMEM2 {
    tag "${sample_id}"
    publishDir "${params.outdir}/alignment", mode: 'copy', pattern: "${sample_id}.aligned.sam"

    input:
    tuple val(sample_id), path(reads), path(reference), path(index_files)

    output:
    tuple val(sample_id), path("${sample_id}.aligned.sam"), emit: aligned_sam

    script:
    def read1 = reads[0]
    def read2 = reads[1]
    """
    set -euo pipefail
    bwa-mem2 mem -t ${task.cpus} ${reference} ${read1} ${read2} > ${sample_id}.aligned.sam
    """

    stub:
    """
    touch ${sample_id}.aligned.sam
    """
}