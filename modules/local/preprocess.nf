process RUN_SAMTOOLS {
    tag "${sample_id}"
    publishDir "${params.outdir}/preprocess", mode: 'copy', pattern: "${sample_id}.sorted.bam*"

    input:
    tuple val(sample_id), path(alignment_sam)

    output:
    tuple val(sample_id), path("${sample_id}.sorted.bam"), path("${sample_id}.sorted.bam.bai"), emit: bam

    script:
    """
    set -euo pipefail
    samtools sort -@ ${task.cpus} -O BAM -o ${sample_id}.sorted.bam ${alignment_sam}
    samtools index -@ ${task.cpus} ${sample_id}.sorted.bam
    """

    stub:
    """
    touch ${sample_id}.sorted.bam ${sample_id}.sorted.bam.bai
    """
}

process RUN_ELPREP {
    tag "${sample_id}"
    publishDir "${params.outdir}/preprocess", mode: 'copy', pattern: "${sample_id}.elprep.bam"

    input:
    tuple val(sample_id), path(alignment_sam), path(reference)

    output:
    tuple val(sample_id), path("${sample_id}.elprep.bam"), emit: bam

    script:
    """
    set -euo pipefail
    elprep sfm ${alignment_sam} ${sample_id}.elprep.bam \
      --sorting-order coordinate \
      --mark-duplicates \
      --nr-of-threads ${task.cpus} \
      --reference ${reference}
    """

    stub:
    """
    touch ${sample_id}.elprep.bam
    """
}

process INDEX_BAM {
    tag "${sample_id}"
    publishDir "${params.outdir}/preprocess", mode: 'copy', pattern: "${sample_id}.*"

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path(bam), path("${bam.name}.bai"), emit: bam

    script:
    """
    set -euo pipefail
    samtools index -@ ${task.cpus} ${bam}
    """

    stub:
    """
    touch ${bam}.bai
    """
}