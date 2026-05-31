process RUN_FASTQC {
    tag "${sample_id}"
    publishDir "${params.outdir}/qc/fastqc", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    path("*.html"), emit: reports
    path("*.zip"), emit: archives

    script:
    def read1 = reads[0]
    def read2 = reads[1]
    """
    set -euo pipefail
    fastqc --threads ${task.cpus} --outdir . ${read1} ${read2}
    """

    stub:
    """
    touch ${sample_id}_R1_fastqc.html ${sample_id}_R2_fastqc.html
    touch ${sample_id}_R1_fastqc.zip ${sample_id}_R2_fastqc.zip
    """
}

process RUN_FASTP {
    tag "${sample_id}"
    publishDir "${params.outdir}/qc/fastp", mode: 'copy', pattern: "${sample_id}.*"

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_R1.trimmed.fastq.gz"), path("${sample_id}_R2.trimmed.fastq.gz"), emit: trimmed_reads
    path("${sample_id}.fastp.json"), emit: json
    path("${sample_id}.fastp.html"), emit: html

    script:
    def read1 = reads[0]
    def read2 = reads[1]
    """
    set -euo pipefail
    fastp \
      --in1 ${read1} \
      --in2 ${read2} \
      --out1 ${sample_id}_R1.trimmed.fastq.gz \
      --out2 ${sample_id}_R2.trimmed.fastq.gz \
      --json ${sample_id}.fastp.json \
      --html ${sample_id}.fastp.html \
      --thread ${task.cpus} \
      --detect_adapter_for_pe
    """

    stub:
    """
    touch ${sample_id}_R1.trimmed.fastq.gz ${sample_id}_R2.trimmed.fastq.gz
    touch ${sample_id}.fastp.json ${sample_id}.fastp.html
    """
}