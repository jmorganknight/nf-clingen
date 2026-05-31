process CALL_GATK {
    tag "${sample_id}"
    publishDir "${params.outdir}/variants", mode: 'copy', pattern: "${sample_id}.gatk.vcf.gz*"

    input:
    tuple val(sample_id), path(bam), path(bai), path(reference), path(fai), path(dict)

    output:
    tuple val(sample_id), path("${sample_id}.gatk.vcf.gz"), path("${sample_id}.gatk.vcf.gz.tbi"), emit: vcf

    script:
    """
    set -euo pipefail
    gatk \
      --java-options "-Xms1g -Xmx${task.memory.toGiga()}g -Djava.io.tmpdir=." \
      HaplotypeCaller \
      -R ${reference} \
      -I ${bam} \
      -O ${sample_id}.gatk.vcf.gz
    gatk IndexFeatureFile -I ${sample_id}.gatk.vcf.gz
    """

    stub:
    """
    touch ${sample_id}.gatk.vcf.gz ${sample_id}.gatk.vcf.gz.tbi
    """
}

process CALL_DEEPVARIANT {
    tag "${sample_id}"
    publishDir "${params.outdir}/variants", mode: 'copy', pattern: "${sample_id}.deepvariant.vcf.gz*"

    input:
    tuple val(sample_id), path(bam), path(bai), path(reference), path(fai)

    output:
    tuple val(sample_id), path("${sample_id}.deepvariant.vcf.gz"), path("${sample_id}.deepvariant.vcf.gz.tbi"), emit: vcf

    script:
    """
    set -euo pipefail
    /opt/deepvariant/bin/run_deepvariant \
      --model_type=WGS \
      --ref=${reference} \
      --reads=${bam} \
      --output_vcf=${sample_id}.deepvariant.vcf.gz \
      --num_shards=${task.cpus}
    tabix -f -p vcf ${sample_id}.deepvariant.vcf.gz
    """

    stub:
    """
    touch ${sample_id}.deepvariant.vcf.gz ${sample_id}.deepvariant.vcf.gz.tbi
    """
}