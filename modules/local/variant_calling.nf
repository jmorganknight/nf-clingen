process CALL_GATK {
    tag "${sample_id}"
    publishDir "${params.outdir}/variants", mode: 'copy', pattern: "${sample_id}.gatk.raw.vcf.gz*"

    input:
    tuple val(sample_id), path(bam), path(bai), path(reference), path(fai), path(dict)

    output:
    tuple val(sample_id), path("${sample_id}.gatk.raw.vcf.gz"), path("${sample_id}.gatk.raw.vcf.gz.tbi"), emit: raw_vcf

    script:
    """
    set -euo pipefail
    gatk \
      --java-options "-Xms1g -Xmx${task.memory.toGiga()}g -Djava.io.tmpdir=." \
      HaplotypeCaller \
      -R ${reference} \
      -I ${bam} \
      -O ${sample_id}.gatk.raw.vcf.gz
    gatk IndexFeatureFile -I ${sample_id}.gatk.raw.vcf.gz
    """

    stub:
    """
    touch ${sample_id}.gatk.raw.vcf.gz ${sample_id}.gatk.raw.vcf.gz.tbi
    """
}

process FILTER_GATK {
    tag "${sample_id}"
    publishDir "${params.outdir}/variants", mode: 'copy', pattern: "${sample_id}.gatk.vcf.gz*"
    container 'broadinstitute/gatk:4.6.1.0'

    input:
    tuple val(sample_id), path(raw_vcf), path(raw_tbi), path(reference), path(fai), path(dict)

    output:
    tuple val(sample_id), path("${sample_id}.gatk.vcf.gz"), path("${sample_id}.gatk.vcf.gz.tbi"), emit: vcf

    script:
    """
    set -euo pipefail
    
    # Hard-filter SNPs and INDELs separately using GATK best-practice thresholds for exome
    # SNP filters: QD < 2.0, FS > 60.0, MQ < 40.0, MQRankSum < -12.5, ReadPosRankSum < -8.0
    # INDEL filters: QD < 2.0, FS > 200.0, ReadPosRankSum < -20.0
    
    # Split into SNPs and INDELs
    gatk SelectVariants -R ${reference} -V ${raw_vcf} --select-type SNP -O snps.vcf.gz
    gatk SelectVariants -R ${reference} -V ${raw_vcf} --select-type INDEL -O indels.vcf.gz
    
    # Filter SNPs
    gatk VariantFiltration \
      -R ${reference} \
      -V snps.vcf.gz \
      --filter-expression "QD < 2.0" --filter-name "LowQD_SNP" \
      --filter-expression "FS > 60.0" --filter-name "HighFS_SNP" \
      --filter-expression "MQ < 40.0" --filter-name "LowMQ_SNP" \
      --filter-expression "MQRankSum < -12.5" --filter-name "LowMQRankSum_SNP" \
      --filter-expression "ReadPosRankSum < -8.0" --filter-name "LowReadPosRankSum_SNP" \
      -O snps.filtered.vcf.gz
    
    # Filter INDELs
    gatk VariantFiltration \
      -R ${reference} \
      -V indels.vcf.gz \
      --filter-expression "QD < 2.0" --filter-name "LowQD_INDEL" \
      --filter-expression "FS > 200.0" --filter-name "HighFS_INDEL" \
      --filter-expression "ReadPosRankSum < -20.0" --filter-name "LowReadPosRankSum_INDEL" \
      -O indels.filtered.vcf.gz
    
    # Merge back together
    gatk MergeVcfs -I snps.filtered.vcf.gz -I indels.filtered.vcf.gz -O ${sample_id}.gatk.vcf.gz
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
    work_dir="\$(pwd)"
    dv_tmp="\${work_dir}/deepvariant_tmp"
    dv_vcf="\${work_dir}/${sample_id}.deepvariant.vcf.gz"
    mkdir -p "\${dv_tmp}"

    # DeepVariant can occasionally return a non-zero code even after writing outputs.
    # Capture the exit code and validate output existence explicitly.
    set +e
    /opt/deepvariant/bin/run_deepvariant \
      --model_type=WES \
      --ref=${reference} \
      --reads=${bam} \
      --output_vcf="\${dv_vcf}" \
      --intermediate_results_dir="\${dv_tmp}" \
      --num_shards=${task.cpus}
    dv_exit=\$?
    set -e

    # Fallback: if DeepVariant wrote output in a nested location, copy it into the task root.
    if [[ ! -s "\${dv_vcf}" ]]; then
      alt_vcf="\$(find "\${work_dir}" -maxdepth 4 -type f -name "${sample_id}.deepvariant.vcf.gz" | head -n 1 || true)"
      if [[ -n "\${alt_vcf}" && -s "\${alt_vcf}" ]]; then
        cp -f "\${alt_vcf}" "\${dv_vcf}"
      fi
    fi

    # Fail fast with a useful message if DeepVariant did not write the expected VCF.
    [[ -s "\${dv_vcf}" ]] || {
      echo "ERROR: DeepVariant did not produce ${sample_id}.deepvariant.vcf.gz" >&2
      echo "run_deepvariant exit code: \${dv_exit}" >&2
      echo "Diagnostics: listing task directory" >&2
      find "\${work_dir}" -maxdepth 3 -type f | sed -n '1,200p' >&2
      exit \${dv_exit}
    }

    if [[ \${dv_exit} -ne 0 ]]; then
      echo "WARN: run_deepvariant exited with code \${dv_exit}, but output VCF exists; continuing." >&2
    fi

    tabix -f -p vcf "\${dv_vcf}"
    """

    stub:
    """
    touch ${sample_id}.deepvariant.vcf.gz ${sample_id}.deepvariant.vcf.gz.tbi
    """
}

process VQSR_GATK {
    tag "${sample_id}"
    publishDir "${params.outdir}/variants", mode: 'copy', pattern: "${sample_id}.gatk.vcf.gz*"
    container 'broadinstitute/gatk:4.6.1.0'

    input:
    tuple val(sample_id), path(raw_vcf), path(raw_tbi), path(reference), path(fai), path(dict), path(vqsr_snp_resource, stageAs: 'vqsr_snp.vcf.gz'), path(vqsr_omni_resource, stageAs: 'vqsr_omni.vcf.gz'), path(vqsr_indel_resource, stageAs: 'vqsr_indel.vcf.gz'), path(vqsr_1kg_resource, stageAs: 'vqsr_1kg.vcf.gz'), path(vqsr_dbsnp_resource, stageAs: 'vqsr_dbsnp.vcf.gz')

    output:
    tuple val(sample_id), path("${sample_id}.gatk.vcf.gz"), path("${sample_id}.gatk.vcf.gz.tbi"), emit: vcf

    script:
    """
    set -euo pipefail

    # If any required training resources are missing, fall back to the input VCF.
    if [[ ! -s "${vqsr_snp_resource}" || ! -s "${vqsr_omni_resource}" || ! -s "${vqsr_indel_resource}" || ! -s "${vqsr_1kg_resource}" || ! -s "${vqsr_dbsnp_resource}" ]]; then
      cp ${raw_vcf} ${sample_id}.gatk.vcf.gz
      cp ${raw_tbi} ${sample_id}.gatk.vcf.gz.tbi
      exit 0
    fi

    # Nextflow stages only declared files; ensure each training resource has an index in the task workdir.
    gatk IndexFeatureFile -I ${vqsr_snp_resource} || true
    gatk IndexFeatureFile -I ${vqsr_omni_resource} || true
    gatk IndexFeatureFile -I ${vqsr_indel_resource} || true
    gatk IndexFeatureFile -I ${vqsr_1kg_resource} || true
    gatk IndexFeatureFile -I ${vqsr_dbsnp_resource} || true

    gatk SelectVariants -R ${reference} -V ${raw_vcf} --select-type SNP -O snps.vcf.gz
    gatk SelectVariants -R ${reference} -V ${raw_vcf} --select-type INDEL -O indels.vcf.gz

    gatk VariantRecalibrator \
      -R ${reference} \
      -V snps.vcf.gz \
      --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${vqsr_dbsnp_resource} \
      --resource:known,known=false,training=true,truth=false,prior=10.0 ${vqsr_1kg_resource} \
      --resource:omni,known=false,training=true,truth=false,prior=12.0 ${vqsr_omni_resource} \
      --resource:hapmap,known=false,training=true,truth=true,prior=15.0 ${vqsr_snp_resource} \
      -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
      --mode SNP \
      -O snps.recal \
      --tranches-file snps.tranches

    gatk ApplyVQSR \
      -R ${reference} \
      -V snps.vcf.gz \
      --recal-file snps.recal \
      --tranches-file snps.tranches \
      --truth-sensitivity-filter-level 99.0 \
      --mode SNP \
      -O snps.filtered.vcf.gz

    gatk VariantRecalibrator \
      -R ${reference} \
      -V indels.vcf.gz \
      --resource:mills,known=false,training=true,truth=true,prior=12.0 ${vqsr_indel_resource} \
      -an QD -an ReadPosRankSum -an FS -an SOR \
      --mode INDEL \
      -O indels.recal \
      --tranches-file indels.tranches

    gatk ApplyVQSR \
      -R ${reference} \
      -V indels.vcf.gz \
      --recal-file indels.recal \
      --tranches-file indels.tranches \
      --truth-sensitivity-filter-level 99.0 \
      --mode INDEL \
      -O indels.filtered.vcf.gz

    gatk MergeVcfs -I snps.filtered.vcf.gz -I indels.filtered.vcf.gz -O ${sample_id}.gatk.vcf.gz
    gatk IndexFeatureFile -I ${sample_id}.gatk.vcf.gz
    """

    stub:
    """
    touch ${sample_id}.gatk.vcf.gz ${sample_id}.gatk.vcf.gz.tbi
    """
}