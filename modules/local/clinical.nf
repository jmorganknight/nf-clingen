process INDEX_CLINVAR_RESOURCE {
    tag 'clinvar-resource'
    publishDir "${params.outdir}/clinical/resources", mode: 'copy', pattern: 'clinvar.vcf.gz*'

    input:
    path(clinvar_vcf)

    output:
    tuple path('clinvar.vcf.gz'), path('clinvar.vcf.gz.tbi'), emit: resource

    script:
    def isGz = clinvar_vcf.name.endsWith('.gz')
    """
    set -euo pipefail
    if [[ '${isGz}' == 'true' ]]; then
      cat ${clinvar_vcf} > clinvar.vcf.gz
    else
      bgzip -c ${clinvar_vcf} > clinvar.vcf.gz
    fi
    tabix -f -p vcf clinvar.vcf.gz
    """

    stub:
    """
    touch clinvar.vcf.gz clinvar.vcf.gz.tbi
    """
}

process RUN_CLINICAL_ANNOTATION {
    tag "${sample_id}"
    publishDir "${params.outdir}/clinical/annotations", mode: 'copy'

    input:
    tuple val(sample_id), path(vcf), path(vcf_index)

    output:
    tuple val(sample_id), path("${sample_id}.clinical.annotated.vcf.gz"), path("${sample_id}.clinical.annotated.vcf.gz.tbi"), path("${sample_id}.clinical_annotation.tsv"), path("${sample_id}.clinical_summary.json"), emit: annotation

    script:
    """
    set -euo pipefail
    ln -sf ${vcf} ${sample_id}.clinical.annotated.vcf.gz
    ln -sf ${vcf_index} ${sample_id}.clinical.annotated.vcf.gz.tbi

    # Build query fields dynamically so missing INFO tags do not crash annotation.
    clnsig_expr='.'
    clndn_expr='.'
    geneinfo_expr='.'
    ann_expr='.'
    vcf_header=\$(bcftools view -h ${sample_id}.clinical.annotated.vcf.gz)
    if grep -q '^##INFO=<ID=CLNSIG,' <<< "\${vcf_header}"; then
      clnsig_expr='%INFO/CLNSIG'
    fi
    if grep -q '^##INFO=<ID=CLNDN,' <<< "\${vcf_header}"; then
      clndn_expr='%INFO/CLNDN'
    fi
    if grep -q '^##INFO=<ID=GENEINFO,' <<< "\${vcf_header}"; then
      geneinfo_expr='%INFO/GENEINFO'
    fi
    if grep -q '^##INFO=<ID=ANN,' <<< "\${vcf_header}"; then
      ann_expr='%INFO/ANN'
    fi

    bcftools query \
      -f "%CHROM\t%POS\t%ID\t%REF\t%ALT\t%FILTER[\t%GT]\t\${clnsig_expr}\t\${clndn_expr}\t\${geneinfo_expr}\t\${ann_expr}\n" \
      ${sample_id}.clinical.annotated.vcf.gz > ${sample_id}.clinical.raw.tsv

    awk 'BEGIN { OFS="\t"; print "rank","chrom","pos","id","ref","alt","filter","genotype","severity","clnsig","disease","gene","annotation" }
         {
             genotype = (\$7 == "" ? "." : \$7)
             if (genotype == "0/0" || genotype == "0|0" || genotype == "./." || genotype == ".|.") next
             clnsig = toupper(\$8)
             disease = (\$9 == "" ? "." : \$9)
             gene = (\$10 == "" ? "." : \$10)
             ann = toupper(\$11)
             severity = "review"
             rank = 4
         if (clnsig ~ /CONFLICTING_CLASSIFICATIONS/) { severity = "review"; rank = 4 }
         else if (clnsig ~ /PATHOGENIC/ && clnsig ~ /LIKELY_PATHOGENIC/) { severity = "high"; rank = 2 }
         else if (clnsig ~ /(^|[|,; ])PATHOGENIC([|,; ]|\$)/) { severity = "critical"; rank = 1 }
             else if (ann ~ /HIGH|STOP_GAINED|FRAMESHIFT|SPLICE/ || clnsig ~ /LIKELY_PATHOGENIC/) { severity = "high"; rank = 2 }
             else if (\$6 == "PASS" || \$6 == ".") { severity = "moderate"; rank = 3 }
             print rank,\$1,\$2,\$3,\$4,\$5,\$6,genotype,severity,(\$8 == "" ? "." : \$8),disease,gene,(\$11 == "" ? "." : \$11)
         }' ${sample_id}.clinical.raw.tsv > ${sample_id}.clinical.scored.tsv

    tail -n +2 ${sample_id}.clinical.scored.tsv | sort -t \$'\t' -k1,1n -k2,2 -k3,3n > ${sample_id}.clinical.sorted.tsv
    {
      head -n 1 ${sample_id}.clinical.scored.tsv | cut -f2-
      head -n ${params.max_report_variants} ${sample_id}.clinical.sorted.tsv | cut -f2-
    } > ${sample_id}.clinical_annotation.tsv

    total_scored=0
    prioritized_count=0
    critical_count=0
    high_count=0
    moderate_count=0
    if [[ -s ${sample_id}.clinical.scored.tsv ]]; then
      total_scored=\$(tail -n +2 ${sample_id}.clinical.scored.tsv | wc -l | tr -d ' ')
      prioritized_count=\$(tail -n +2 ${sample_id}.clinical_annotation.tsv | wc -l | tr -d ' ')
      critical_count=\$(awk -F '\t' 'NR>1 && \$8=="critical" {c++} END {print c+0}' ${sample_id}.clinical_annotation.tsv)
      high_count=\$(awk -F '\t' 'NR>1 && \$8=="high" {c++} END {print c+0}' ${sample_id}.clinical_annotation.tsv)
      moderate_count=\$(awk -F '\t' 'NR>1 && \$8=="moderate" {c++} END {print c+0}' ${sample_id}.clinical_annotation.tsv)
    fi

    printf '{\n  "sample_id": "%s",\n  "source": "intrinsic_vcf_annotations",\n  "prioritized_variants": %s,\n  "critical": %s,\n  "high": %s,\n  "moderate": %s\n}\n' \
      "${sample_id}" "\${prioritized_count}" "\${critical_count}" "\${high_count}" "\${moderate_count}" > ${sample_id}.clinical_summary.json
    """

    stub:
    """
    touch ${sample_id}.clinical.annotated.vcf.gz ${sample_id}.clinical.annotated.vcf.gz.tbi
    echo -e "chrom\tpos\tid\tref" > ${sample_id}.clinical_annotation.tsv
    echo '{"sample_id": "${sample_id}", "prioritized_variants": 0}' > ${sample_id}.clinical_summary.json
    """
}

process RUN_CLINVAR_ANNOTATION {
    tag "${sample_id}"
    publishDir "${params.outdir}/clinical/annotations", mode: 'copy'

    input:
    tuple val(sample_id), path(vcf), path(vcf_index), path(clinvar_vcf), path(clinvar_tbi)

    output:
    tuple val(sample_id), path("${sample_id}.clinical.annotated.vcf.gz"), path("${sample_id}.clinical.annotated.vcf.gz.tbi"), path("${sample_id}.clinical_annotation.tsv"), path("${sample_id}.clinical_summary.json"), emit: annotation

    script:
    """
    set -euo pipefail
    # Debug: list input files
    echo "DEBUG: Provided input files:"
    ls -lh ${vcf} ${vcf_index} ${clinvar_vcf} ${clinvar_tbi} 2>&1 || true
    
    # Link input files into work directory for consistent access within container
    ln -sf ${vcf} input.vcf.gz
    ln -sf ${vcf_index} input.vcf.gz.tbi
    ln -sf ${clinvar_vcf} clinvar_input.vcf.gz
    ln -sf ${clinvar_tbi} clinvar_input.vcf.gz.tbi
    
    echo "DEBUG: Linked files created:"
    ls -lh input.* clinvar_input.* 2>&1 || true
    
    # Test bcftools access
    echo "DEBUG: Testing bcftools access..."
    bcftools --version || { echo "FATAL: bcftools not available"; exit 1; }
    
    query_contig=\$(bcftools view -H input.vcf.gz 2>&1 | head -n 1 | cut -f1 || echo "ERROR")
    clinvar_contig=\$(bcftools view -H clinvar_input.vcf.gz 2>&1 | head -n 1 | cut -f1 || echo "ERROR")
    
    if [[ "\${query_contig}" == "ERROR" ]] || [[ "\${clinvar_contig}" == "ERROR" ]]; then
      echo "FATAL: Failed to read VCF headers."
      echo "query_contig=\${query_contig}"
      echo "clinvar_contig=\${clinvar_contig}"
      exit 1
    fi
    
    echo "DEBUG: Contigs detected - query=\${query_contig}, clinvar=\${clinvar_contig}"

    clinvar_for_annot=clinvar_input.vcf.gz
    if [[ "\${query_contig}" == chr* && "\${clinvar_contig}" != chr* ]]; then
      echo "DEBUG: Chromosome renaming needed (query=\${query_contig}, clinvar=\${clinvar_contig})"
      : > rename_chrs.tsv
      for c in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y MT; do
        printf '%s\tchr%s\n' "\${c}" "\${c}" >> rename_chrs.tsv
      done
      echo "DEBUG: Running bcftools annotate --rename-chrs..."
      bcftools annotate --rename-chrs rename_chrs.tsv -Oz -o clinvar.norm.vcf.gz clinvar_input.vcf.gz || \
        { echo "FATAL: bcftools annotate rename failed"; exit 1; }
      echo "DEBUG: Running tabix..."
      tabix -f -p vcf clinvar.norm.vcf.gz || \
        { echo "FATAL: tabix failed"; exit 1; }
      clinvar_for_annot=clinvar.norm.vcf.gz
    else
      echo "DEBUG: No chromosome renaming needed"
    fi

    echo "DEBUG: Running main bcftools annotate..."
    bcftools annotate \
      -a \${clinvar_for_annot} \
      -c ID,INFO/CLNSIG,INFO/CLNDN,INFO/GENEINFO \
      -Oz \
      -o ${sample_id}.clinical.annotated.vcf.gz \
      input.vcf.gz || \
      { echo "FATAL: bcftools annotate failed"; exit 1; }
    
    echo "DEBUG: Running final tabix..."
    tabix -f -p vcf ${sample_id}.clinical.annotated.vcf.gz || \
      { echo "FATAL: final tabix failed"; exit 1; }
    
    echo "DEBUG: Files created after annotation:"
    ls -lh ${sample_id}.clinical.annotated.vcf.gz* 2>&1 || true

    # Build query fields dynamically so missing INFO tags do not crash annotation.
    clnsig_expr='.'
    clndn_expr='.'
    geneinfo_expr='.'
    ann_expr='.'
    vcf_header=\$(bcftools view -h ${sample_id}.clinical.annotated.vcf.gz)
    if grep -q '^##INFO=<ID=CLNSIG,' <<< "\${vcf_header}"; then
      clnsig_expr='%INFO/CLNSIG'
    fi
    if grep -q '^##INFO=<ID=CLNDN,' <<< "\${vcf_header}"; then
      clndn_expr='%INFO/CLNDN'
    fi
    if grep -q '^##INFO=<ID=GENEINFO,' <<< "\${vcf_header}"; then
      geneinfo_expr='%INFO/GENEINFO'
    fi
    if grep -q '^##INFO=<ID=ANN,' <<< "\${vcf_header}"; then
      ann_expr='%INFO/ANN'
    fi

    bcftools query \
      -f "%CHROM\t%POS\t%ID\t%REF\t%ALT\t%FILTER[\t%GT]\t\${clnsig_expr}\t\${clndn_expr}\t\${geneinfo_expr}\t\${ann_expr}\n" \
      ${sample_id}.clinical.annotated.vcf.gz > ${sample_id}.clinical.raw.tsv

    awk 'BEGIN { OFS="\t"; print "rank","chrom","pos","id","ref","alt","filter","genotype","severity","clnsig","disease","gene","annotation" }
         {
             genotype = (\$7 == "" ? "." : \$7)
             if (genotype == "0/0" || genotype == "0|0" || genotype == "./." || genotype == ".|.") next
             clnsig = toupper(\$8)
             disease = (\$9 == "" ? "." : \$9)
             gene = (\$10 == "" ? "." : \$10)
             ann = toupper(\$11)
             severity = "review"
             rank = 4
         if (clnsig ~ /CONFLICTING_CLASSIFICATIONS/) { severity = "review"; rank = 4 }
         else if (clnsig ~ /PATHOGENIC/ && clnsig ~ /LIKELY_PATHOGENIC/) { severity = "high"; rank = 2 }
         else if (clnsig ~ /(^|[|,; ])PATHOGENIC([|,; ]|\$)/) { severity = "critical"; rank = 1 }
             else if (ann ~ /HIGH|STOP_GAINED|FRAMESHIFT|SPLICE/ || clnsig ~ /LIKELY_PATHOGENIC/) { severity = "high"; rank = 2 }
             else if (\$6 == "PASS" || \$6 == ".") { severity = "moderate"; rank = 3 }
             print rank,\$1,\$2,\$3,\$4,\$5,\$6,genotype,severity,(\$8 == "" ? "." : \$8),disease,gene,(\$11 == "" ? "." : \$11)
         }' ${sample_id}.clinical.raw.tsv > ${sample_id}.clinical.scored.tsv

    tail -n +2 ${sample_id}.clinical.scored.tsv | sort -t \$'\t' -k1,1n -k2,2 -k3,3n > ${sample_id}.clinical.sorted.tsv
    {
      head -n 1 ${sample_id}.clinical.scored.tsv | cut -f2-
      head -n ${params.max_report_variants} ${sample_id}.clinical.sorted.tsv | cut -f2-
    } > ${sample_id}.clinical_annotation.tsv

    total_scored=0
    prioritized_count=0
    critical_count=0
    high_count=0
    moderate_count=0
    if [[ -s ${sample_id}.clinical.scored.tsv ]]; then
      total_scored=\$(tail -n +2 ${sample_id}.clinical.scored.tsv | wc -l | tr -d ' ')
      prioritized_count=\$(tail -n +2 ${sample_id}.clinical_annotation.tsv | wc -l | tr -d ' ')
      critical_count=\$(awk -F '\t' 'NR>1 && \$8=="critical" {c++} END {print c+0}' ${sample_id}.clinical_annotation.tsv)
      high_count=\$(awk -F '\t' 'NR>1 && \$8=="high" {c++} END {print c+0}' ${sample_id}.clinical_annotation.tsv)
      moderate_count=\$(awk -F '\t' 'NR>1 && \$8=="moderate" {c++} END {print c+0}' ${sample_id}.clinical_annotation.tsv)
    fi

    printf '{\n  "sample_id": "%s",\n  "source": "clinvar_overlay",\n  "prioritized_variants": %s,\n  "critical": %s,\n  "high": %s,\n  "moderate": %s\n}\n' \
      "${sample_id}" "\${prioritized_count}" "\${critical_count}" "\${high_count}" "\${moderate_count}" > ${sample_id}.clinical_summary.json
    """

    stub:
    """
    touch ${sample_id}.clinical.annotated.vcf.gz ${sample_id}.clinical.annotated.vcf.gz.tbi
    echo -e "chrom\tpos\tid\tref" > ${sample_id}.clinical_annotation.tsv
    echo '{"sample_id": "${sample_id}", "prioritized_variants": 0}' > ${sample_id}.clinical_summary.json
    """
}

process BUILD_CLINICAL_REPORT {
    tag "${sample_id}"
    publishDir "${params.outdir}/clinical/reports", mode: 'copy'

    input:
    tuple val(sample_id), path(vcf), path(vcf_index), path(annotation_tsv), path(summary_json), path(report_builder_script)

    output:
    tuple val(sample_id), path("${sample_id}.critical_variants.tsv"), path("${sample_id}.clinical_report.html"), path("${sample_id}.clinical_report.pdf"), emit: deliverables

    script:
    """
    set -euo pipefail

    apt-get update
    apt-get install -y --no-install-recommends \
      procps \
      libcairo2 \
      libpango-1.0-0 \
      libpangocairo-1.0-0 \
      libgdk-pixbuf-2.0-0 \
      libffi8 \
      libjpeg62-turbo \
      libopenjp2-7 \
      libxml2 \
      libxslt1.1 \
      shared-mime-info \
      fonts-dejavu-core \
      ca-certificates
    rm -rf /var/lib/apt/lists/*

    python -m pip install --no-cache-dir weasyprint==61.2 pydyf==0.10.0

    python ${report_builder_script} \
      --sample-id ${sample_id} \
      --annotation-tsv ${annotation_tsv} \
      --summary-json ${summary_json} \
      --annotated-vcf ${vcf} \
      --aligner "${params.aligner}" \
      --preprocess "${params.preprocess}" \
      --caller "${params.caller}"

    weasyprint ${sample_id}.clinical_report.html ${sample_id}.clinical_report.pdf
    """

    stub:
    """
    touch ${sample_id}.critical_variants.tsv ${sample_id}.clinical_report.html ${sample_id}.clinical_report.pdf
    """
}
