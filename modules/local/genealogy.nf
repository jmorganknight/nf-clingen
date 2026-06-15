// Phase variants with Eagle2 then impute with Beagle 5.4.
// Required inputs:
//   params.eagle2_genetic_map  — path to Eagle2 genetic map directory (e.g. Eagle_v2.4.1/tables/)
//   params.beagle_ref_panel    — path to per-chromosome Beagle reference panels,
//                                either VCFs (*.chr{CHROM}.vcf.gz) or bref3 (*.chr{CHROM}.bref3)
//   params.beagle_genetic_map  — path to per-chromosome Beagle genetic maps, named *.chr{CHROM}.map.gz (optional; improves accuracy)
//   params.beagle_jar          — path to beagle.*.jar (defaults to container path)
// Container images used:
//   phase  : indapa/indapa-eagle:latest
//   impute : quay.io/biocontainers/beagle:5.4_22Jul22.46e--hdfd78af_0
process PHASE_EAGLE2 {
    tag "${sample_id}"
    container 'indapa/indapa-eagle:latest'
    publishDir "${params.outdir}/genealogy/phased", mode: 'copy'

    input:
    tuple val(sample_id), path(vcf), path(vcf_index)
    path(genetic_map_dir)
    path(ref_panel_dir)

    output:
    tuple val(sample_id), path("${sample_id}.phased.vcf.gz"), path("${sample_id}.phased.vcf.gz.tbi"), emit: phased

    script:
    """
    set -euo pipefail

    # Extract chromosomes present in this VCF using portable numeric ordering (1..22,X,Y)
    chroms=\$(bcftools view -h ${vcf} | grep '^##contig' | sed -n 's/.*ID=\\([^,>]*\\).*/\\1/p' | grep -E '^(chr)?[0-9XY]+\$' | awk '{ c=\$0; sub(/^chr/,"",c); if(c=="X")o=23; else if(c=="Y")o=24; else o=c+0; print o"\\t"\$0 }' | sort -n | cut -f2 | tr '\\n' ' ')

    # Phase each chromosome separately then merge
    phased_vcfs=""
    for chrom in \$chroms; do
        # Locate per-chromosome reference panel — prefer hg38, fall back to any match
        ref_vcf=\$(find ${ref_panel_dir} -name "*hg38*.\${chrom}.vcf.gz" 2>/dev/null | sort | head -1)
        if [[ -z "\${ref_vcf}" ]]; then ref_vcf=\$(find ${ref_panel_dir} -name "*.\${chrom}.vcf.gz" 2>/dev/null | sort | head -1); fi
        if [[ -z "\${ref_vcf}" ]]; then
            echo "WARNING: no ref panel for \${chrom}; skipping Eagle2 phasing for this chromosome" >&2
            continue
        fi

        /opt/Eagle_v2.4.1/eagle \\
            --vcfTarget=${vcf} \\
            --vcfRef=\${ref_vcf} \\
            --geneticMapFile=${genetic_map_dir}/genetic_map_hg38_withX.txt.gz \\
            --chrom=\${chrom} \\
            --outPrefix=${sample_id}.\${chrom}.eagle2 \\
            --numThreads=${task.cpus} \\
            2>&1 | tee -a eagle2_\${chrom}.log

        bcftools index --tbi ${sample_id}.\${chrom}.eagle2.vcf.gz
        phased_vcfs="\${phased_vcfs} ${sample_id}.\${chrom}.eagle2.vcf.gz"
    done

    if [[ -z "\${phased_vcfs// }" ]]; then
        echo "ERROR: no chromosomes were phased; check contig naming and reference panel availability" >&2
        exit 1
    fi

    # Merge per-chromosome phased VCFs back into a single file
    bcftools concat --allow-overlaps --rm-dups all \\
        \${phased_vcfs} \\
        -Oz -o ${sample_id}.phased.vcf.gz
    bcftools index --tbi ${sample_id}.phased.vcf.gz
    """

    stub:
    """
    touch ${sample_id}.phased.vcf.gz ${sample_id}.phased.vcf.gz.tbi
    """
}

process IMPUTE_BEAGLE {
    tag "${sample_id}"
    container 'quay.io/biocontainers/beagle:5.4_22Jul22.46e--hdfd78af_0'
    publishDir "${params.outdir}/genealogy/imputed", mode: 'copy', pattern: '*.imputed.vcf.gz'

    input:
    tuple val(sample_id), path(phased_vcf), path(phased_vcf_index)
    path(ref_panel_dir)
    path(beagle_map_dir)

    output:
    tuple val(sample_id), path('*.imputed.vcf.gz'), emit: chrom_vcfs

    script:
    // Beagle jar location: prefer params, fall back to container path
    def beagle_jar = params.beagle_jar ?: '/usr/local/share/beagle-5.4_22Jul22.46e-0/beagle.jar'
    """
    set -euo pipefail

    # Extract chromosome list from VCF header in portable order (1..22,X,Y)
    chroms=\$(zcat ${phased_vcf} | grep '^##contig' | sed -n 's/.*ID=\\([^,>]*\\).*/\\1/p' | grep -E '^(chr)?[0-9XY]+\$' | awk '{ c=\$0; sub(/^chr/,"",c); if(c=="X")o=23; else if(c=="Y")o=24; else o=c+0; print o"\\t"\$0 }' | sort -n | cut -f2)

    if [[ -z "\${chroms}" ]]; then
        echo "WARNING: no ##contig headers; scanning records for chromosome names" >&2
        chroms=\$(zcat ${phased_vcf} | grep -v '^#' | cut -f1 | grep -E '^(chr)?[0-9XY]+\$' | awk '{ c=\$0; sub(/^chr/,"",c); if(c=="X")o=23; else if(c=="Y")o=24; else o=c+0; print o"\\t"\$0 }' | sort -n | cut -f2 | awk '!seen[\$0]++')
    fi

    for chrom in \${chroms}; do
        # Prepare chromosome-specific target VCF and remove duplicate markers
        # by key (CHROM,POS,REF,ALT), keeping first occurrence.
        chrom_target_dedup="${sample_id}.\${chrom}.target.dedup.vcf.gz"
        zcat "${phased_vcf}" | awk -v C="\${chrom}" 'BEGIN{FS=OFS="\t"} /^#/{print; next} \$1==C{key=\$1":"\$2":"\$4":"\$5; if(!seen[key]++){print}}' | gzip > "\${chrom_target_dedup}"

      # Skip chromosomes that have no variant records in the target VCF.
      # Beagle errors on header-only VCFs (common for chrY in XX samples).
      if [[ \$(zcat "\${chrom_target_dedup}" | grep -vc '^#') -eq 0 ]]; then
        echo "INFO: no target variants for \${chrom}; skipping imputation for this chromosome" >&2
        continue
      fi

        # Locate reference panel — prefer hg38 VCF, then hg38 bref3, then any build
        ref=\$(find ${ref_panel_dir} -name "*hg38*.\${chrom}.vcf.gz" 2>/dev/null | sort | head -1)
        if [[ -z "\${ref}" ]]; then ref=\$(find ${ref_panel_dir} -name "*hg38*.\${chrom}.bref3" 2>/dev/null | sort | head -1); fi
        if [[ -z "\${ref}" ]]; then ref=\$(find ${ref_panel_dir} -name "*.\${chrom}.vcf.gz" 2>/dev/null | sort | head -1); fi
        if [[ -z "\${ref}" ]]; then ref=\$(find ${ref_panel_dir} -name "*.\${chrom}.bref3" 2>/dev/null | sort | head -1); fi

        # If the reference is a VCF, drop duplicate markers by key
        # (CHROM,POS,REF,ALT), keeping first occurrence.
        ref_arg=""
        if [[ -n "\${ref}" && "\${ref}" == *.vcf.gz ]]; then
            ref_dedup="${sample_id}.\${chrom}.ref.dedup.vcf.gz"
            zcat "\${ref}" | awk -v C="\${chrom}" 'BEGIN{FS=OFS="\t"} /^#/{print; next} \$1==C{key=\$1":"\$2":"\$4":"\$5; if(!seen[key]++){print}}' | gzip > "\${ref_dedup}"
            ref_arg="ref=\${ref_dedup}"
        elif [[ -n "\${ref}" ]]; then
            ref_arg="ref=\${ref}"
        fi

        # Locate genetic map
        map_arg=""
        map_file=\$(find ${beagle_map_dir} -name "*.\${chrom}.map.gz" 2>/dev/null | sort | head -1)
        # Beagle map files may use bare chrom names (e.g. "1") while VCFs use "chr1".
        # Re-prefix on the fly so Beagle can match chromosomes correctly.
        if [[ -n "\${map_file}" ]]; then
            zcat "\${map_file}" | awk 'BEGIN{OFS="\t"} \$1!~/^chr/{sub(/^/,"chr",\$1)} {print}' | gzip > tmp_map_\${chrom}.gz
            map_arg="map=tmp_map_\${chrom}.gz"
        fi

        if [[ -n "\${ref_arg}" ]]; then
            java -Xmx\${JAVA_HEAP:-16g} -jar ${beagle_jar} \\
            gt=\${chrom_target_dedup} \\
            \${ref_arg} \\
                chrom=\${chrom} \\
                \${map_arg} \\
                out=${sample_id}.\${chrom}.imputed \\
                nthreads=${task.cpus} \\
                2>&1 | tee -a beagle_\${chrom}.log
        else
            echo "WARNING: no ref panel for \${chrom}; running phasing refinement only" >&2
            java -Xmx\${JAVA_HEAP:-16g} -jar ${beagle_jar} \\
            gt=\${chrom_target_dedup} \\
                chrom=\${chrom} \\
                \${map_arg} \\
                out=${sample_id}.\${chrom}.imputed \\
                nthreads=${task.cpus} \\
                2>&1 | tee -a beagle_\${chrom}.log
        fi
    done
    """

    stub:
    """
    touch ${sample_id}.chr1.imputed.vcf.gz ${sample_id}.chr22.imputed.vcf.gz
    """
}

process CONCAT_IMPUTED {
    tag "${sample_id}"
    publishDir "${params.outdir}/genealogy/imputed", mode: 'copy'

    input:
    tuple val(sample_id), path(chrom_vcfs)

    output:
    tuple val(sample_id), path("${sample_id}.imputed.vcf.gz"), path("${sample_id}.imputed.vcf.gz.tbi"), emit: imputed

    script:
    """
    set -euo pipefail

    # Build deterministic chromosome order from Nextflow-provided inputs.
    # Avoid relying on runtime globbing, which can be fragile across task setups.
    ordered_files=""
    for c in \$(seq 1 22) X Y; do
      for f in ${chrom_vcfs}; do
        base=\$(basename "\${f}")
        if [[ "\${base}" == *.chr\${c}.imputed.vcf.gz || "\${base}" == *.\${c}.imputed.vcf.gz ]]; then
          ordered_files="\${ordered_files} \${f}"
        fi
        done
    done

    # Append any remaining files not captured above (defensive)
    for f in ${chrom_vcfs}; do
        [[ " \${ordered_files} " == *" \${f} "* ]] || ordered_files="\${ordered_files} \${f}"
    done

    if [[ -z "\${ordered_files// }" ]]; then
        echo "ERROR: no per-chromosome imputed VCFs found for concat" >&2
        exit 1
    fi

    for f in \${ordered_files}; do
        bcftools index --tbi "\${f}"
    done

    bcftools concat --allow-overlaps --rm-dups all \\
        \${ordered_files} \\
        -Oz -o ${sample_id}.imputed.vcf.gz
    bcftools index --tbi ${sample_id}.imputed.vcf.gz
    """

    stub:
    """
    touch ${sample_id}.imputed.vcf.gz ${sample_id}.imputed.vcf.gz.tbi
    """
}

process BUILD_GENEALOGY_MANIFEST {
    tag "${sample_id}"
    publishDir "${params.outdir}/genealogy", mode: 'copy'

    input:
    tuple val(sample_id), path(phased_vcf), path(phased_vcf_index), path(imputed_vcf), path(imputed_vcf_index)

    output:
    tuple val(sample_id), path("${sample_id}.genealogy.manifest.txt"), emit: manifest
    tuple val(sample_id), path("${sample_id}.genealogy.report.html"), emit: report
    tuple val(sample_id), path("${sample_id}.genealogy.audit.json"), emit: audit

    script:
    """
    set -euo pipefail

    phased_variants=\$(bcftools stats ${phased_vcf} | awk '/^SN.*number of SNPs/ {print \$NF}')
    imputed_variants=\$(bcftools stats ${imputed_vcf} | awk '/^SN.*number of SNPs/ {print \$NF}')
    phased_records=\$(bcftools view -H ${phased_vcf} | wc -l)
    imputed_records=\$(bcftools view -H ${imputed_vcf} | wc -l)

    phased_sha256=\$(sha256sum ${phased_vcf} | awk '{print \$1}')
    phased_tbi_sha256=\$(sha256sum ${phased_vcf_index} | awk '{print \$1}')
    imputed_sha256=\$(sha256sum ${imputed_vcf} | awk '{print \$1}')
    imputed_tbi_sha256=\$(sha256sum ${imputed_vcf_index} | awk '{print \$1}')

    phased_size=\$(stat -c%s ${phased_vcf})
    imputed_size=\$(stat -c%s ${imputed_vcf})

    phased_has_index=\$([[ -s ${phased_vcf_index} ]] && echo "PASS" || echo "FAIL")
    imputed_has_index=\$([[ -s ${imputed_vcf_index} ]] && echo "PASS" || echo "FAIL")
    imputed_growth=\$([[ \${imputed_records} -ge \${phased_records} ]] && echo "PASS" || echo "REVIEW")

    run_date_utc=\$(date -u +%Y-%m-%dT%H:%M:%SZ)

    cat > ${sample_id}.genealogy.manifest.txt <<EOF_MAN
sample_id=${sample_id}
phased_vcf=${phased_vcf}
imputed_vcf=${imputed_vcf}
phase_engine=Eagle2
imputation_engine=Beagle 5.4
phased_snps=\${phased_variants}
imputed_snps=\${imputed_variants}
eagle2_genetic_map=${params.eagle2_genetic_map}
beagle_ref_panel=${params.beagle_ref_panel}
nf_clingen_version=${workflow.manifest.version}
run_name=${workflow.runName}
run_session_id=${workflow.sessionId}
run_date=\${run_date_utc}
EOF_MAN

    cat > ${sample_id}.genealogy.audit.json <<EOF_AUDIT
{
  "sample_id": "${sample_id}",
  "run": {
    "workflow_name": "${workflow.runName}",
    "session_id": "${workflow.sessionId}",
    "pipeline_version": "${workflow.manifest.version}",
    "generated_at_utc": "\${run_date_utc}"
  },
  "inputs": {
    "reference_fasta": "${params.reference}",
    "eagle2_genetic_map": "${params.eagle2_genetic_map}",
    "beagle_reference_panel": "${params.beagle_ref_panel}",
    "beagle_genetic_map": "${params.beagle_genetic_map ?: 'not_provided'}"
  },
  "metrics": {
    "phased_snp_count": \${phased_variants},
    "imputed_snp_count": \${imputed_variants},
    "phased_record_count": \${phased_records},
    "imputed_record_count": \${imputed_records},
    "phased_size_bytes": \${phased_size},
    "imputed_size_bytes": \${imputed_size}
  },
  "hashes": {
    "phased_vcf_sha256": "\${phased_sha256}",
    "phased_index_sha256": "\${phased_tbi_sha256}",
    "imputed_vcf_sha256": "\${imputed_sha256}",
    "imputed_index_sha256": "\${imputed_tbi_sha256}"
  },
  "quality_checks": {
    "phased_index_present": "\${phased_has_index}",
    "imputed_index_present": "\${imputed_has_index}",
    "imputed_records_gte_phased": "\${imputed_growth}"
  },
  "intended_use": "Research use only. Not a standalone diagnostic report."
}
EOF_AUDIT

    cat > ${sample_id}.genealogy.report.html <<EOF_HTML
<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1" />
  <title>Genealogy Imputation Report - ${sample_id}</title>
  <style>
    body { margin:0; background:#f4f7fb; color:#1d2a3a; font:15px/1.5 "Segoe UI", Tahoma, sans-serif; }
    .wrap { max-width:980px; margin:0 auto; padding:24px 18px 36px; }
    .card { background:#fff; border:1px solid #d7e0ea; border-radius:12px; padding:16px 18px; margin-bottom:14px; }
    h1 { margin:0 0 6px; color:#083b66; }
    h2 { margin:0 0 8px; color:#0f4b7f; font-size:18px; }
    .muted { color:#5a6b80; }
    table { width:100%; border-collapse:collapse; }
    th, td { border-bottom:1px solid #d7e0ea; padding:8px 6px; text-align:left; vertical-align:top; }
    th { width:34%; color:#35506d; }
    code { background:#f3f7fb; border:1px solid #dce7f2; border-radius:5px; padding:1px 5px; font-size:12px; }
  </style>
</head>
<body>
  <div class="wrap">
    <div class="card">
      <h1>Genealogy Imputation Report</h1>
      <p class="muted">Clinician-facing summary with run traceability and audit metrics.</p>
      <p><b>Sample:</b> ${sample_id}</p>
      <p><b>Run:</b> ${workflow.runName} | <b>Session:</b> ${workflow.sessionId} | <b>Generated (UTC):</b> \${run_date_utc}</p>
    </div>

    <div class="card">
      <h2>Interpretation Summary</h2>
      <p>Population phasing (Eagle2) and genotype imputation (Beagle 5.4) completed for genealogy-oriented downstream analysis.</p>
      <p>Imputed record count is <b>\${imputed_records}</b> versus phased record count <b>\${phased_records}</b>.</p>
      <p class="muted">Research use only. Not a standalone diagnostic medical interpretation.</p>
    </div>

    <div class="card">
      <h2>Core Metrics</h2>
      <table>
        <tr><th>Phased SNP Count</th><td>\${phased_variants}</td></tr>
        <tr><th>Imputed SNP Count</th><td>\${imputed_variants}</td></tr>
        <tr><th>Phased Record Count</th><td>\${phased_records}</td></tr>
        <tr><th>Imputed Record Count</th><td>\${imputed_records}</td></tr>
        <tr><th>Phased Index Present</th><td>\${phased_has_index}</td></tr>
        <tr><th>Imputed Index Present</th><td>\${imputed_has_index}</td></tr>
        <tr><th>Imputed >= Phased Records</th><td>\${imputed_growth}</td></tr>
      </table>
    </div>

    <div class="card">
      <h2>Traceability</h2>
      <table>
        <tr><th>Reference FASTA</th><td><code>${params.reference}</code></td></tr>
        <tr><th>Eagle2 Genetic Map</th><td><code>${params.eagle2_genetic_map}</code></td></tr>
        <tr><th>Beagle Reference Panel</th><td><code>${params.beagle_ref_panel}</code></td></tr>
        <tr><th>Beagle Genetic Map</th><td><code>${params.beagle_genetic_map ?: 'not_provided'}</code></td></tr>
        <tr><th>Phased VCF SHA-256</th><td><code>\${phased_sha256}</code></td></tr>
        <tr><th>Imputed VCF SHA-256</th><td><code>\${imputed_sha256}</code></td></tr>
        <tr><th>Phased Index SHA-256</th><td><code>\${phased_tbi_sha256}</code></td></tr>
        <tr><th>Imputed Index SHA-256</th><td><code>\${imputed_tbi_sha256}</code></td></tr>
      </table>
      <p class="muted">Machine-readable audit file: <code>${sample_id}.genealogy.audit.json</code></p>
    </div>
  </div>
</body>
</html>
EOF_HTML
    """

    stub:
    """
    touch ${sample_id}.genealogy.manifest.txt
    touch ${sample_id}.genealogy.report.html
    touch ${sample_id}.genealogy.audit.json
    """
}
