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

    # Extract chromosomes present in this VCF (use grep+sed; awk capture groups are gawk-only)
    chroms=\$(bcftools view -h ${vcf} | grep '^##contig' | sed -n 's/.*ID=\\([^,>]*\\).*/\\1/p' | grep -E '^(chr)?[0-9XY]+\$' | sort -V | tr '\\n' ' ')

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

    # Extract chromosome list from VCF header
    chroms=\$(zcat ${phased_vcf} | grep '^##contig' | sed -n 's/.*ID=\\([^,>]*\\).*/\\1/p' | grep -E '^(chr)?[0-9XY]+\$' | sort -V)

    if [[ -z "\${chroms}" ]]; then
        echo "WARNING: no ##contig headers; scanning records for chromosome names" >&2
        chroms=\$(zcat ${phased_vcf} | grep -v '^#' | cut -f1 | sort -uV | grep -E '^(chr)?[0-9XY]+\$')
    fi

    for chrom in \${chroms}; do
        # Locate reference panel — prefer hg38 VCF, then hg38 bref3, then any build
        ref=\$(find ${ref_panel_dir} -name "*hg38*.\${chrom}.vcf.gz" 2>/dev/null | sort | head -1)
        if [[ -z "\${ref}" ]]; then ref=\$(find ${ref_panel_dir} -name "*hg38*.\${chrom}.bref3" 2>/dev/null | sort | head -1); fi
        if [[ -z "\${ref}" ]]; then ref=\$(find ${ref_panel_dir} -name "*.\${chrom}.vcf.gz" 2>/dev/null | sort | head -1); fi
        if [[ -z "\${ref}" ]]; then ref=\$(find ${ref_panel_dir} -name "*.\${chrom}.bref3" 2>/dev/null | sort | head -1); fi

        # Locate genetic map
        map_arg=""
        map_file=\$(find ${beagle_map_dir} -name "*.\${chrom}.map.gz" 2>/dev/null | sort | head -1)
        if [[ -n "\${map_file}" ]]; then map_arg="map=\${map_file}"; fi

        if [[ -n "\${ref}" ]]; then
            java -Xmx\${JAVA_HEAP:-16g} -jar ${beagle_jar} \\
                gt=${phased_vcf} \\
                ref=\${ref} \\
                chrom=\${chrom} \\
                \${map_arg} \\
                out=${sample_id}.\${chrom}.imputed \\
                nthreads=${task.cpus} \\
                2>&1 | tee -a beagle_\${chrom}.log
        else
            echo "WARNING: no ref panel for \${chrom}; running phasing refinement only" >&2
            java -Xmx\${JAVA_HEAP:-16g} -jar ${beagle_jar} \\
                gt=${phased_vcf} \\
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

    # Index each per-chromosome VCF then concat in sorted order
    for f in \$(ls *.imputed.vcf.gz | sort -V); do
        bcftools index --tbi "\${f}"
    done

    bcftools concat --allow-overlaps --rm-dups all \\
        \$(ls *.imputed.vcf.gz | sort -V | tr '\\n' ' ') \\
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

    script:
    """
    set -euo pipefail

    phased_variants=\$(bcftools stats ${phased_vcf} | awk '/^SN.*number of SNPs/ {print \$NF}')
    imputed_variants=\$(bcftools stats ${imputed_vcf} | awk '/^SN.*number of SNPs/ {print \$NF}')

    cat > ${sample_id}.genealogy.manifest.txt <<EOF
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
run_date=\$(date -u +%Y-%m-%dT%H:%M:%SZ)
EOF
    """

    stub:
    """
    touch ${sample_id}.genealogy.manifest.txt
    """
}