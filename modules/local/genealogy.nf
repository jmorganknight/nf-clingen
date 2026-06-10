// Phase variants with Eagle2 then impute with Beagle 5.4.
// Required inputs:
//   params.eagle2_genetic_map  — path to Eagle2 genetic map directory (e.g. Eagle_v2.4.1/tables/)
//   params.beagle_ref_panel    — path to per-chromosome Beagle reference panel VCFs, named *.chr{CHROM}.vcf.gz
//   params.beagle_genetic_map  — path to per-chromosome Beagle genetic maps, named *.chr{CHROM}.map.gz (optional; improves accuracy)
//   params.beagle_jar          — path to beagle.*.jar (default: /opt/beagle/beagle.jar inside container)
// Container images used:
//   phase  : quay.io/biocontainers/eagle:2.4.1--h9ee0642_4
//   impute : quay.io/biocontainers/beagle:5.4_22Jul22.46e--hdfd78af_0
process PHASE_EAGLE2 {
    tag "${sample_id}"
    container 'quay.io/biocontainers/eagle:2.4.1--h9ee0642_4'
    publishDir "${params.outdir}/genealogy/phased", mode: 'copy'

    input:
    tuple val(sample_id), path(vcf), path(vcf_index)
    path(genetic_map_dir)

    output:
    tuple val(sample_id), path("${sample_id}.phased.vcf.gz"), path("${sample_id}.phased.vcf.gz.tbi"), emit: phased

    script:
    """
    set -euo pipefail

    # Extract chromosomes present in this VCF
    chroms=\$(bcftools view -h ${vcf} | awk '/^#CHROM/ {next} /^##contig/ {
        match(\$0, /ID=([^,>]+)/, a); print a[1]
    }' | grep -E '^(chr)?[0-9XY]+\$' | sort -V | tr '\\n' ' ')

    # Phase each chromosome separately then merge
    phased_vcfs=""
    for chrom in \$chroms; do
        eagle \\
            --vcfTarget=${vcf} \\
            --geneticMapFile=${genetic_map_dir}/genetic_map_hg38_withX.txt.gz \\
            --chrom=\${chrom} \\
            --outPrefix=${sample_id}.\${chrom}.eagle2 \\
            --numThreads=${task.cpus} \\
            2>&1 | tee -a eagle2_\${chrom}.log
        phased_vcfs="\${phased_vcfs} ${sample_id}.\${chrom}.eagle2.vcf.gz"
    done

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
    publishDir "${params.outdir}/genealogy/imputed", mode: 'copy'

    input:
    tuple val(sample_id), path(phased_vcf), path(phased_vcf_index)
    path(ref_panel_dir)
    path(beagle_map_dir)

    output:
    tuple val(sample_id), path("${sample_id}.imputed.vcf.gz"), path("${sample_id}.imputed.vcf.gz.tbi"), emit: imputed

    script:
    // Beagle jar location: prefer params, fall back to container default
    def beagle_jar = params.beagle_jar ?: '/opt/beagle/beagle.jar'
    """
    set -euo pipefail

    chroms=\$(bcftools view -h ${phased_vcf} | awk '/^##contig/ {
        match(\$0, /ID=([^,>]+)/, a); print a[1]
    }' | grep -E '^(chr)?[0-9XY]+\$' | sort -V)

    imputed_vcfs=""
    for chrom in \$chroms; do
        # Locate reference panel for this chromosome (supports chr-prefixed and bare names)
        ref=\$(ls ${ref_panel_dir}/*.chr\${chrom}.vcf.gz ${ref_panel_dir}/*.\${chrom}.vcf.gz 2>/dev/null | head -1 || true)
        if [[ -z "\$ref" ]]; then
            echo "WARNING: no reference panel found for \${chrom}, skipping imputation for this chromosome." >&2
            bcftools view ${phased_vcf} \${chrom} -Oz -o ${sample_id}.\${chrom}.imputed.vcf.gz
            bcftools index --tbi ${sample_id}.\${chrom}.imputed.vcf.gz
            imputed_vcfs="\${imputed_vcfs} ${sample_id}.\${chrom}.imputed.vcf.gz"
            continue
        fi

        # Optionally use genetic map for this chromosome
        map_arg=""
        map_file=\$(ls ${beagle_map_dir}/*.chr\${chrom}.map.gz ${beagle_map_dir}/*.\${chrom}.map.gz 2>/dev/null | head -1 || true)
        if [[ -n "\$map_file" ]]; then map_arg="map=\$map_file"; fi

        java -Xmx\${JAVA_HEAP:-16g} -jar ${beagle_jar} \\
            gt=${phased_vcf} \\
            ref=\$ref \\
            chrom=\${chrom} \\
            \${map_arg} \\
            out=${sample_id}.\${chrom}.imputed \\
            nthreads=${task.cpus} \\
            2>&1 | tee -a beagle_\${chrom}.log

        bcftools index --tbi ${sample_id}.\${chrom}.imputed.vcf.gz
        imputed_vcfs="\${imputed_vcfs} ${sample_id}.\${chrom}.imputed.vcf.gz"
    done

    bcftools concat --allow-overlaps --rm-dups all \\
        \${imputed_vcfs} \\
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