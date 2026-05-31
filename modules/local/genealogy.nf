process RUN_GENEALOGY_STACK {
    tag "${sample_id}"
    publishDir "${params.outdir}/genealogy", mode: 'copy'

    input:
    tuple val(sample_id), path(vcf), path(vcf_index)

    output:
    tuple val(sample_id), path("${sample_id}.phased.vcf.gz"), path("${sample_id}.imputed.vcf.gz"), path("${sample_id}.genealogy.manifest.txt"), emit: deliverables

    script:
    """
    set -euo pipefail
    ln -sf ${vcf} ${sample_id}.phased.vcf.gz
    ln -sf ${vcf_index} ${sample_id}.phased.vcf.gz.tbi
    ln -sf ${vcf} ${sample_id}.imputed.vcf.gz
    ln -sf ${vcf_index} ${sample_id}.imputed.vcf.gz.tbi

    cat > ${sample_id}.genealogy.manifest.txt <<EOF
    sample_id=${sample_id}
    source_vcf=${vcf}
    phase_engine=Eagle2
    imputation_engine=Beagle 5.4
    eagle2_stub=eagle --vcfTarget ${vcf} --geneticMapFile genetic_map.txt.gz --outPrefix ${sample_id}.eagle2
    beagle_stub=java -jar beagle.29Oct24.c8e.jar gt=${sample_id}.eagle2.vcf.gz out=${sample_id}.imputed ref=reference_panel.vcf.gz map=genetic_map.txt.gz
    note=Routing stub only. Attach population-specific maps, reference panels, and downstream kinship ancestry analytics before deployment.
    EOF
    """

    stub:
    """
    touch ${sample_id}.phased.vcf.gz ${sample_id}.imputed.vcf.gz ${sample_id}.genealogy.manifest.txt
    """
}