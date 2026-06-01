process FILTER_BY_PHENOTYPE {
    tag "${sample_id}"
    publishDir "${params.outdir}/clinical/phenotype_filtered", mode: 'copy'

    input:
    tuple val(sample_id), path(annotation_tsv), val(phenotype_id), path(phenotype_map)

    output:
    tuple val(sample_id), path("${sample_id}.phenotype_filtered.tsv"), emit: filtered

    script:
    """
    set -euo pipefail

    # Extract genes relevant to patient phenotype
    phenotype_genes=\$(grep "^${phenotype_id}\t" "${phenotype_map}" | cut -f3)
    
    if [[ -z "\$phenotype_genes" ]]; then
        echo "WARNING: Phenotype ${phenotype_id} not found in gene map. Returning unfiltered variants." >&2
        cp "${annotation_tsv}" "${sample_id}.phenotype_filtered.tsv"
        exit 0
    fi

    # Convert comma-separated genes to regex for filtering
    gene_regex=\$(echo "\$phenotype_genes" | sed 's/,/|/g')

    # Filter annotation TSV to only genes relevant to patient phenotype
    # Keep header, then filter rows where gene field matches phenotype-relevant genes
    {
        head -n 1 "${annotation_tsv}"
        awk -F '\t' -v regex="\$gene_regex" '
            NR > 1 {
                gene = \$11  # gene is field 11 in clinical_annotation.tsv
                # Extract gene name (may contain multiple genes, separated by /)
                genes = \$11
                if (genes ~ regex) {
                    print \$0
                }
            }' "${annotation_tsv}"
    } > "${sample_id}.phenotype_filtered.tsv"

    lines=\$(tail -n +2 "${sample_id}.phenotype_filtered.tsv" | wc -l)
    echo "Phenotype filter (${phenotype_id}): retained \$lines variants" >&2
    """

    stub:
    """
    head -n 5 "${annotation_tsv}" > "${sample_id}.phenotype_filtered.tsv"
    """
}
