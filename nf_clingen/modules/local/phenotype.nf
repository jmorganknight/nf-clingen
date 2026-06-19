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

    # Extract genes for all space-separated HP terms
    phenotype_genes=""
    matched_terms=""
    for hp_term in ${phenotype_id}; do
        genes=\$(grep "^\${hp_term}\t" "${phenotype_map}" | cut -f3 || true)
        if [[ -n "\$genes" ]]; then
            phenotype_genes="\${phenotype_genes},\${genes}"
            matched_terms="\${matched_terms} \${hp_term}"
        else
            echo "WARNING: HP term \${hp_term} not found in gene map, skipping." >&2
        fi
    done
    phenotype_genes="\${phenotype_genes#,}"  # strip leading comma

    if [[ -z "\$phenotype_genes" ]]; then
        echo "WARNING: No HP terms matched in gene map for '${phenotype_id}'. Returning unfiltered variants." >&2
        cp "${annotation_tsv}" "${sample_id}.phenotype_filtered.tsv"
        exit 0
    fi
    echo "Matched phenotype terms:\${matched_terms}" >&2

    # Build lookup tables used for explicit HPO-to-variant correlation.
    # 1) phenotype_genes.list: all genes across requested HPO terms
    # 2) phenotype_gene_to_term.tsv: gene -> one or more matched HPO terms
    : > phenotype_gene_to_term.tsv
    for hp_term in ${phenotype_id}; do
        genes=\$(grep "^\${hp_term}\t" "${phenotype_map}" | cut -f3 || true)
        if [[ -n "\$genes" ]]; then
            echo "\$genes" | tr ',' '\n' | sed 's/^ *//;s/ *\$//' | awk -v term="\${hp_term}" 'NF>0{print toupper(\$0)"\t"term}' >> phenotype_gene_to_term.tsv
        fi
    done
    cut -f1 phenotype_gene_to_term.tsv | sort -u > phenotype_genes.list

    # Filter annotation TSV to variants whose genes intersect requested HPO terms.
    # Adds explicit correlation columns:
    # - matched_hpo_terms: HPO IDs implicated by matched genes
    # - matched_hpo_genes: matched gene symbols in this variant row
    {
        awk -F '\t' 'NR==1{print \$0"\tmatched_hpo_terms\tmatched_hpo_genes"; exit}' "${annotation_tsv}"
        awk -F '\t' '
            BEGIN {
                while ((getline line < "phenotype_gene_to_term.tsv") > 0) {
                    if (line == "") continue
                    split(line, p, "\t")
                    g = p[1]
                    t = p[2]
                    if (g != "" && t != "") {
                        if (!(g in gene_terms)) {
                            gene_terms[g] = t
                        } else if (index("," gene_terms[g] ",", "," t ",") == 0) {
                            gene_terms[g] = gene_terms[g] "," t
                        }
                    }
                }
                close("phenotype_gene_to_term.tsv")

                while ((getline line < "phenotype_genes.list") > 0) {
                    if (line != "") keep[line] = 1
                }
                close("phenotype_genes.list")
            }
            NR > 1 {
                genes_field = \$11  # gene column in clinical_annotation.tsv
                n = split(genes_field, items, /[|,;]/)
                matched = 0
                matched_genes = ""
                matched_terms = ""
                for (i = 1; i <= n; i++) {
                    g = items[i]
                    sub(/:.*/, "", g)  # convert GENE:ENTREZ to GENE
                    gsub(/ /, "", g)
                    g = toupper(g)
                    if (g != "" && (g in keep)) {
                        matched = 1
                        if (matched_genes == "") matched_genes = g
                        else if (index("," matched_genes ",", "," g ",") == 0) matched_genes = matched_genes "," g

                        terms = gene_terms[g]
                        m = split(terms, tarr, /,/) 
                        for (j = 1; j <= m; j++) {
                            t = tarr[j]
                            if (t == "") continue
                            if (matched_terms == "") matched_terms = t
                            else if (index("," matched_terms ",", "," t ",") == 0) matched_terms = matched_terms "," t
                        }
                    }
                }
                if (matched) print \$0 "\t" matched_terms "\t" matched_genes
            }
        ' "${annotation_tsv}"
    } > "${sample_id}.phenotype_filtered.tsv"

    lines=\$(tail -n +2 "${sample_id}.phenotype_filtered.tsv" | wc -l)
    echo "Phenotype filter (${phenotype_id}): retained \$lines variants" >&2
    """

    stub:
    """
    head -n 5 "${annotation_tsv}" > "${sample_id}.phenotype_filtered.tsv"
    """
}
