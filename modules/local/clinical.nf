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

    bcftools query \
      -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%FILTER[\t%GT]\t%INFO/CLNSIG\t%INFO/CLNDN\t%INFO/GENEINFO\t%INFO/ANN\n' \
      ${sample_id}.clinical.annotated.vcf.gz > ${sample_id}.clinical.raw.tsv || true

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
             if (clnsig ~ /PATHOGENIC/) { severity = "critical"; rank = 1 }
             else if (ann ~ /HIGH|STOP_GAINED|FRAMESHIFT|SPLICE/ || clnsig ~ /LIKELY_PATHOGENIC/) { severity = "high"; rank = 2 }
             else if (\$6 == "PASS" || \$6 == ".") { severity = "moderate"; rank = 3 }
             print rank,\$1,\$2,\$3,\$4,\$5,\$6,genotype,severity,(\$8 == "" ? "." : \$8),disease,gene,(\$11 == "" ? "." : \$11)
         }' ${sample_id}.clinical.raw.tsv > ${sample_id}.clinical.scored.tsv

    {
      head -n 1 ${sample_id}.clinical.scored.tsv | cut -f2-
      tail -n +2 ${sample_id}.clinical.scored.tsv | sort -t \$'\t' -k1,1n -k2,2 -k3,3n | head -n ${params.max_report_variants} | cut -f2-
    } > ${sample_id}.clinical_annotation.tsv

    total_scored=0
    critical_count=0
    high_count=0
    moderate_count=0
    if [[ -s ${sample_id}.clinical.scored.tsv ]]; then
      total_scored=\$(tail -n +2 ${sample_id}.clinical.scored.tsv | wc -l | tr -d ' ')
      critical_count=\$(awk -F '\t' 'NR>1 && \$9=="critical" {c++} END {print c+0}' ${sample_id}.clinical_annotation.tsv)
      high_count=\$(awk -F '\t' 'NR>1 && \$9=="high" {c++} END {print c+0}' ${sample_id}.clinical_annotation.tsv)
      moderate_count=\$(awk -F '\t' 'NR>1 && \$9=="moderate" {c++} END {print c+0}' ${sample_id}.clinical_annotation.tsv)
    fi

    cat > ${sample_id}.clinical_summary.json <<EOF
    {
      "sample_id": "${sample_id}",
      "source": "intrinsic_vcf_annotations",
      "prioritized_variants": ${total_scored},
      "critical": ${critical_count},
      "high": ${high_count},
      "moderate": ${moderate_count}
    }
    EOF
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
    bcftools annotate \
      -a ${clinvar_vcf} \
      -c ID,INFO/CLNSIG,INFO/CLNDN,INFO/GENEINFO \
      -Oz \
      -o ${sample_id}.clinical.annotated.vcf.gz \
      ${vcf}
    tabix -f -p vcf ${sample_id}.clinical.annotated.vcf.gz

    bcftools query \
      -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%FILTER[\t%GT]\t%INFO/CLNSIG\t%INFO/CLNDN\t%INFO/GENEINFO\t%INFO/ANN\n' \
      ${sample_id}.clinical.annotated.vcf.gz > ${sample_id}.clinical.raw.tsv || true

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
             if (clnsig ~ /PATHOGENIC/) { severity = "critical"; rank = 1 }
             else if (ann ~ /HIGH|STOP_GAINED|FRAMESHIFT|SPLICE/ || clnsig ~ /LIKELY_PATHOGENIC/) { severity = "high"; rank = 2 }
             else if (\$6 == "PASS" || \$6 == ".") { severity = "moderate"; rank = 3 }
             print rank,\$1,\$2,\$3,\$4,\$5,\$6,genotype,severity,(\$8 == "" ? "." : \$8),disease,gene,(\$11 == "" ? "." : \$11)
         }' ${sample_id}.clinical.raw.tsv > ${sample_id}.clinical.scored.tsv

    {
      head -n 1 ${sample_id}.clinical.scored.tsv | cut -f2-
      tail -n +2 ${sample_id}.clinical.scored.tsv | sort -t \$'\t' -k1,1n -k2,2 -k3,3n | head -n ${params.max_report_variants} | cut -f2-
    } > ${sample_id}.clinical_annotation.tsv

    total_scored=0
    critical_count=0
    high_count=0
    moderate_count=0
    if [[ -s ${sample_id}.clinical.scored.tsv ]]; then
      total_scored=\$(tail -n +2 ${sample_id}.clinical.scored.tsv | wc -l | tr -d ' ')
      critical_count=\$(awk -F '\t' 'NR>1 && \$9=="critical" {c++} END {print c+0}' ${sample_id}.clinical_annotation.tsv)
      high_count=\$(awk -F '\t' 'NR>1 && \$9=="high" {c++} END {print c+0}' ${sample_id}.clinical_annotation.tsv)
      moderate_count=\$(awk -F '\t' 'NR>1 && \$9=="moderate" {c++} END {print c+0}' ${sample_id}.clinical_annotation.tsv)
    fi

    cat > ${sample_id}.clinical_summary.json <<EOF
    {
      "sample_id": "${sample_id}",
      "source": "clinvar_overlay",
      "prioritized_variants": ${total_scored},
      "critical": ${critical_count},
      "high": ${high_count},
      "moderate": ${moderate_count}
    }
    EOF
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
    tuple val(sample_id), path(vcf), path(vcf_index), path(annotation_tsv), path(summary_json)

    output:
    tuple val(sample_id), path("${sample_id}.critical_variants.tsv"), path("${sample_id}.clinical_report.html"), path("${sample_id}.clinical_report.pdf"), emit: deliverables

    script:
    """
    set -euo pipefail

    python - <<'PY'
    import csv
    import html
    import json
    from pathlib import Path

    sample_id = "${sample_id}"
    annotation_path = Path("${annotation_tsv}")
    summary_path = Path("${summary_json}")
    html_path = Path(f"{sample_id}.clinical_report.html")
    critical_path = Path(f"{sample_id}.critical_variants.tsv")

    with summary_path.open() as handle:
        summary = json.load(handle)

    with annotation_path.open() as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        rows = list(reader)

    critical_rows = [row for row in rows if row['severity'] in {'critical', 'high'}]

    with critical_path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=reader.fieldnames, delimiter='\t')
        writer.writeheader()
        writer.writerows(critical_rows)

    def render_table(records):
        if not records:
            return '<tr><td colspan="8">No critical or high-priority variants passed local triage.</td></tr>'
        return '\n'.join(
            '<tr>'
            f"<td>{html.escape(record['chrom'])}</td>"
            f"<td>{html.escape(record['pos'])}</td>"
            f"<td>{html.escape(record['gene'])}</td>"
            f"<td>{html.escape(record['ref'])}</td>"
            f"<td>{html.escape(record['alt'])}</td>"
            f"<td>{html.escape(record['genotype'])}</td>"
            f"<td>{html.escape(record['severity'])}</td>"
            f"<td>{html.escape(record['clnsig'])}</td>"
            '</tr>'
            for record in records
        )

    variant_rows = render_table(critical_rows)
    method_rows = ''.join(
        [
            '<li>FASTQC and FASTP are applied unless the run uses --skip_qc.</li>',
            f'<li>Alignment engine: {html.escape("${params.aligner}")}; preprocess engine: {html.escape("${params.preprocess}")}.</li>',
            f'<li>Variant caller: {html.escape("${params.caller}")}.</li>',
            f'<li>Annotation stack: bcftools query/filter with source {html.escape(summary["source"])} and WeasyPrint PDF rendering.</li>'
        ]
    )

    html_path.write_text(f'''
    <!doctype html>
    <html lang=\"en\">
    <head>
      <meta charset=\"utf-8\">
      <title>nf-prism clinical report: {html.escape(sample_id)}</title>
      <style>
        @page {{ size: A4; margin: 16mm; }}
        body {{ font-family: DejaVu Sans, sans-serif; color: #1e293b; font-size: 10pt; line-height: 1.45; }}
        h1, h2, h3 {{ margin-bottom: 6px; }}
        .hero {{ background: linear-gradient(120deg, #e0f2fe, #fef3c7); border: 1px solid #cbd5e1; border-radius: 12px; padding: 18px; margin-bottom: 18px; }}
        .grid {{ display: grid; grid-template-columns: repeat(4, 1fr); gap: 8px; margin: 14px 0 18px; }}
        .card {{ background: #f8fafc; border: 1px solid #cbd5e1; border-radius: 10px; padding: 10px; }}
        .section {{ margin-bottom: 18px; }}
        table {{ width: 100%; border-collapse: collapse; }}
        th, td {{ border: 1px solid #cbd5e1; padding: 6px; vertical-align: top; }}
        th {{ background: #e2e8f0; text-align: left; }}
        .small {{ color: #475569; font-size: 8.5pt; }}
        ul {{ margin: 0; padding-left: 16px; }}
      </style>
    </head>
    <body>
      <div class=\"hero\">
        <h1>nf-prism clinical variant triage</h1>
        <div>Sample: <strong>{html.escape(sample_id)}</strong></div>
        <div>Annotated VCF: {html.escape(Path("${vcf}").name)}</div>
        <div class=\"small\">Local automated triage for review workflows. This is not a substitute for molecular pathology sign-out or regulated ACMG adjudication.</div>
      </div>

      <div class=\"grid\">
        <div class=\"card\"><strong>Prioritized</strong><br>{summary['prioritized_variants']}</div>
        <div class=\"card\"><strong>Critical</strong><br>{summary['critical']}</div>
        <div class=\"card\"><strong>High</strong><br>{summary['high']}</div>
        <div class=\"card\"><strong>Moderate</strong><br>{summary['moderate']}</div>
      </div>

      <div class=\"section\">
        <h2>Executive summary</h2>
        <p>The callset was routed through the nf-prism clinical branch and prioritized with a local open-source annotation stack. Variants carrying pathogenic ClinVar evidence or high-impact functional consequences are promoted into the review table below.</p>
      </div>

      <div class=\"section\">
        <h2>Flagged variants</h2>
        <table>
          <thead>
            <tr><th>Chr</th><th>Pos</th><th>Gene</th><th>Ref</th><th>Alt</th><th>GT</th><th>Priority</th><th>ClinVar</th></tr>
          </thead>
          <tbody>
            {variant_rows}
          </tbody>
        </table>
      </div>

      <div class=\"section\">
        <h2>Interpretation notes</h2>
        <p>Priority tiers are derived from ClinVar significance when available, followed by high-impact consequence terms present in the source VCF annotations. Moderate findings retain passing variants with non-reference genotypes for downstream review but are not included in the critical deliverable TSV unless promoted.</p>
      </div>

      <div class=\"section\">
        <h2>Methods</h2>
        <ul>{method_rows}</ul>
      </div>

      <div class=\"section\">
        <h2>Limitations</h2>
        <p class=\"small\">Population frequency filters, phenotype-specific gene panels, transcript curation, inheritance models, CNV/SV interpretation, and final ACMG criteria scoring remain institution-specific responsibilities outside this boilerplate workflow.</p>
      </div>
    </body>
    </html>
    '''.strip())
    PY

    weasyprint ${sample_id}.clinical_report.html ${sample_id}.clinical_report.pdf
    """

    stub:
    """
    touch ${sample_id}.critical_variants.tsv ${sample_id}.clinical_report.html ${sample_id}.clinical_report.pdf
    """
}