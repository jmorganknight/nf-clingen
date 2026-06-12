#!/usr/bin/env python3

import argparse
import csv
import html
import json
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build clinical HTML/PDF report inputs")
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--annotation-tsv", required=True)
    parser.add_argument("--summary-json", required=True)
    parser.add_argument("--annotated-vcf", required=True)
    parser.add_argument("--aligner", required=True)
    parser.add_argument("--preprocess", required=True)
    parser.add_argument("--caller", required=True)
    parser.add_argument("--reference", required=True)
    parser.add_argument("--clinvar-resource", default="disabled")
    parser.add_argument("--gnomad-resource", default="disabled")
    parser.add_argument("--pipeline-version", default="unknown")
    parser.add_argument("--nextflow-version", default="unknown")
    parser.add_argument("--aligner-version", default="unknown")
    parser.add_argument("--preprocess-version", default="unknown")
    parser.add_argument("--caller-version", default="unknown")
    parser.add_argument("--annotation-version", default="unknown")
    parser.add_argument("--reporting-version", default="unknown")
    parser.add_argument("--patient-phenotype", default="none")
    parser.add_argument("--run-id", default="unknown")
    parser.add_argument("--run-timestamp", default="unknown")
    parser.add_argument("--run-parameters", default="")
    return parser.parse_args()


def render_table(records: list[dict[str, str]]) -> str:
    if not records:
        return '<tr><td colspan="8">No critical or high-priority variants passed local triage.</td></tr>'

    return "\n".join(
        "<tr>"
        f"<td>{html.escape(record['chrom'])}</td>"
        f"<td>{html.escape(record['pos'])}</td>"
        f"<td>{html.escape(record['gene'])}</td>"
        f"<td>{html.escape(record['ref'])}</td>"
        f"<td>{html.escape(record['alt'])}</td>"
        f"<td>{html.escape(record['genotype'])}</td>"
        f"<td>{html.escape(record['severity'])}</td>"
        f"<td>{html.escape(record['clnsig'])}</td>"
        "</tr>"
        for record in records
    )


def main() -> None:
    args = parse_args()

    sample_id = args.sample_id
    annotation_path = Path(args.annotation_tsv)
    # Sanitize sample_id for filesystem use (strip characters unsafe in filenames)
    import re as _re
    safe_sample_id = _re.sub(r"[^\w\-.]", "_", sample_id)
    summary_path = Path(args.summary_json)
    html_path = Path(f"{safe_sample_id}.clinical_report.html")
    critical_path = Path(f"{safe_sample_id}.critical_variants.tsv")

    with summary_path.open() as handle:
        summary = json.load(handle)

    with annotation_path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)

    critical_rows = [row for row in rows if row.get("severity") in {"critical", "high"}]

    with critical_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=reader.fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(critical_rows)

    variant_rows = render_table(critical_rows)
    clinvar_resource = args.clinvar_resource
    annotation_db_rows = [
      f"<li>Primary annotation source: {html.escape(str(summary.get('source', '.')))}</li>",
      (
        f"<li>ClinVar resource: {html.escape(Path(clinvar_resource).name)}"
        if clinvar_resource not in {"", "disabled", "null", "none"}
        else "<li>ClinVar resource: not provided for this run</li>"
      ),
      f"<li>Reference genome: {html.escape(Path(args.reference).name)}</li>",
    ]

    transparency_rows = [
      ("Pipeline", f"nf-clingen {html.escape(args.pipeline_version)}"),
      ("Workflow engine", f"Nextflow {html.escape(args.nextflow_version)}"),
      ("Aligner", f"{html.escape(args.aligner)} ({html.escape(args.aligner_version)})"),
      ("Preprocess", f"{html.escape(args.preprocess)} ({html.escape(args.preprocess_version)})"),
      ("Variant caller", f"{html.escape(args.caller)} ({html.escape(args.caller_version)})"),
      ("Clinical annotation", html.escape(args.annotation_version)),
      ("Report rendering", html.escape(args.reporting_version)),
      ("Patient phenotype input", html.escape(args.patient_phenotype)),
    ]

    transparency_table = "\n".join(
      f"<tr><th>{label}</th><td>{value}</td></tr>" for label, value in transparency_rows
    )

    # Build run execution audit subsection
    run_audit_rows = [
      ("Run ID", html.escape(args.run_id)),
      ("Run timestamp", html.escape(args.run_timestamp)),
    ]
    run_audit_table = "\n".join(
      f"<tr><th>{label}</th><td>{value}</td></tr>" for label, value in run_audit_rows
    )
    
    # Parse and format run parameters
    run_params_html = "<li>No parameters captured</li>"
    if args.run_parameters and args.run_parameters.strip():
      params_list = []
      for line in args.run_parameters.split("\n"):
        if line.strip():
          # Split by '=' if it looks like a key=value pair
          if "=" in line:
            key, value = line.split("=", 1)
            params_list.append(f"<li><strong>{html.escape(key)}:</strong> {html.escape(value)}</li>")
          else:
            params_list.append(f"<li>{html.escape(line)}</li>")
      run_params_html = "\n".join(params_list) if params_list else "<li>No parameters captured</li>"

    method_rows = "".join(
        [
            "<li>FASTQC and FASTP are applied unless the run uses --skip_qc.</li>",
            f"<li>Alignment engine: {html.escape(args.aligner)}; preprocess engine: {html.escape(args.preprocess)}.</li>",
            f"<li>Variant caller: {html.escape(args.caller)}.</li>",
            f"<li>Annotation stack: bcftools query/filter with source {html.escape(str(summary.get('source', '.')))} and WeasyPrint PDF rendering.</li>",
        ]
    )

    html_path.write_text(
        f"""
<!doctype html>
<html lang=\"en\">
<head>
  <meta charset=\"utf-8\">
  <title>nf-clingen clinical report: {html.escape(sample_id)}</title>
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
    <h1>nf-clingen clinical variant triage</h1>
    <div>Sample: <strong>{html.escape(sample_id)}</strong></div>
    <div>Annotated VCF: {html.escape(Path(args.annotated_vcf).name)}</div>
    <div class=\"small\">Local automated triage for review workflows. This is not a substitute for molecular pathology sign-out or regulated ACMG adjudication.</div>
  </div>

  <div class=\"grid\">
    <div class=\"card\"><strong>Prioritized</strong><br>{summary.get('prioritized_variants', 0)}</div>
    <div class=\"card\"><strong>Critical</strong><br>{summary.get('critical', 0)}</div>
    <div class=\"card\"><strong>High</strong><br>{summary.get('high', 0)}</div>
    <div class=\"card\"><strong>Moderate</strong><br>{summary.get('moderate', 0)}</div>
  </div>

  <div class=\"section\">
    <h2>Executive summary</h2>
    <p>The callset was routed through the nf-clingen clinical branch and prioritized with a local open-source annotation stack. Variants carrying pathogenic ClinVar evidence or high-impact functional consequences are promoted into the review table below.</p>
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
    <h2>Provenance and transparency</h2>
    
    <h3>Run execution audit</h3>
    <table>
      <tbody>
        {run_audit_table}
      </tbody>
    </table>
    <p class=\"small\"><strong>Run parameters:</strong></p>
    <ul>
      {run_params_html}
    </ul>

    <h3>Software and tool versions</h3>
    <table>
      <tbody>
        {transparency_table}
      </tbody>
    </table>
    
    <h3>Reference and annotation databases</h3>
    <ul>
      {''.join(annotation_db_rows)}
    </ul>
    <p class=\"small\">This section provides run-level software and data provenance to support auditability, reproducibility, and regulatory/quality review workflows.</p>
  </div>

  <div class=\"section\">
    <h2>Limitations</h2>
    <p class=\"small\">Population frequency filters, phenotype-specific gene panels, transcript curation, inheritance models, CNV/SV interpretation, and final ACMG criteria scoring remain institution-specific responsibilities outside this boilerplate workflow.</p>
  </div>
</body>
</html>
""".strip()
    )


if __name__ == "__main__":
    main()
