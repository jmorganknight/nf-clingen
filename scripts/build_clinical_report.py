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
    summary_path = Path(args.summary_json)
    html_path = Path(f"{sample_id}.clinical_report.html")
    critical_path = Path(f"{sample_id}.critical_variants.tsv")

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
    <h2>Limitations</h2>
    <p class=\"small\">Population frequency filters, phenotype-specific gene panels, transcript curation, inheritance models, CNV/SV interpretation, and final ACMG criteria scoring remain institution-specific responsibilities outside this boilerplate workflow.</p>
  </div>
</body>
</html>
""".strip()
    )


if __name__ == "__main__":
    main()
