"""Tests for scripts/build_clinical_report.py"""

import csv
import json
import sys
from pathlib import Path

import pytest

# Make scripts importable
sys.path.insert(0, str(Path(__file__).parent.parent / "scripts"))

import build_clinical_report as bcr


# ── fixtures ────────────────────────────────────────────────────────────────

@pytest.fixture()
def summary_json(tmp_path):
    data = {
        "prioritized_variants": 5,
        "critical": 2,
        "high": 3,
        "moderate": 0,
        "source": "bcftools-triage",
    }
    p = tmp_path / "summary.json"
    p.write_text(json.dumps(data))
    return p


@pytest.fixture()
def annotation_tsv(tmp_path):
    rows = [
        {"chrom": "chr1", "pos": "100", "gene": "BRCA1", "ref": "A", "alt": "T",
         "genotype": "0/1", "severity": "critical", "clnsig": "Pathogenic"},
        {"chrom": "chr2", "pos": "200", "gene": "TP53",  "ref": "G", "alt": "C",
         "genotype": "1/1", "severity": "high",     "clnsig": "Likely_pathogenic"},
        {"chrom": "chr3", "pos": "300", "gene": "MYH7",  "ref": "C", "alt": "A",
         "genotype": "0/1", "severity": "moderate", "clnsig": "VUS"},
    ]
    p = tmp_path / "annotations.tsv"
    with p.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    return p


@pytest.fixture()
def base_args(tmp_path, summary_json, annotation_tsv, monkeypatch):
    monkeypatch.chdir(tmp_path)
    return [
        "--sample-id", "SAMPLE001",
        "--annotation-tsv", str(annotation_tsv),
        "--summary-json", str(summary_json),
        "--annotated-vcf", "sample.vcf.gz",
        "--aligner", "minimap2",
        "--preprocess", "samtools",
        "--caller", "deepvariant",
        "--reference", "GRCh38.fasta",
    ]


# ── render_table ────────────────────────────────────────────────────────────

class TestRenderTable:
    def test_empty_returns_no_variants_message(self):
        result = bcr.render_table([])
        assert "No critical" in result

    def test_single_row_renders_all_fields(self):
        row = {"chrom": "chr1", "pos": "100", "gene": "BRCA1",
               "ref": "A", "alt": "T", "genotype": "0/1",
               "severity": "critical", "clnsig": "Pathogenic"}
        result = bcr.render_table([row])
        assert "BRCA1" in result
        assert "chr1" in result
        assert "Pathogenic" in result

    def test_xss_is_escaped(self):
        row = {"chrom": "<script>alert(1)</script>", "pos": "1", "gene": "G",
               "ref": "A", "alt": "T", "genotype": "0/1",
               "severity": "critical", "clnsig": "."}
        result = bcr.render_table([row])
        assert "<script>" not in result
        assert "&lt;script&gt;" in result

    def test_multiple_rows_all_rendered(self):
        rows = [
            {"chrom": f"chr{i}", "pos": str(i), "gene": f"G{i}",
             "ref": "A", "alt": "T", "genotype": "0/1",
             "severity": "critical", "clnsig": "."}
            for i in range(1, 6)
        ]
        result = bcr.render_table(rows)
        for i in range(1, 6):
            assert f"chr{i}" in result


# ── main() integration ───────────────────────────────────────────────────────

class TestMain:
    def test_creates_html_output(self, tmp_path, base_args, monkeypatch):
        monkeypatch.chdir(tmp_path)
        sys.argv = ["build_clinical_report.py"] + base_args
        bcr.main()
        assert (tmp_path / "SAMPLE001.clinical_report.html").exists()

    def test_creates_critical_tsv(self, tmp_path, base_args, monkeypatch):
        monkeypatch.chdir(tmp_path)
        sys.argv = ["build_clinical_report.py"] + base_args
        bcr.main()
        tsv = tmp_path / "SAMPLE001.critical_variants.tsv"
        assert tsv.exists()
        with tsv.open() as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            rows = list(reader)
        # Only critical and high rows should appear
        assert all(r["severity"] in {"critical", "high"} for r in rows)
        assert len(rows) == 2

    def test_html_contains_sample_id(self, tmp_path, base_args, monkeypatch):
        monkeypatch.chdir(tmp_path)
        sys.argv = ["build_clinical_report.py"] + base_args
        bcr.main()
        html = (tmp_path / "SAMPLE001.clinical_report.html").read_text()
        assert "SAMPLE001" in html

    def test_html_contains_variant_gene(self, tmp_path, base_args, monkeypatch):
        monkeypatch.chdir(tmp_path)
        sys.argv = ["build_clinical_report.py"] + base_args
        bcr.main()
        html = (tmp_path / "SAMPLE001.clinical_report.html").read_text()
        assert "BRCA1" in html

    def test_clinvar_disabled_label(self, tmp_path, base_args, monkeypatch):
        monkeypatch.chdir(tmp_path)
        sys.argv = ["build_clinical_report.py"] + base_args + ["--clinvar-resource", "disabled"]
        bcr.main()
        html = (tmp_path / "SAMPLE001.clinical_report.html").read_text()
        assert "not provided" in html

    def test_clinvar_path_shown(self, tmp_path, base_args, monkeypatch):
        monkeypatch.chdir(tmp_path)
        sys.argv = ["build_clinical_report.py"] + base_args + [
            "--clinvar-resource", "/data/clinvar/clinvar_GRCh38.vcf.gz"
        ]
        bcr.main()
        html = (tmp_path / "SAMPLE001.clinical_report.html").read_text()
        assert "clinvar_GRCh38.vcf.gz" in html

    def test_run_parameters_rendered(self, tmp_path, base_args, monkeypatch):
        monkeypatch.chdir(tmp_path)
        sys.argv = ["build_clinical_report.py"] + base_args + [
            "--run-parameters", "caller=deepvariant\nmax_cpus=16"
        ]
        bcr.main()
        html = (tmp_path / "SAMPLE001.clinical_report.html").read_text()
        assert "max_cpus" in html

    def test_xss_in_sample_id_escaped(self, tmp_path, summary_json, annotation_tsv, monkeypatch):
        """sample_id with special chars: filename is sanitized, HTML content is escaped."""
        monkeypatch.chdir(tmp_path)
        sys.argv = [
            "build_clinical_report.py",
            "--sample-id", "<script>evil</script>",
            "--annotation-tsv", str(annotation_tsv),
            "--summary-json", str(summary_json),
            "--annotated-vcf", "s.vcf.gz",
            "--aligner", "minimap2",
            "--preprocess", "samtools",
            "--caller", "deepvariant",
            "--reference", "ref.fa",
        ]
        bcr.main()
        # Filename must be sanitized (no raw < > chars in filename)
        html_files = list(tmp_path.glob("*.clinical_report.html"))
        assert len(html_files) == 1
        fname = html_files[0].name
        assert "<" not in fname and ">" not in fname
        # HTML content must escape the tag
        content = html_files[0].read_text()
        assert "<script>evil</script>" not in content
        assert "&lt;script&gt;" in content
