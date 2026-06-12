"""Tests for the genealogy build-compatibility guard in main.nf.

We test the Groovy helper logic by re-implementing the same inference rules
in Python and verifying they match the documented behaviour. This gives fast
deterministic coverage without running Nextflow.
"""

import re


# ── pure-Python reimplementation of inferBuildTag / inferPanelBuildTag ───────

def infer_build_tag(value: str) -> str:
    """Mirror of inferBuildTag() in main.nf"""
    v = value.lower()
    if re.search(r"grch38|hg38|b38", v):
        return "grch38"
    if re.search(r"grch37|hg37|b37", v):
        return "grch37"
    return "unknown"


def infer_panel_build_tag(files: list[str]) -> str:
    """Mirror of inferPanelBuildTag() in main.nf — input is list of filenames."""
    tags = {infer_build_tag(f) for f in files}
    if "grch38" in tags:
        return "grch38"
    if "grch37" in tags:
        return "grch37"
    return "unknown"


# ── inferBuildTag ─────────────────────────────────────────────────────────────

class TestInferBuildTag:
    def test_grch38_exact(self):
        assert infer_build_tag("GRCh38") == "grch38"

    def test_hg38_variant(self):
        assert infer_build_tag("hg38") == "grch38"

    def test_b38_variant(self):
        assert infer_build_tag("b38") == "grch38"

    def test_grch37_exact(self):
        assert infer_build_tag("GRCh37") == "grch37"

    def test_hg37_variant(self):
        assert infer_build_tag("hg37") == "grch37"

    def test_b37_variant(self):
        assert infer_build_tag("b37") == "grch37"

    def test_embedded_in_path(self):
        assert infer_build_tag("/data/ref/GRCh38.fasta") == "grch38"

    def test_embedded_in_panel_filename(self):
        assert infer_build_tag("1000G_phase3.hg38.chr1.vcf.gz") == "grch38"

    def test_unknown_returns_unknown(self):
        assert infer_build_tag("mm10") == "unknown"

    def test_empty_string(self):
        assert infer_build_tag("") == "unknown"


# ── inferPanelBuildTag ────────────────────────────────────────────────────────

class TestInferPanelBuildTag:
    def test_all_hg38_files(self):
        files = [
            "1000G_phase3.hg38.chr1.vcf.gz",
            "1000G_phase3.hg38.chr2.vcf.gz",
        ]
        assert infer_panel_build_tag(files) == "grch38"

    def test_all_b37_files(self):
        files = [
            "1000G_phase3.b37.chr1.vcf.gz",
            "1000G_phase3.b37.chr2.vcf.gz",
        ]
        assert infer_panel_build_tag(files) == "grch37"

    def test_mixed_prefers_grch38(self):
        """When both builds present, hg38 wins (setup_genealogy_resources.sh
        generates both; guard should accept the run as hg38)."""
        files = [
            "1000G_phase3.hg38.chr1.vcf.gz",
            "1000G_phase3.b37.chr1.vcf.gz",
        ]
        assert infer_panel_build_tag(files) == "grch38"

    def test_empty_list_returns_unknown(self):
        assert infer_panel_build_tag([]) == "unknown"

    def test_unknown_filenames_return_unknown(self):
        assert infer_panel_build_tag(["panel.chr1.vcf.gz"]) == "unknown"


# ── build compatibility guard logic ──────────────────────────────────────────

class TestBuildCompatibilityGuard:
    """Simulate the guard block in main.nf that blocks mismatched runs."""

    @staticmethod
    def guard(reference: str, map_dir: str, panel_files: list[str]) -> str | None:
        """Returns None if compatible, error string if mismatch."""
        ref_tag = infer_build_tag(reference)
        map_tag = infer_build_tag(map_dir)
        panel_tag = infer_panel_build_tag(panel_files)
        tags = {t for t in [ref_tag, map_tag, panel_tag] if t != "unknown"}
        if len(tags) > 1:
            return f"Build mismatch: ref={ref_tag}, map={map_tag}, panel={panel_tag}"
        return None

    def test_all_grch38_passes(self):
        assert self.guard(
            reference="data/GRCh38.fasta",
            map_dir="data/eagle2/tables_hg38",
            panel_files=["1000G_phase3.hg38.chr1.vcf.gz"],
        ) is None

    def test_ref_grch38_panel_b37_blocked(self):
        result = self.guard(
            reference="data/GRCh38.fasta",
            map_dir="data/eagle2/tables_hg38",
            panel_files=["1000G_phase3.b37.chr1.vcf.gz"],
        )
        assert result is not None
        assert "mismatch" in result.lower()

    def test_ref_grch37_map_grch38_blocked(self):
        result = self.guard(
            reference="data/GRCh37.fasta",
            map_dir="data/eagle2/tables_hg38",
            panel_files=["1000G_phase3.b37.chr1.vcf.gz"],
        )
        assert result is not None

    def test_all_unknown_passes(self):
        """If no build can be inferred from any name, guard is permissive."""
        assert self.guard(
            reference="data/ref.fa",
            map_dir="data/maps",
            panel_files=["panel.chr1.vcf.gz"],
        ) is None

    def test_mixed_hg38_and_b37_panel_passes(self):
        """Both hg38 + b37 panels in dir → inferred as grch38 → matches ref."""
        assert self.guard(
            reference="data/GRCh38.fasta",
            map_dir="data/eagle2/tables_hg38",
            panel_files=[
                "1000G_phase3.hg38.chr1.vcf.gz",
                "1000G_phase3.b37.chr1.vcf.gz",  # legacy file still present
            ],
        ) is None
