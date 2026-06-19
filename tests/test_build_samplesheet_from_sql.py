"""Tests for scripts/build_samplesheet_from_sql.py"""

import csv
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent / "scripts"))

import build_samplesheet_from_sql as bss


# ── helpers ──────────────────────────────────────────────────────────────────

def make_sqlite_db(tmp_path, rows):
    """Create a tiny SQLite DB with a fastq_manifest table."""
    import sqlite3
    db = tmp_path / "test.db"
    conn = sqlite3.connect(db)
    conn.execute(
        "CREATE TABLE fastq_manifest (sample_id TEXT, r1 TEXT, r2 TEXT)"
    )
    conn.executemany("INSERT INTO fastq_manifest VALUES (?, ?, ?)", rows)
    conn.commit()
    conn.close()
    return db


# ── load_query ───────────────────────────────────────────────────────────────

class TestLoadQuery:
    def test_inline_query_returned(self, tmp_path):
        args = bss.parse_args.__wrapped__ if hasattr(bss.parse_args, "__wrapped__") else None
        # test via namespace directly
        import argparse
        ns = argparse.Namespace(query="SELECT 1", query_file=None)
        assert bss.load_query(ns) == "SELECT 1"

    def test_file_query_loaded(self, tmp_path):
        import argparse
        qfile = tmp_path / "q.sql"
        qfile.write_text("SELECT * FROM t")
        ns = argparse.Namespace(query=None, query_file=str(qfile))
        assert bss.load_query(ns) == "SELECT * FROM t"


# ── main() integration ───────────────────────────────────────────────────────

class TestMain:
    def test_basic_sqlite_write(self, tmp_path):
        pytest.importorskip("sqlalchemy")
        db = make_sqlite_db(tmp_path, [
            ("SAMPLE_A", "s3://bucket/a_R1.fastq.gz", "s3://bucket/a_R2.fastq.gz"),
            ("SAMPLE_B", "s3://bucket/b_R1.fastq.gz", "s3://bucket/b_R2.fastq.gz"),
        ])
        out = tmp_path / "sheet.csv"
        sys.argv = [
            "build_samplesheet_from_sql.py",
            "--db-url", f"sqlite:///{db}",
            "--query", "SELECT sample_id, r1, r2 FROM fastq_manifest",
            "--output", str(out),
        ]
        rc = bss.main()
        assert rc == 0
        assert out.exists()
        with out.open() as fh:
            rows = list(csv.DictReader(fh))
        assert len(rows) == 2
        assert rows[0]["sample_id"] == "SAMPLE_A"

    def test_output_has_correct_header(self, tmp_path):
        pytest.importorskip("sqlalchemy")
        db = make_sqlite_db(tmp_path, [("S1", "r1.fastq.gz", "r2.fastq.gz")])
        out = tmp_path / "sheet.csv"
        sys.argv = [
            "build_samplesheet_from_sql.py",
            "--db-url", f"sqlite:///{db}",
            "--query", "SELECT sample_id, r1, r2 FROM fastq_manifest",
            "--output", str(out),
        ]
        bss.main()
        with out.open() as fh:
            header = fh.readline().strip().split(",")
        assert header == ["sample_id", "r1", "r2"]

    def test_custom_column_names(self, tmp_path):
        pytest.importorskip("sqlalchemy")
        import sqlite3
        db = tmp_path / "custom.db"
        conn = sqlite3.connect(db)
        conn.execute("CREATE TABLE t (sid TEXT, read1 TEXT, read2 TEXT)")
        conn.execute("INSERT INTO t VALUES ('S1', 'a.fq', 'b.fq')")
        conn.commit()
        conn.close()
        out = tmp_path / "out.csv"
        sys.argv = [
            "build_samplesheet_from_sql.py",
            "--db-url", f"sqlite:///{db}",
            "--query", "SELECT sid, read1, read2 FROM t",
            "--sample-col", "sid",
            "--r1-col", "read1",
            "--r2-col", "read2",
            "--output", str(out),
        ]
        rc = bss.main()
        assert rc == 0

    def test_empty_query_returns_error(self, tmp_path):
        pytest.importorskip("sqlalchemy")
        db = make_sqlite_db(tmp_path, [])
        out = tmp_path / "sheet.csv"
        sys.argv = [
            "build_samplesheet_from_sql.py",
            "--db-url", f"sqlite:///{db}",
            "--query", "SELECT sample_id, r1, r2 FROM fastq_manifest",
            "--output", str(out),
        ]
        rc = bss.main()
        assert rc == 1

    def test_missing_column_returns_error(self, tmp_path):
        pytest.importorskip("sqlalchemy")
        db = make_sqlite_db(tmp_path, [("S1", "r1.fq", "r2.fq")])
        out = tmp_path / "sheet.csv"
        sys.argv = [
            "build_samplesheet_from_sql.py",
            "--db-url", f"sqlite:///{db}",
            # query returns 'sample_id' but we ask for 'sid' column
            "--query", "SELECT sample_id, r1, r2 FROM fastq_manifest",
            "--sample-col", "sid",
            "--r1-col", "r1",
            "--r2-col", "r2",
            "--output", str(out),
        ]
        rc = bss.main()
        assert rc == 1

    def test_parent_dir_created(self, tmp_path):
        pytest.importorskip("sqlalchemy")
        db = make_sqlite_db(tmp_path, [("S1", "r1.fq", "r2.fq")])
        out = tmp_path / "nested" / "deep" / "sheet.csv"
        sys.argv = [
            "build_samplesheet_from_sql.py",
            "--db-url", f"sqlite:///{db}",
            "--query", "SELECT sample_id, r1, r2 FROM fastq_manifest",
            "--output", str(out),
        ]
        rc = bss.main()
        assert rc == 0
        assert out.exists()

    def test_bad_db_url_returns_nonzero(self, tmp_path):
        pytest.importorskip("sqlalchemy")
        out = tmp_path / "sheet.csv"
        sys.argv = [
            "build_samplesheet_from_sql.py",
            "--db-url", "sqlite:////nonexistent/path/db.sqlite",
            "--query", "SELECT 1",
            "--output", str(out),
        ]
        rc = bss.main()
        assert rc != 0
