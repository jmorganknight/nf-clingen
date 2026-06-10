#!/usr/bin/env python3
"""Build a Nextflow samplesheet (sample_id,r1,r2) from SQL metadata.

This helper is intentionally generic so teams can adapt it to PostgreSQL,
MySQL, SQLite, or warehouse engines with minimal changes.

Examples:
  python scripts/build_samplesheet_from_sql.py \
    --db-url sqlite:///data/demo_samples.db \
    --query "SELECT sample_id, r1_uri, r2_uri FROM fastq_manifest" \
    --output data/samplesheet.csv

  python scripts/build_samplesheet_from_sql.py \
    --db-url postgresql+psycopg2://user:pass@host:5432/dbname \
    --query-file sql/manifest_query.sql \
    --output s3_manifest.csv
"""

from __future__ import annotations

import argparse
import csv
import pathlib
import sys


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate Nextflow samplesheet CSV from SQL query output"
    )
    parser.add_argument(
        "--db-url",
        required=True,
        help="SQLAlchemy DB URL, e.g. sqlite:///file.db or postgresql+psycopg2://...",
    )
    query_group = parser.add_mutually_exclusive_group(required=True)
    query_group.add_argument("--query", help="Inline SQL query")
    query_group.add_argument("--query-file", help="Path to SQL query file")
    parser.add_argument(
        "--output",
        required=True,
        help="Output CSV path with columns sample_id,r1,r2",
    )
    parser.add_argument(
        "--sample-col",
        default="sample_id",
        help="Column name for sample id in query result (default: sample_id)",
    )
    parser.add_argument(
        "--r1-col",
        default="r1",
        help="Column name for read1 path/URI in query result (default: r1)",
    )
    parser.add_argument(
        "--r2-col",
        default="r2",
        help="Column name for read2 path/URI in query result (default: r2)",
    )
    return parser.parse_args()


def load_query(args: argparse.Namespace) -> str:
    if args.query:
        return args.query
    assert args.query_file
    return pathlib.Path(args.query_file).read_text(encoding="utf-8")


def main() -> int:
    args = parse_args()
    query = load_query(args)

    try:
        from sqlalchemy import create_engine, text
    except Exception as exc:  # pragma: no cover
        print(
            "ERROR: sqlalchemy is required. Install with: pip install sqlalchemy",
            file=sys.stderr,
        )
        print(f"Import error: {exc}", file=sys.stderr)
        return 2

    engine = create_engine(args.db_url)

    try:
        with engine.connect() as conn:
            rows = conn.execute(text(query)).mappings().all()
    except Exception as exc:
        print(f"ERROR: SQL execution failed: {exc}", file=sys.stderr)
        return 1

    if not rows:
        print("ERROR: query returned no rows", file=sys.stderr)
        return 1

    required_cols = [args.sample_col, args.r1_col, args.r2_col]
    first_keys = set(rows[0].keys())
    missing = [col for col in required_cols if col not in first_keys]
    if missing:
        print(
            f"ERROR: query result missing required columns: {', '.join(missing)}",
            file=sys.stderr,
        )
        return 1

    output_path = pathlib.Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with output_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["sample_id", "r1", "r2"])
        for row in rows:
            sample_id = str(row[args.sample_col]).strip()
            r1 = str(row[args.r1_col]).strip()
            r2 = str(row[args.r2_col]).strip()
            if not sample_id or not r1 or not r2:
                print(
                    f"ERROR: empty value in row: sample_id={sample_id!r}, r1={r1!r}, r2={r2!r}",
                    file=sys.stderr,
                )
                return 1
            writer.writerow([sample_id, r1, r2])

    print(f"Wrote samplesheet: {output_path} ({len(rows)} samples)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
