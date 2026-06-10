#!/usr/bin/env bash
set -euo pipefail

# Run hap.py with clinically relevant exome evaluation boundaries.

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
QUERY_VCF="${1:-}"
OUTDIR="${2:-${ROOT_DIR}/results/benchmark_clinical_eval/happy}"

if [[ -z "${QUERY_VCF}" ]]; then
  echo "ERROR: Provide a query VCF as the first argument." >&2
  exit 1
fi

bash "${ROOT_DIR}/scripts/run_giab_happy.sh" "${QUERY_VCF}" "${OUTDIR}"