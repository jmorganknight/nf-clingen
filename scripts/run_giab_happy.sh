#!/usr/bin/env bash
set -euo pipefail

# GIAB truth-set scoring using hap.py in Docker.
# Usage:
#   scripts/run_giab_happy.sh [query_vcf] [output_dir]
# If query_vcf is not provided, the script auto-selects the newest VCF in
# results/variants whose filename contains GIAB_SAMPLE (default: HG002).

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
GIAB_SAMPLE="${GIAB_SAMPLE:-HG002}"
TRUTH_VCF="${TRUTH_VCF:-${ROOT_DIR}/data/benchmark/giab_hg002/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz}"
CONF_BED="${CONF_BED:-${ROOT_DIR}/data/benchmark/giab_hg002/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed}"
TARGET_BED="${ROOT_DIR}/data/nextera_expanded_exome_target_regions.bed"
REF_FASTA="${ROOT_DIR}/data/GRCh38.fasta"

QUERY_VCF="${1:-}"
OUTDIR="${2:-${ROOT_DIR}/results/benchmark_eval/happy}"

if [[ -z "${QUERY_VCF}" ]]; then
  QUERY_VCF="$(find "${ROOT_DIR}/results" -type f \( -name '*.gatk.vcf.gz' -o -name '*.deepvariant.vcf.gz' -o -name '*.vcf.gz' \) -path '*/variants/*' -name "*${GIAB_SAMPLE}*" -printf '%T@ %p\n' 2>/dev/null | sort -nr | awk 'NR==1{print $2}')"
fi

if [[ -z "${QUERY_VCF}" ]]; then
  echo "ERROR: No query VCF found for GIAB sample '${GIAB_SAMPLE}'." >&2
  echo "Hint: provide query VCF explicitly as first argument, or set GIAB_SAMPLE to match your VCF naming." >&2
  exit 1
fi

if [[ "${ALLOW_SAMPLE_MISMATCH:-0}" != "1" ]]; then
  query_base="$(basename "${QUERY_VCF}")"
  if [[ "${query_base}" != *"${GIAB_SAMPLE}"* ]]; then
    echo "ERROR: Query VCF '${query_base}' does not match GIAB_SAMPLE='${GIAB_SAMPLE}'." >&2
    echo "Set GIAB_SAMPLE to the query sample name or set ALLOW_SAMPLE_MISMATCH=1 to override." >&2
    exit 1
  fi
fi

for f in "${QUERY_VCF}" "${TRUTH_VCF}" "${CONF_BED}" "${REF_FASTA}"; do
  if [[ ! -f "${f}" ]]; then
    echo "ERROR: Missing required file: ${f}" >&2
    exit 1
  fi
done

mkdir -p "${OUTDIR}"
PREFIX="${OUTDIR}/${GIAB_SAMPLE}_giab"

EVAL_BED="${CONF_BED}"
if [[ -f "${TARGET_BED}" ]]; then
  EVAL_BED="${OUTDIR}/${GIAB_SAMPLE}_giab_eval_regions.bed"
  docker run --rm \
    -v "${ROOT_DIR}:${ROOT_DIR}" \
    -w "${ROOT_DIR}" \
    quay.io/biocontainers/bedtools:2.31.1--hf5e1c6e_2 \
    sh -lc "bedtools intersect -a '${CONF_BED}' -b '${TARGET_BED}' > '${EVAL_BED}'"
fi

if [[ ! -f "${QUERY_VCF}.tbi" ]]; then
  echo "Indexing query VCF with tabix: ${QUERY_VCF}"
  docker run --rm \
    -v "${ROOT_DIR}:${ROOT_DIR}" \
    -w "${ROOT_DIR}" \
    quay.io/biocontainers/bcftools:1.20--h8b25389_0 \
    sh -lc "tabix -f -p vcf '${QUERY_VCF}'"
fi

# Restrict to high-confidence intervals and compare query calls to GIAB truth.
docker run --rm \
  -v "${ROOT_DIR}:${ROOT_DIR}" \
  -w "${ROOT_DIR}" \
  jmcdani20/hap.py:v0.3.12 \
  /opt/hap.py/bin/hap.py \
  "${TRUTH_VCF}" \
  "${QUERY_VCF}" \
  -f "${EVAL_BED}" \
  -r "${REF_FASTA}" \
  -o "${PREFIX}" \
  --threads 8

echo "hap.py complete. Key outputs:"
echo "  ${PREFIX}.summary.csv"
echo "  ${PREFIX}.metrics.json.gz"
echo "  ${PREFIX}.vcf.gz"

if [[ -f "${PREFIX}.summary.csv" ]]; then
  echo
  echo "Summary (all rows):"
  cat "${PREFIX}.summary.csv"
fi
