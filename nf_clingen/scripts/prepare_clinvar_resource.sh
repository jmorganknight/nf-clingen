#!/usr/bin/env bash
set -euo pipefail

# Download and prepare ClinVar VCF matching the active reference genome.
# Auto-detects genome version from nextflow.config and downloads matching ClinVar.
# Usage: bash scripts/prepare_clinvar_resource.sh [outdir]

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUTDIR="${1:-${ROOT_DIR}/data/clinvar}"

mkdir -p "${OUTDIR}"

# Read reference from nextflow.config and detect genome version
REF_CONFIG=$(grep "params.reference" "${ROOT_DIR}/nextflow.config" | head -n1 | grep -v "test" || echo "")

if [[ -z "${REF_CONFIG}" ]]; then
  echo "ERROR: Could not determine reference genome from nextflow.config" >&2
  exit 1
fi

# Extract the reference filename to detect genome version
REF_FILE=$(echo "${REF_CONFIG}" | sed -E 's/.*\/([^"]+)".*/\1/')

# Detect genome version from filename
if [[ "${REF_FILE}" =~ GRCh38|hg38 ]]; then
  GRCH_VERSION="GRCh38"
  CLINVAR_VCF_URL="https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz"
  CLINVAR_TBI_URL="https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi"
elif [[ "${REF_FILE}" =~ GRCh37|hg19 ]]; then
  GRCH_VERSION="GRCh37"
  CLINVAR_VCF_URL="https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz"
  CLINVAR_TBI_URL="https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz.tbi"
else
  echo "ERROR: Could not detect genome version from reference: ${REF_FILE}" >&2
  echo "Expected filename to contain GRCh38, hg38, GRCh37, or hg19" >&2
  exit 1
fi

echo "[INFO] Detected reference: ${REF_FILE}"
echo "[INFO] Downloading ClinVar for: ${GRCH_VERSION}"

CLINVAR_VCF="${OUTDIR}/clinvar_${GRCH_VERSION}.vcf.gz"
CLINVAR_TBI="${OUTDIR}/clinvar_${GRCH_VERSION}.vcf.gz.tbi"

CLINVAR_CHR_VCF="${OUTDIR}/clinvar_${GRCH_VERSION}_chr.vcf.gz"
CLINVAR_CHR_TBI="${OUTDIR}/clinvar_${GRCH_VERSION}_chr.vcf.gz.tbi"

echo "[INFO] Downloading ClinVar VCF for ${GRCH_VERSION}..."
curl -L --fail -o "${CLINVAR_VCF}" "${CLINVAR_VCF_URL}"
curl -L --fail -o "${CLINVAR_TBI}" "${CLINVAR_TBI_URL}"

echo "[INFO] Creating chr-prefixed copy (handles both chr1 and 1 naming conventions)..."
bcftools annotate \
  -I +'%CHROM' \
  --rename-chrs <(cat <<'CHRMAP'
1	chr1
2	chr2
3	chr3
4	chr4
5	chr5
6	chr6
7	chr7
8	chr8
9	chr9
10	chr10
11	chr11
12	chr12
13	chr13
14	chr14
15	chr15
16	chr16
17	chr17
18	chr18
19	chr19
20	chr20
21	chr21
22	chr22
X	chrX
Y	chrY
MT	chrMT
CHRMAP
) \
  "${CLINVAR_VCF}" \
  -O z -o "${CLINVAR_CHR_VCF}"

tabix -f -p vcf "${CLINVAR_CHR_VCF}"

echo "[INFO] ClinVar resources prepared:"
echo "  - ${CLINVAR_VCF}"
echo "  - ${CLINVAR_TBI}"
echo "  - ${CLINVAR_CHR_VCF}"
echo "  - ${CLINVAR_CHR_TBI}"
echo ""
echo "[INFO] Use in pipeline:"
echo "  --clinvar_vcf ${CLINVAR_CHR_VCF}"
echo ""
echo "[INFO] Last modified on NCBI: $(curl -s -I ${CLINVAR_VCF_URL} | grep -i last-modified || echo 'unknown')"
