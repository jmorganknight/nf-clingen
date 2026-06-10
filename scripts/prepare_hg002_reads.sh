#!/usr/bin/env bash
set -euo pipefail

# Download and stage HG002 Illumina 2x250 FASTQs from GIAB.
#
# Modes:
#   mini (default): download the first paired chunk only (~12 GB total)
#   full          : download all paired chunks listed in the GIAB index
#
# Usage:
#   scripts/prepare_hg002_reads.sh [mini|full]
#
# Outputs:
#   data/benchmark/giab_hg002/HG002_NIST_Illumina2x250_R1.fastq.gz
#   data/benchmark/giab_hg002/HG002_NIST_Illumina2x250_R2.fastq.gz

MODE="${1:-mini}"
if [[ "${MODE}" != "mini" && "${MODE}" != "full" ]]; then
  echo "ERROR: Mode must be 'mini' or 'full'." >&2
  exit 1
fi

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUTDIR="${ROOT_DIR}/data/benchmark/giab_hg002"
RAWDIR="${OUTDIR}/raw_illumina_2x250"
INDEX_URL="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data_indexes/AshkenazimTrio/sequence.index.AJtrio_Illumina_2x250bps_06012016_updated.HG002"

mkdir -p "${RAWDIR}"

tmp_index="$(mktemp)"
trap 'rm -f "${tmp_index}"' EXIT

curl -fsSL "${INDEX_URL}" -o "${tmp_index}"

pairs_tsv="${RAWDIR}/selected_pairs.tsv"
awk -F '\t' 'NR>1 && NF>=3 {print $1"\t"$3}' "${tmp_index}" > "${pairs_tsv}"

pair_count="$(wc -l < "${pairs_tsv}" | tr -d ' ')"
if [[ "${pair_count}" -eq 0 ]]; then
  echo "ERROR: No FASTQ pairs found in index: ${INDEX_URL}" >&2
  exit 1
fi

if [[ "${MODE}" == "mini" ]]; then
  head -n 1 "${pairs_tsv}" > "${pairs_tsv}.tmp"
  mv "${pairs_tsv}.tmp" "${pairs_tsv}"
fi

echo "Selected $(wc -l < "${pairs_tsv}" | tr -d ' ') HG002 FASTQ pair(s) in mode=${MODE}"

while IFS=$'\t' read -r r1_url r2_url; do
  r1_url="${r1_url/ftp:\/\//https://}"
  r2_url="${r2_url/ftp:\/\//https://}"

  r1_name="$(basename "${r1_url}")"
  r2_name="$(basename "${r2_url}")"

  if [[ ! -s "${RAWDIR}/${r1_name}" ]]; then
    echo "Downloading ${r1_name}"
    curl -fL --retry 3 --retry-delay 2 "${r1_url}" -o "${RAWDIR}/${r1_name}"
  fi

  if [[ ! -s "${RAWDIR}/${r2_name}" ]]; then
    echo "Downloading ${r2_name}"
    curl -fL --retry 3 --retry-delay 2 "${r2_url}" -o "${RAWDIR}/${r2_name}"
  fi
done < "${pairs_tsv}"

r1_out="${OUTDIR}/HG002_NIST_Illumina2x250_R1.fastq.gz"
r2_out="${OUTDIR}/HG002_NIST_Illumina2x250_R2.fastq.gz"

awk -F '\t' '{print $1}' "${pairs_tsv}" | while read -r u; do basename "${u}"; done > "${RAWDIR}/r1.list"
awk -F '\t' '{print $2}' "${pairs_tsv}" | while read -r u; do basename "${u}"; done > "${RAWDIR}/r2.list"

# Concatenate gzip members; resulting .gz remains valid and streams in order.
: > "${r1_out}"
while read -r fn; do
  cat "${RAWDIR}/${fn}" >> "${r1_out}"
done < "${RAWDIR}/r1.list"

: > "${r2_out}"
while read -r fn; do
  cat "${RAWDIR}/${fn}" >> "${r2_out}"
done < "${RAWDIR}/r2.list"

echo "Prepared HG002 read inputs:"
echo "  ${r1_out}"
echo "  ${r2_out}"
ls -lh "${r1_out}" "${r2_out}"
