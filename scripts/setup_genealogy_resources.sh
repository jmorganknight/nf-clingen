#!/usr/bin/env bash
set -euo pipefail

# Prepare local directories and map resources for the genealogy workflow.
# This script does not download reference panels because those are large and
# source-specific; it creates the expected directory layout and validates what
# is present so runs can be prepared quickly.

# Optional: set DOWNLOAD_1000G_BUILD=hg38|b37 (default: hg38).
# Optional: set DOWNLOAD_1000G_FORMAT=vcf|bref3 (default: vcf).

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
EAGLE_DIR="${ROOT_DIR}/data/eagle2/tables"
BEAGLE_MAP_DIR="${ROOT_DIR}/data/beagle/genetic_maps"
BEAGLE_REF_DIR="${ROOT_DIR}/data/beagle/ref_panels"
TMP_DIR="${ROOT_DIR}/.tmp_genealogy_setup"

EAGLE_MAP_URL="https://alkesgroup.broadinstitute.org/Eagle/downloads/tables/genetic_map_hg38_withX.txt.gz"
BEAGLE_MAP_ZIP_URL="https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh38.map.zip"

log() { echo "[genealogy-setup] $*"; }

mkdir -p "${EAGLE_DIR}" "${BEAGLE_MAP_DIR}" "${BEAGLE_REF_DIR}" "${TMP_DIR}"

if [[ ! -s "${EAGLE_DIR}/genetic_map_hg38_withX.txt.gz" ]]; then
  log "Downloading Eagle2 genetic map (GRCh38)..."
  curl -fL "${EAGLE_MAP_URL}" -o "${EAGLE_DIR}/genetic_map_hg38_withX.txt.gz"
else
  log "Eagle2 genetic map already present."
fi

# Download and normalize Beagle map files to the naming expected by the pipeline:
#   *.chr{N}.map.gz (for 1..22 and X)
if ! ls "${BEAGLE_MAP_DIR}"/*.chr1.map.gz >/dev/null 2>&1; then
  log "Downloading Beagle GRCh38 genetic maps..."
  rm -rf "${TMP_DIR}"/*
  curl -fL "${BEAGLE_MAP_ZIP_URL}" -o "${TMP_DIR}/plink.GRCh38.map.zip"
  unzip -o "${TMP_DIR}/plink.GRCh38.map.zip" -d "${TMP_DIR}/maps" >/dev/null

  SRC_DIR="${TMP_DIR}/maps/no_chr_in_chrom_field"
  if [[ ! -d "${SRC_DIR}" ]]; then
    log "ERROR: Unexpected map archive layout: ${SRC_DIR} missing"
    exit 1
  fi

  for chrom in {1..22} X; do
    src="${SRC_DIR}/plink.chr${chrom}.GRCh38.map"
    dst="${BEAGLE_MAP_DIR}/plink.GRCh38.chr${chrom}.map.gz"
    if [[ -s "${src}" ]]; then
      gzip -c "${src}" > "${dst}"
    fi
  done
  log "Beagle genetic maps prepared in ${BEAGLE_MAP_DIR}."
else
  log "Beagle genetic maps already present."
fi

# Validate what exists for reference panels.
panel_hg38_count=$(find "${BEAGLE_REF_DIR}" -maxdepth 1 -type f \( -name '*hg38*.vcf.gz' -o -name '*hg38*.bref3' \) | wc -l | tr -d ' ')
panel_b37_count=$(find "${BEAGLE_REF_DIR}" -maxdepth 1 -type f \( -name '*b37*.vcf.gz' -o -name '*b37*.bref3' \) | wc -l | tr -d ' ')
map_count=$(find "${BEAGLE_MAP_DIR}" -maxdepth 1 -type f -name '*.map.gz' | wc -l | tr -d ' ')

log "Summary"
log "  Eagle map:   ${EAGLE_DIR}/genetic_map_hg38_withX.txt.gz"
log "  Beagle maps: ${map_count} files in ${BEAGLE_MAP_DIR}"
log "  Ref panels:  ${panel_hg38_count} hg38 files in ${BEAGLE_REF_DIR}"
log "  Ref panels:  ${panel_b37_count} b37 files in ${BEAGLE_REF_DIR}"

if [[ "${panel_hg38_count}" -eq 0 ]]; then
  cat <<EOF
[genealogy-setup] WARNING: No Beagle reference panel files found.
[genealogy-setup] The workflow can start, but imputation will skip chromosomes
[genealogy-setup] without matching files in:
[genealogy-setup]   ${BEAGLE_REF_DIR}
[genealogy-setup] Expected naming examples:
[genealogy-setup]   panel.hg38.chr22.vcf.gz
[genealogy-setup]   panel.hg38.chr22.vcf.gz.tbi
EOF

  if [[ "${DOWNLOAD_1000G_BUILD:-auto}" == "auto" || "${DOWNLOAD_1000G_BUILD:-hg38}" == "hg38" ]]; then
    log "No GRCh38 panels found -> auto-downloading GRCh38 1000G VCF panels"
    PANEL_BUILD=hg38 PANEL_FORMAT=vcf bash "${ROOT_DIR}/scripts/prepare_beagle_ref_panel_1000g.sh" "${BEAGLE_REF_DIR}"
    panel_hg38_count=$(find "${BEAGLE_REF_DIR}" -maxdepth 1 -type f \( -name '*hg38*.vcf.gz' -o -name '*hg38*.bref3' \) | wc -l | tr -d ' ')
    log "Ref panels:  ${panel_hg38_count} hg38 files in ${BEAGLE_REF_DIR}"
  fi
fi

if [[ "${DOWNLOAD_1000G_BUILD:-}" == "b37" ]]; then
  log "DOWNLOAD_1000G_BUILD=b37 -> downloading 1000G b37 panels"
  PANEL_BUILD=b37 PANEL_FORMAT=${DOWNLOAD_1000G_FORMAT:-vcf} bash "${ROOT_DIR}/scripts/prepare_beagle_ref_panel_1000g.sh" "${BEAGLE_REF_DIR}"
fi

cat <<EOF
[genealogy-setup] Ready-to-run command template:
nextflow run . -profile docker,adaptive_local \
  --workflow genealogy \
  --caller haplotypecaller \
  --eagle2_genetic_map ${EAGLE_DIR} \
  --beagle_ref_panel ${BEAGLE_REF_DIR} \
  --beagle_genetic_map ${BEAGLE_MAP_DIR} \
  -resume
EOF
