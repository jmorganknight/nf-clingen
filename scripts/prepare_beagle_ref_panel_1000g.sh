#!/usr/bin/env bash
set -euo pipefail

# Download 1000 Genomes Phase 3 reference panels and, by default, liftover them
# into GRCh38-compatible VCFs for the genealogy workflow.
#
# Upstream source:
#   /beagle/1000_Genomes_phase3_v5a/b37.vcf/chr{N}.1kg.phase3.v5a.vcf.gz(+.tbi)
#
# Default behavior:
#   source build = b37 (GRCh37)
#   target build = hg38 (GRCh38)
#   output format = VCF/VCF.GZ
#
# This gives a panel that the current genealogy path can consume safely with a
# GRCh38 reference and GRCh38 genetic maps.

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT_DIR="${1:-${ROOT_DIR}/data/beagle/ref_panels}"
TMP_DIR="${ROOT_DIR}/.tmp_1000g_beagle"
BASE_URL="https://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a"
PANEL_BUILD="${PANEL_BUILD:-hg38}"   # hg38 | b37
PANEL_FORMAT="${PANEL_FORMAT:-vcf}"  # vcf | bref3 (bref3 is b37-only here)
PANEL_CHROMS="${PANEL_CHROMS:-1..22,X}"
REF_FASTA="${REF_FASTA:-${ROOT_DIR}/data/GRCh38.fasta}"
CHAIN_FILE="${CHAIN_FILE:-${ROOT_DIR}/data/gatk_bundle/vqsr_hg38/hg19ToHg38.primary.over.chain.gz}"

log() { echo "[beagle-1000g] $*"; }

mkdir -p "${OUT_DIR}" "${TMP_DIR}"

CHR_MAP_FILE="${TMP_DIR}/b37_to_hg19_chr.map"
cat > "${CHR_MAP_FILE}" <<'EOF'
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
MT	chrM
EOF

if [[ "${PANEL_BUILD}" == "hg38" ]]; then
  if [[ ! -f "${REF_FASTA}" ]]; then
    echo "ERROR: GRCh38 reference FASTA not found: ${REF_FASTA}" >&2
    exit 1
  fi
  if [[ ! -f "${CHAIN_FILE}" ]]; then
    echo "ERROR: liftover chain file not found: ${CHAIN_FILE}" >&2
    exit 1
  fi
fi

expand_chroms() {
  local spec="$1"
  local expanded=""
  IFS=',' read -r -a parts <<< "$spec"
  for p in "${parts[@]}"; do
    if [[ "$p" =~ ^([0-9]+)\.\.([0-9]+)$ ]]; then
      for i in $(seq "${BASH_REMATCH[1]}" "${BASH_REMATCH[2]}"); do
        expanded+=" ${i}"
      done
    else
      expanded+=" ${p}"
    fi
  done
  echo "$expanded"
}

chroms=$(expand_chroms "$PANEL_CHROMS")

for chrom in $chroms; do
  src_vcf="${BASE_URL}/b37.vcf/chr${chrom}.1kg.phase3.v5a.vcf.gz"
  src_tbi="${src_vcf}.tbi"
  tmp_src_vcf="${TMP_DIR}/chr${chrom}.1kg.phase3.v5a.b37.vcf.gz"
  tmp_src_tbi="${tmp_src_vcf}.tbi"
  tmp_hg19_vcf="${TMP_DIR}/chr${chrom}.1kg.phase3.v5a.hg19.vcf.gz"

  if [[ "${PANEL_BUILD}" == "b37" && "${PANEL_FORMAT}" == "bref3" ]]; then
    dst="${OUT_DIR}/1000G_phase3.b37.chr${chrom}.bref3"
    if [[ -s "${dst}" ]]; then
      log "chr${chrom}: bref3 already present"
      continue
    fi
    log "chr${chrom}: downloading bref3"
    curl -fL "${BASE_URL}/b37.bref3/chr${chrom}.1kg.phase3.v5a.b37.bref3" -o "${dst}"
    continue
  fi

  dst_vcf="${OUT_DIR}/1000G_phase3.hg38.chr${chrom}.vcf.gz"
  dst_tbi="${dst_vcf}.tbi"
  if [[ -s "${dst_vcf}" && -s "${dst_tbi}" ]]; then
    log "chr${chrom}: hg38 vcf+tbi already present"
    continue
  fi

  log "chr${chrom}: downloading b37 source VCF"
  curl -fL "${src_vcf}" -o "${tmp_src_vcf}"
  curl -fL "${src_tbi}" -o "${tmp_src_tbi}"

  log "chr${chrom}: renaming b37 contigs to hg19-style (chr*)"
  docker run --rm \
    -v "${ROOT_DIR}:${ROOT_DIR}" \
    -w "${ROOT_DIR}" \
    quay.io/biocontainers/bcftools:1.20--h8b25389_0 \
    sh -lc "bcftools annotate --rename-chrs ${CHR_MAP_FILE} -Oz -o ${tmp_hg19_vcf} ${tmp_src_vcf} && tabix -f -p vcf ${tmp_hg19_vcf}"

  log "chr${chrom}: liftover to GRCh38"
  docker run --rm \
    -v "${ROOT_DIR}:${ROOT_DIR}" \
    -w "${ROOT_DIR}" \
    broadinstitute/gatk:4.6.1.0 \
    gatk LiftoverVcf \
      -I "${tmp_hg19_vcf}" \
      -O "${dst_vcf}" \
      -CHAIN "${CHAIN_FILE}" \
      -REJECT "${TMP_DIR}/chr${chrom}.rejected.vcf.gz" \
      -R "${REF_FASTA}" \
      --RECOVER_SWAPPED_REF_ALT true \
      --CREATE_INDEX true \
      --MAX_RECORDS_IN_RAM 100000
done

cat <<EOF
[beagle-1000g] Download complete.
[beagle-1000g] Source build: b37 (GRCh37)
[beagle-1000g] Target build: ${PANEL_BUILD}
[beagle-1000g] Selected format: ${PANEL_FORMAT}
[beagle-1000g] Selected chromosomes: ${PANEL_CHROMS}
[beagle-1000g] Panels are ready for direct use by the genealogy module.
[beagle-1000g] Supported formats in current code path: .bref3 and .vcf.gz
[beagle-1000g] Files:
[beagle-1000g]   ${OUT_DIR}/1000G_phase3.hg38.chr*.vcf.gz(+.tbi)
EOF
