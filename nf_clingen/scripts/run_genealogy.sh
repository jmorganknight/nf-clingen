#!/usr/bin/env bash
# Run genealogy workflow (Eagle2 phasing → Beagle imputation)
# Pre-requisites: bash scripts/setup_genealogy_resources.sh (downloads genetic maps)
# Optional: Populate data/beagle/ref_panels with per-chromosome reference VCFs
# Usage: bash scripts/run_genealogy.sh [mode] [extra-nextflow-args...]
#   mode options: stub (fast smoke test), mini (small FASTQ), full (production)

set -euo pipefail

MODE="${1:-mini}"
shift $(( $# > 0 ? 1 : 0 ))
NF_EXTRA_ARGS=("$@")

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUTDIR="${ROOT_DIR}/results/genealogy_${MODE}"
NF_PROFILES="${NF_PROFILES:-docker,adaptive_local}"
TEST_INPUT="${ROOT_DIR}/data/test_ref.fa"

# Genealogy resource paths (auto-resolved by setup_genealogy_resources.sh)
EAGLE_MAP_DIR="${ROOT_DIR}/data/eagle2/tables"
BEAGLE_REF_DIR="${ROOT_DIR}/data/beagle/ref_panels"
BEAGLE_MAP_DIR="${ROOT_DIR}/data/beagle/genetic_maps"

# Bootstrap resources by default so the workflow is immediately consumable.
if [[ "${SETUP_GENEALOGY:-true}" == "true" ]]; then
  echo "SETUP_GENEALOGY=true -> running setup_genealogy_resources.sh"
  bash "${ROOT_DIR}/scripts/setup_genealogy_resources.sh"
fi

# Validate critical resource after setup.
if [[ ! -f "${EAGLE_MAP_DIR}/genetic_map_hg38_withX.txt.gz" ]]; then
  echo "ERROR: Eagle2 genetic map missing even after setup. Check network or rerun setup script." >&2
  exit 1
fi

echo "========================================"
echo "  nf-clingen genealogy workflow"
echo "========================================"
echo "Mode:     ${MODE}"
echo "Profiles: ${NF_PROFILES}"
echo "Output:   ${OUTDIR}"
echo "Eagle map: ${EAGLE_MAP_DIR}"
echo "Beagle ref: ${BEAGLE_REF_DIR}"
echo "Beagle map: ${BEAGLE_MAP_DIR}"
echo

# Select input mode
READS_ARG=""
case "${MODE}" in
  stub)
    echo "Running stub-mode smoke test (not a real run)..."
    READS_ARG="--reads ${TEST_INPUT}"
    NF_PROFILES="${NF_PROFILES},test"
    STUB_ARG="-stub-run"
    ;;
  mini)
    echo "Downloading mini HG002 dataset for testing..."
    bash "${ROOT_DIR}/scripts/prepare_hg002_reads.sh" mini
    READS_ARG="--reads ${ROOT_DIR}/data/benchmark/giab_hg002/HG002_NIST_Illumina2x250_R{1,2}.fastq.gz"
    STUB_ARG=""
    ;;
  full)
    echo "Downloading full HG002 dataset (large, slow)..."
    bash "${ROOT_DIR}/scripts/prepare_hg002_reads.sh" full
    READS_ARG="--reads ${ROOT_DIR}/data/benchmark/giab_hg002/HG002_NIST_Illumina2x250_R{1,2}.fastq.gz"
    STUB_ARG=""
    ;;
  *)
    echo "ERROR: Invalid mode '${MODE}'. Use: stub, mini, or full" >&2
    exit 1
    ;;
esac

echo "Starting Nextflow..."
set -x
nextflow run "${ROOT_DIR}" \
  -profile "${NF_PROFILES}" \
  ${STUB_ARG:-} \
  ${READS_ARG} \
  --reference "${ROOT_DIR}/data/GRCh38.fasta" \
  --workflow genealogy \
  --caller haplotypecaller \
  --skip_qc false \
  --eagle2_genetic_map "${EAGLE_MAP_DIR}" \
  --beagle_ref_panel "${BEAGLE_REF_DIR}" \
  --beagle_genetic_map "${BEAGLE_MAP_DIR}" \
  --outdir "${OUTDIR}" \
  "${NF_EXTRA_ARGS[@]}"
