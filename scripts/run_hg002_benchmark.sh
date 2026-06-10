#!/usr/bin/env bash
set -euo pipefail

# End-to-end HG002 benchmark run:
# 1) Ensure HG002 reads are present (mini by default)
# 2) Run nf-clingen clinical route with DeepVariant
# 3) Score against HG002 GIAB truth with hap.py helper
#
# Usage:
#   scripts/run_hg002_benchmark.sh [mini|full]
#
# Optional environment overrides:
#   NF_PROFILES="docker,adaptive_local" Nextflow profile chain
#   NF_SCRATCH_DIR="/path/to/scratch" Scratch root for Nextflow work dir
#   NF_RESERVE_CORES="4"              Keep this many host cores free (auto-computes max_cpus)

MODE="${1:-mini}"
shift $(( $# > 0 ? 1 : 0 ))
NF_EXTRA_ARGS=("$@")
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUTDIR="${ROOT_DIR}/results/benchmark_hg002_dv"
NF_PROFILES="${NF_PROFILES:-docker,adaptive_local}"
SCRATCH_ARG=()
RESERVE_ARG=()

if [[ -n "${NF_SCRATCH_DIR:-}" ]]; then
  SCRATCH_ARG=(--scratch_dir "${NF_SCRATCH_DIR}")
fi

if [[ -n "${NF_RESERVE_CORES:-}" ]]; then
  if ! [[ "${NF_RESERVE_CORES}" =~ ^[0-9]+$ ]]; then
    echo "ERROR: NF_RESERVE_CORES must be a non-negative integer" >&2
    exit 1
  fi

  if command -v nproc >/dev/null 2>&1; then
    HOST_CPUS="$(nproc --all)"
  else
    HOST_CPUS="$(getconf _NPROCESSORS_ONLN)"
  fi

  TARGET_CPUS=$(( HOST_CPUS - NF_RESERVE_CORES ))
  if (( TARGET_CPUS < 1 )); then
    TARGET_CPUS=1
  fi
  RESERVE_ARG=(--max_cpus "${TARGET_CPUS}")
fi

"${ROOT_DIR}/scripts/prepare_hg002_reads.sh" "${MODE}"

echo "Running Nextflow with profiles: ${NF_PROFILES}"
if [[ -n "${NF_SCRATCH_DIR:-}" ]]; then
  echo "Using scratch directory: ${NF_SCRATCH_DIR}"
fi
if [[ -n "${NF_RESERVE_CORES:-}" ]]; then
  echo "Reserving host CPU cores: ${NF_RESERVE_CORES} (max_cpus=${TARGET_CPUS})"
fi

nextflow run "${ROOT_DIR}" -profile "${NF_PROFILES}" \
  --reads "${ROOT_DIR}/data/benchmark/giab_hg002/HG002_NIST_Illumina2x250_R{1,2}.fastq.gz" \
  --reference "${ROOT_DIR}/data/GRCh38.fasta" \
  --workflow clinical \
  --caller deepvariant \
  --outdir "${OUTDIR}" \
  "${SCRATCH_ARG[@]}" \
  "${RESERVE_ARG[@]}" \
  "${NF_EXTRA_ARGS[@]}"

GIAB_SAMPLE=HG002 bash "${ROOT_DIR}/scripts/run_giab_happy.sh" \
  "${OUTDIR}/variants/HG002_NIST_Illumina2x250_R.deepvariant.vcf.gz" \
  "${OUTDIR}/happy"

echo "HG002 benchmark complete. Key outputs:"
echo "  ${OUTDIR}/variants/HG002_NIST_Illumina2x250_R.deepvariant.vcf.gz"
echo "  ${OUTDIR}/happy/HG002_giab.summary.csv"
