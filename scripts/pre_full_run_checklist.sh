#!/usr/bin/env bash
set -euo pipefail

# Pre-flight checks before launching a full HG002 benchmark run.
# This script is intentionally read-only except for creating no files.

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SUMMARY_CSV="${ROOT_DIR}/results/benchmark_hg002_dv/happy/HG002_giab.summary.csv"
RUNINFO_JSON="${ROOT_DIR}/results/benchmark_hg002_dv/happy/HG002_giab.runinfo.json"
CLINVAR_CHR_VCF="${ROOT_DIR}/data/clinvar/clinvar_GRCh38_chr.vcf.gz"
QUERY_VCF="${ROOT_DIR}/results/benchmark_hg002_dv/variants/HG002_NIST_Illumina2x250_R.deepvariant.vcf.gz"
BAM="${ROOT_DIR}/results/benchmark_hg002_dv/preprocess/HG002_NIST_Illumina2x250_R.sorted.bam"
TARGET_BED="${ROOT_DIR}/data/nextera_expanded_exome_target_regions.bed"

SCRATCH_DIR="${NF_SCRATCH_DIR:-}"
if [[ -z "${SCRATCH_DIR}" ]]; then
  SCRATCH_DIR="${ROOT_DIR}"
fi

echo "== nf-prism pre-full-run checklist =="
echo "workspace: ${ROOT_DIR}"
echo

# 1) GIAB confidence + target masking check via hap.py runinfo
if [[ -f "${RUNINFO_JSON}" ]]; then
  echo "[check] hap.py masking flags in latest runinfo"
  python3 - <<PY
import json
p = "${RUNINFO_JSON}"
with open(p) as fh:
    d = json.load(fh)
fa = d.get("final_args", {})
print("  fp_bedfile:", fa.get("fp_bedfile"))
print("  targets_bedfile:", fa.get("targets_bedfile"))
ok = bool(fa.get("fp_bedfile")) and bool(fa.get("targets_bedfile"))
print("  status:", "OK" if ok else "WARN")
PY
else
  echo "[check] runinfo missing: ${RUNINFO_JSON}"
fi

echo
# 2) Current benchmark metrics snapshot
if [[ -f "${SUMMARY_CSV}" ]]; then
  echo "[check] current mini/full summary metrics"
  awk -F',' 'NR==1 || $1=="SNP" || $1=="INDEL" {print "  "$1","$2",Recall=" $11 ",Precision=" $12 ",F1=" $14}' "${SUMMARY_CSV}"
else
  echo "[check] summary missing: ${SUMMARY_CSV}"
fi

echo
# 3) ClinVar contig normalization readiness
if [[ -f "${CLINVAR_CHR_VCF}" ]]; then
  echo "[check] ClinVar contig-normalized VCF present"
  ls -lh "${CLINVAR_CHR_VCF}" "${CLINVAR_CHR_VCF}.tbi" 2>/dev/null | sed 's/^/  /'
else
  echo "[check] WARN missing contig-normalized ClinVar VCF: ${CLINVAR_CHR_VCF}"
fi

echo
# 4) Basic disk headroom checks
if command -v df >/dev/null 2>&1; then
  echo "[check] disk free space"
  df -h "${ROOT_DIR}" | sed -n '1,2p' | sed 's/^/  /'
  df -h "${SCRATCH_DIR}" | sed -n '1,2p' | sed 's/^/  /'
fi

echo
# 5) Input existence checks for a full run
missing=0
for f in "${QUERY_VCF}" "${BAM}" "${TARGET_BED}"; do
  if [[ ! -f "${f}" ]]; then
    echo "[check] MISSING: ${f}"
    missing=1
  fi
done
if [[ "${missing}" -eq 0 ]]; then
  echo "[check] required local HG002 artifacts present"
fi

echo
# 6) Suggested full-run command (portable, adaptive tuning)
cat <<'CMD'
Suggested full run command:

  export NF_SCRATCH_DIR="/path/to/7TB_scratch"
  export NF_PROFILES="docker,adaptive_local"
  # Optional: tune interactive headroom
  # export NF_RESERVE_CORES="4"
  bash scripts/run_hg002_benchmark.sh full

Notes:
  - adaptive_local is portable; set NF_RESERVE_CORES to auto-cap max_cpus from host CPU count.
  - Leave NF_SCRATCH_DIR unset to keep working files in the project workspace.
CMD

echo
echo "Checklist complete."
