#!/usr/bin/env bash
# Comprehensive preflight checks before running nf-clingen
# Validates: Java, Nextflow, Docker, reference files, disk space, memory, parameters
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ERRORS=()
WARNINGS=()

function log_info() {
  echo "[INFO] $*"
}

function log_pass() {
  echo "✓ $*"
}

function log_warn() {
  echo "⚠ $*"
  WARNINGS+=("$*")
}

function log_error() {
  echo "✗ $*"
  ERRORS+=("$*")
}

echo "========================================"
echo "  nf-clingen Preflight Validation"
echo "========================================"
echo "workspace: ${ROOT_DIR}"
echo

log_info "Checking runtime environment..."
echo

# 1. Java
if command -v java &>/dev/null; then
  java_version=$(java -version 2>&1 | head -n1)
  log_pass "Java installed: $java_version"
else
  log_error "Java not found. Install Java 11+ (required for Nextflow)"
fi

# 2. Nextflow
if command -v nextflow &>/dev/null; then
  nf_version=$(nextflow -version 2>&1 | head -n1)
  log_pass "Nextflow installed: $nf_version"
else
  log_error "Nextflow not found. Run: curl -fsSL https://get.nextflow.io | bash"
fi

# 3. Docker or alternative container runtime
if command -v docker &>/dev/null; then
  docker_version=$(docker --version)
  log_pass "Docker installed: $docker_version"
  if ! docker ps &>/dev/null; then
    log_warn "Docker daemon may not be running. Test with: docker ps"
  fi
elif command -v singularity &>/dev/null; then
  log_pass "Singularity available as alternative container runtime"
elif command -v apptainer &>/dev/null; then
  log_pass "Apptainer available as alternative container runtime"
else
  log_error "No container runtime found (Docker, Singularity, or Apptainer required)"
fi

echo
log_info "Checking Nextflow configuration..."
echo

# 4. Nextflow config validation
if nextflow config &>/dev/null; then
  log_pass "Nextflow config validates"
else
  log_error "Nextflow config has syntax errors"
fi

# 5. Parameter validation
if nextflow run "$ROOT_DIR" --validate_only true -profile docker &>/dev/null; then
  log_pass "Pipeline parameters validate"
else
  log_warn "Parameter validation failed. Review: params.yaml, nextflow_schema.json"
fi

echo
log_info "Checking input data and resources..."
echo

# 6. Reference genome
if [[ -f "${ROOT_DIR}/data/reference.fasta" ]]; then
  ref_size=$(du -h "${ROOT_DIR}/data/reference.fasta" | cut -f1)
  log_pass "Reference genome found: $ref_size"
else
  log_warn "Reference genome not found at ${ROOT_DIR}/data/reference.fasta"
fi

# 7. Input FASTQ data
fastq_pattern="${1:-${ROOT_DIR}/data/*_R{1,2}.fastq.gz}"
fastq_count=$(eval "ls -1 $fastq_pattern 2>/dev/null | wc -l" || echo "0")
if [[ $fastq_count -gt 0 ]]; then
  log_pass "Input FASTQ files found: $fastq_count files"
else
  log_warn "No FASTQ files found matching: $fastq_pattern"
fi

# 8. ClinVar (if using clinical workflow)
if [[ -f "${ROOT_DIR}/data/clinvar/clinvar_GRCh38_chr.vcf.gz" ]]; then
  clinvar_size=$(du -h "${ROOT_DIR}/data/clinvar/clinvar_GRCh38_chr.vcf.gz" | cut -f1)
  log_pass "ClinVar resource found: $clinvar_size"
else
  log_warn "ClinVar resource missing. To prepare, run: bash ${ROOT_DIR}/scripts/prepare_clinvar_resource.sh"
fi

# 9. VQSR resources (if using HaplotypeCaller)
vqsr_path="${ROOT_DIR}/data/gatk_bundle/vqsr_hg38"
if [[ -d "$vqsr_path" ]]; then
  vqsr_count=$(find "$vqsr_path" -name "*.vcf.gz" | wc -l)
  log_pass "VQSR resources found: $vqsr_count VCF files"
else
  log_warn "VQSR resources missing. To prepare, run: bash ${ROOT_DIR}/scripts/prepare_vqsr_resources.sh"
fi

echo
log_info "Checking disk space and memory..."
echo

# 10. Disk space (work directory)
work_dir="${ROOT_DIR}/work"
if [[ -d "$work_dir" ]]; then
  work_size=$(du -sh "$work_dir" 2>/dev/null | cut -f1)
  log_info "Existing work directory: $work_size"
fi

available_disk=$(df -h "$ROOT_DIR" | tail -n1 | awk '{print $4}')
log_info "Available disk on workspace: $available_disk"
if [[ $available_disk > 50G ]]; then
  log_pass "Sufficient disk space (>50 GB available)"
else
  log_warn "Limited disk space. Full pipeline may require 50+ GB. Consider --scratch_dir"
fi

# 11. Memory
if command -v free &>/dev/null; then
  total_mem=$(free -h | grep "^Mem:" | awk '{print $2}')
  available_mem=$(free -h | grep "^Mem:" | awk '{print $7}')
  log_info "Total system memory: $total_mem (available: $available_mem)"
  if [[ $available_mem > 32G ]]; then
    log_pass "Sufficient memory for DeepVariant (≥32 GB available)"
  else
    log_warn "Limited memory. DeepVariant may require --deepvariant_memory tuning"
  fi
fi

echo
log_info "Checking pipeline structure..."
echo

# 12. Critical files
for file in main.nf nextflow.config nextflow_schema.json params.yaml; do
  if [[ -f "${ROOT_DIR}/${file}" ]]; then
    log_pass "Found: $file"
  else
    log_error "Missing critical file: $file"
  fi
done

echo
log_info "Suggestion: Run sample command..."
echo

echo "Example clinical run:"
echo "  nextflow run . -profile docker,adaptive_local \\"
echo "    --workflow clinical --caller deepvariant \\"
echo "    --max_cpus 30 --deepvariant_cpus 30 \\"
echo "    --deepvariant_memory '96 GB' -resume"
echo

echo "========================================"
echo "  Preflight Summary"
echo "========================================"

if [[ ${#ERRORS[@]} -gt 0 ]]; then
  echo "Errors (${#ERRORS[@]}):"
  for e in "${ERRORS[@]}"; do
    echo "  ✗ $e"
  done
  echo
fi

if [[ ${#WARNINGS[@]} -gt 0 ]]; then
  echo "Warnings (${#WARNINGS[@]}):"
  for w in "${WARNINGS[@]}"; do
    echo "  ⚠ $w"
  done
  echo
fi

if [[ ${#ERRORS[@]} -eq 0 ]]; then
  log_pass "All critical checks passed!"
  exit 0
else
  log_error "Fix errors above before running pipeline"
  exit 1
fi
