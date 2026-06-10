#!/usr/bin/env bash
set -euo pipefail

# Prepare GRCh38 VQSR resources from the public legacy b37 bundle.
# This is a fallback path for environments where the official hg38 bundle is unavailable.

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUTDIR="${1:-${ROOT_DIR}/data/gatk_bundle/vqsr_hg38}"
REF_FASTA="${2:-${ROOT_DIR}/data/GRCh38.fasta}"
PRIMARY_DICT="${OUTDIR}/tmp/GRCh38.primary.dict"

mkdir -p "${OUTDIR}" "${OUTDIR}/tmp"

if [[ ! -f "${REF_FASTA}" ]]; then
  echo "ERROR: Reference FASTA not found: ${REF_FASTA}" >&2
  exit 1
fi

if [[ ! -f "${REF_FASTA}.fai" ]]; then
  echo "ERROR: Reference FASTA index not found: ${REF_FASTA}.fai" >&2
  exit 1
fi

REF_DICT=""
for candidate in "${REF_FASTA}.dict" "${REF_FASTA%.fasta}.dict" "${REF_FASTA%.fa}.dict"; do
  if [[ -f "${candidate}" ]]; then
    REF_DICT="${candidate}"
    break
  fi
done

if [[ -z "${REF_DICT}" ]]; then
  REF_DICT="${REF_FASTA%.fasta}.dict"
  docker run --rm \
    -v "${ROOT_DIR}:${ROOT_DIR}" \
    -w "${ROOT_DIR}" \
    broadinstitute/gatk:4.6.1.0 \
    gatk CreateSequenceDictionary \
      -R "${REF_FASTA}" \
      -O "${REF_DICT}"
fi

    awk 'BEGIN{OFS="\t"} /^@HD/ {print; next} /^@SQ/ {sn=""; ln=""; for(i=1;i<=NF;i++){if($i ~ /^SN:/){sn=substr($i,4)} else if($i ~ /^LN:/){ln=substr($i,4)}} if (sn ~ /^chr([1-9]$|1[0-9]$|2[0-2]$|X$|Y$|M$)/ && ln != "") print} /^@PG/ {next} /^@CO/ {next}' "${REF_DICT}" > "${PRIMARY_DICT}"

CHAIN_URL="https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz"
CHAIN_FILE="${OUTDIR}/hg19ToHg38.over.chain.gz"
CHAIN_PRIMARY_FILE="${OUTDIR}/hg19ToHg38.primary.over.chain.gz"

if [[ ! -s "${CHAIN_FILE}" ]]; then
  curl -L --fail -o "${CHAIN_FILE}" "${CHAIN_URL}"
fi

zcat "${CHAIN_FILE}" | awk '
  /^chain / {
    keep = ($3 ~ /^chr([1-9]|1[0-9]|2[0-2]|X|Y|M)$/ && $8 ~ /^chr([1-9]|1[0-9]|2[0-2]|X|Y|M)$/)
  }
  keep { print }
' | gzip -c > "${CHAIN_PRIMARY_FILE}"

cat > "${OUTDIR}/tmp/b37_to_hg19_chr.map" <<'EOF'
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

prepare_resource() {
  local source_url="$1"
  local source_name="$2"
  local output_name="$3"

  local source_vcf="${OUTDIR}/tmp/${source_name}"
  local source_norm_txt="${OUTDIR}/tmp/${source_name%.vcf.gz}.norm.vcf"
  local source_norm_vcf="${OUTDIR}/tmp/${source_name%.vcf.gz}.norm.vcf.gz"
  local renamed_vcf="${OUTDIR}/tmp/${source_name%.vcf.gz}.ucsc.vcf.gz"
  local primary_vcf="${OUTDIR}/tmp/${source_name%.vcf.gz}.ucsc.primary.vcf.gz"
  local lifted_vcf="${OUTDIR}/${output_name}"
  local rejected_vcf="${OUTDIR}/tmp/${output_name%.vcf.gz}.rejected.vcf.gz"

  echo "Preparing ${output_name}"
  if [[ ! -s "${source_vcf}" ]]; then
    curl -L --fail -o "${source_vcf}" "${source_url}"
  else
    echo "Using existing source: ${source_vcf}"
  fi

  # Normalize legacy inputs (space-delimited or incomplete contig headers) into strict VCF.
  zcat "${source_vcf}" | awk 'BEGIN{OFS="\t"}
    /^##contig=/ {seen_contig=1; print; next}
    /^##/ {print; next}
    /^#CHROM/ {
      if (!seen_contig) {
        for (i=1;i<=22;i++) print "##contig=<ID=" i ">";
        print "##contig=<ID=X>";
        print "##contig=<ID=Y>";
        print "##contig=<ID=MT>";
      }
      print "#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO";
      next
    }
    !/^#/ {
      gsub(/[ ]+/,"\t");
      print;
      next
    }' > "${source_norm_txt}"

  docker run --rm \
    -v "${OUTDIR}:/data" \
    quay.io/biocontainers/bcftools:1.20--h8b25389_0 \
    sh -lc "bgzip -f -c /data/tmp/$(basename "${source_norm_txt}") > /data/tmp/$(basename "${source_norm_vcf}") && tabix -f -p vcf /data/tmp/$(basename "${source_norm_vcf}")"

  docker run --rm \
    -v "${OUTDIR}:/data" \
    quay.io/biocontainers/bcftools:1.20--h8b25389_0 \
    sh -lc "bcftools annotate --rename-chrs /data/tmp/b37_to_hg19_chr.map -Oz -o /data/tmp/$(basename "${renamed_vcf}") /data/tmp/$(basename "${source_norm_vcf}") && tabix -f -p vcf /data/tmp/$(basename "${renamed_vcf}") && bcftools view -r chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrM -Oz -o /data/tmp/$(basename "${primary_vcf}") /data/tmp/$(basename "${renamed_vcf}") && tabix -f -p vcf /data/tmp/$(basename "${primary_vcf}")"

  rm -f "${lifted_vcf}" "${lifted_vcf}.tbi" "${rejected_vcf}" "${rejected_vcf}.tbi"

  docker run --rm \
    -v "${ROOT_DIR}:${ROOT_DIR}" \
    -v "${OUTDIR}:${OUTDIR}" \
    -w "${ROOT_DIR}" \
    broadinstitute/gatk:4.6.1.0 \
    gatk UpdateVcfSequenceDictionary \
      -I "${primary_vcf}" \
      -SD "${PRIMARY_DICT}" \
      -O "${primary_vcf%.vcf.gz}.headerfix.vcf.gz"

  mv -f "${primary_vcf%.vcf.gz}.headerfix.vcf.gz" "${primary_vcf}"
  if [[ -f "${primary_vcf%.vcf.gz}.headerfix.vcf.gz.tbi" ]]; then
    mv -f "${primary_vcf%.vcf.gz}.headerfix.vcf.gz.tbi" "${primary_vcf}.tbi"
  else
    docker run --rm \
      -v "${OUTDIR}:/data" \
      quay.io/biocontainers/bcftools:1.20--h8b25389_0 \
      sh -lc "tabix -f -p vcf /data/tmp/$(basename "${primary_vcf}")"
  fi

  docker run --rm \
    -v "${ROOT_DIR}:${ROOT_DIR}" \
    -v "${OUTDIR}:${OUTDIR}" \
    -w "${ROOT_DIR}" \
    broadinstitute/gatk:4.6.1.0 \
    gatk LiftoverVcf \
      -I "${primary_vcf}" \
      -O "${lifted_vcf}" \
      -CHAIN "${CHAIN_PRIMARY_FILE}" \
      -REJECT "${rejected_vcf}" \
      -R "${REF_FASTA}" \
      --RECOVER_SWAPPED_REF_ALT true \
      --CREATE_INDEX true \
      --MAX_RECORDS_IN_RAM 100000
}

prepare_resource \
  "https://storage.googleapis.com/gatk-legacy-bundles/b37/hapmap_3.3.b37.vcf.gz" \
  "hapmap_3.3.b37.vcf.gz" \
  "hapmap_3.3.hg38.sites.vcf.gz"

prepare_resource \
  "https://storage.googleapis.com/gatk-legacy-bundles/b37/1000G_omni2.5.b37.vcf.gz" \
  "1000G_omni2.5.b37.vcf.gz" \
  "1000G_omni2.5.hg38.sites.vcf.gz"

prepare_resource \
  "https://storage.googleapis.com/gatk-legacy-bundles/b37/1000G_phase1.snps.high_confidence.b37.vcf.gz" \
  "1000G_phase1.snps.high_confidence.b37.vcf.gz" \
  "1000G_phase1.snps.high_confidence.hg38.vcf.gz"

prepare_resource \
  "https://storage.googleapis.com/gatk-legacy-bundles/b37/Mills_and_1000G_gold_standard.indels.b37.vcf.gz" \
  "Mills_and_1000G_gold_standard.indels.b37.vcf.gz" \
  "Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"

echo
echo "Prepared VQSR resources in ${OUTDIR}:"
ls -1 "${OUTDIR}"/*.vcf.gz 2>/dev/null || true
ls -1 "${OUTDIR}"/*.vcf.gz "${OUTDIR}"/*.vcf.gz.tbi