<p align="center">
  <img src="https://img.shields.io/badge/nextflow-%E2%89%A525.04.3-brightgreen?logo=nextflow&logoColor=white" alt="Nextflow">
  <img src="https://img.shields.io/badge/docker-enabled-2496ED?logo=docker&logoColor=white" alt="Docker">
  <img src="https://img.shields.io/badge/python-3.10%2B-3776AB?logo=python&logoColor=white" alt="Python">
  <img src="https://img.shields.io/badge/tests-41%20passed-brightgreen?logo=pytest" alt="Tests">
  <img src="https://img.shields.io/badge/genome-GRCh38-blueviolet" alt="GRCh38">
  <img src="https://img.shields.io/badge/license-MIT-blue" alt="License">
</p>

<h1 align="center">nf-clingen</h1>

<p align="center">
  <b>Modular Nextflow DSL2 pipeline for clinical variant reporting and genealogy imputation</b><br>
  Containerized · Reproducible · GRCh38 · Works locally or on AWS Batch
</p>

---

## What it does

nf-clingen routes paired-end human sequencing data through a fully containerized stack, branching at variant calling into two downstream endpoints:

| Endpoint | What it produces |
|---|---|
| **`--workflow clinical`** | Annotated VCF, HTML + PDF clinical report, critical variants TSV, run audit trail |
| **`--workflow genealogy`** | Eagle2-phased VCF, Beagle 5.4-imputed GRCh38 VCF, genealogy manifest |

The same upstream stages (QC → alignment → preprocessing → variant calling) feed both endpoints, and every stage is swappable:

```
FASTQ R1/R2
  └─▶ [QC: FastQC + Fastp]        (--skip_qc true to bypass)
        └─▶ [Align: minimap2 | bwamem2]
              └─▶ [Preprocess: samtools | elprep]
                    └─▶ [Call: HaplotypeCaller | DeepVariant]
                          ├─▶ clinical  → annotation + report
                          └─▶ genealogy → Eagle2 phasing → Beagle imputation
```

You can also skip straight to imputation by passing an existing VCF with `--input_vcf`.

---

## Validated performance (HG002, GRCh38)

Evaluated against GIAB v4 truth set, DeepVariant caller, full-genome:

| Type | Recall | Precision | F1 |
|---|---|---|---|
| SNP | **99.40%** | **99.90%** | 0.9965 |
| INDEL | **98.13%** | **99.46%** | 0.9879 |

---

## Quick start

### Requirements
- Java 17+
- Nextflow ≥ 25.04.3
- Docker (or Singularity/Apptainer)
- Reference FASTA: GRCh38 (`data/GRCh38.fasta`)

```bash
# Install Nextflow
curl -s https://get.nextflow.io | bash

# Stub smoke test (no data needed)
nextflow run . -profile test -stub-run --workflow clinical --caller deepvariant
```

### Clinical run

```bash
# Optionally refresh ClinVar (auto-detects GRCh38)
bash scripts/prepare_clinvar_resource.sh

nextflow run . \
  -profile docker,adaptive_local \
  --workflow clinical \
  --caller deepvariant \
  --reference data/GRCh38.fasta \
  --clinvar_vcf data/clinvar/clinvar_GRCh38_chr.vcf.gz \
  --outdir results/clinical_run
```

### Genealogy run

```bash
# Download + liftover 1000G Phase 3 panels to GRCh38 (first run only, ~3 hrs)
bash scripts/setup_genealogy_resources.sh

nextflow run . \
  -profile docker,adaptive_local \
  --workflow genealogy \
  --caller haplotypecaller \
  --reference data/GRCh38.fasta \
  --eagle2_genetic_map data/eagle2/tables \
  --beagle_ref_panel data/beagle/ref_panels \
  --beagle_genetic_map data/beagle/genetic_maps \
  --outdir results/genealogy_run
```

> **Build safety guard**: the pipeline blocks runs where `--reference`, `--eagle2_genetic_map`, and `--beagle_ref_panel` point to different genome builds and exits with a clear error message.

### Input VCF shortcut (skip alignment/calling)

```bash
nextflow run . \
  -profile docker,adaptive_local \
  --workflow genealogy \
  --caller deepvariant \
  --input_vcf results/variants/sample.vcf.gz \
  --reference data/GRCh38.fasta \
  --eagle2_genetic_map data/eagle2/tables \
  --beagle_ref_panel data/beagle/ref_panels \
  --outdir results/genealogy_from_vcf
```

---

## Key parameters

| Parameter | Default | Description |
|---|---|---|
| `--workflow` | `clinical` | Downstream endpoint: `clinical` or `genealogy` |
| `--caller` | `haplotypecaller` | Variant caller: `haplotypecaller` or `deepvariant` |
| `--aligner` | `minimap2` | Aligner: `minimap2` or `bwamem2` |
| `--preprocess` | `samtools` | Preprocessor: `samtools` or `elprep` |
| `--reference` | *(required)* | GRCh38 FASTA path |
| `--input_vcf` | `null` | Skip to endpoint with existing VCF |
| `--clinvar_vcf` | `null` | ClinVar VCF for clinical annotation overlay |
| `--patient_phenotype` | `null` | HPO term (e.g. `HP:0002664`) for phenotype filtering |
| `--eagle2_genetic_map` | `null` | Eagle2 tables directory |
| `--beagle_ref_panel` | `null` | Directory of 1000G GRCh38 per-chromosome VCFs |
| `--beagle_genetic_map` | `null` | Directory of Beagle `.map.gz` files (optional) |
| `--outdir` | `results` | Output directory |
| `--skip_qc` | `false` | Skip FastQC + Fastp |
| `--max_cpus` | auto | CPU ceiling for local executor |
| `--max_ram_gb` | auto | RAM ceiling (GB) for local executor |
| `--scratch_dir` | `/tmp` | Redirect work files to fast scratch |

Full parameter reference: [`nextflow_schema.json`](nextflow_schema.json) · [`params.yaml`](params.yaml)

---

## GRCh38 reference panels for genealogy

The 1000G Phase 3 v5a panels are only published in b37 format. nf-clingen ships a conversion script that automatically downloads and lifts them to GRCh38:

```bash
# First run (downloads + converts chr1–22, X; ~3 hours):
bash scripts/setup_genealogy_resources.sh

# Or run conversion directly for a subset:
PANEL_BUILD=hg38 PANEL_FORMAT=vcf PANEL_CHROMS=1..22,X \
  bash scripts/prepare_beagle_ref_panel_1000g.sh data/beagle/ref_panels
```

**Conversion details:**
1. Download b37 VCF from Beagle host (Washington)
2. `bcftools annotate --rename-chrs` — b37 `1,2,...,X` → hg19-style `chr1,chr2,...,chrX`
3. `GATK LiftoverVcf` with `hg19ToHg38.over.chain.gz` → GRCh38 coordinates
4. `tabix` index output

Validated: **99.7% variant liftover success** on chr1 (7,523 / 2,428,653 unmapped — mostly centromeric gaps).

---

## Clinical phenotype filtering

Filter reported variants to genes relevant to a patient's HPO phenotype:

```bash
nextflow run . \
  -profile docker,adaptive_local \
  --workflow clinical \
  --caller deepvariant \
  --clinvar_vcf data/clinvar/clinvar_GRCh38_chr.vcf.gz \
  --patient_phenotype "HP:0002664"   # Neoplasm
```

Common HPO terms: `HP:0002664` Neoplasm · `HP:0001674` Cardiomyopathy · `HP:0001250` Seizures · `HP:0000365` Hearing impairment

Customize the gene mapping by editing [`data/phenotype_gene_map.tsv`](data/phenotype_gene_map.tsv).

---

## AWS Batch

```bash
nextflow run . \
  -profile awsbatch \
  --execution_mode aws \
  --aws_region us-east-1 \
  --aws_batch_queue nf-clingen-batch-queue \
  --aws_workdir s3://my-nf-bucket/nf-clingen-work \
  --reads 's3://my-raw-bucket/fastq/*_R{1,2}.fastq.gz' \
  --reference s3://my-ref-bucket/GRCh38.fasta \
  --outdir s3://my-nf-bucket/results/run_001 \
  -resume
```

### SQL metadata → samplesheet → pipeline

```bash
python scripts/build_samplesheet_from_sql.py \
  --db-url sqlite:///data/demo_samples.db \
  --query "SELECT sample_id, r1_uri AS r1, r2_uri AS r2 FROM fastq_manifest" \
  --output data/samplesheet.csv

nextflow run . -profile awsbatch \
  --execution_mode aws \
  --samplesheet data/samplesheet.csv \
  --reference s3://my-ref-bucket/GRCh38.fasta \
  --outdir s3://my-nf-bucket/results/cohort_001
```

---

## Configuration

### Quick presets

```bash
# Shared 32-core machine
nextflow run . -profile docker,adaptive_local \
  --max_cpus 30 --deepvariant_cpus 30 --deepvariant_memory '96 GB' \
  --workflow clinical --caller deepvariant

# Dedicated high-RAM host
nextflow run . -profile docker,aggressive_local \
  --workflow clinical --caller deepvariant \
  --clinvar_vcf data/clinvar/clinvar_GRCh38_chr.vcf.gz
```

### params.yaml

Edit [`params.yaml`](params.yaml) to set defaults without repeating flags:

```yaml
workflow: clinical
caller: deepvariant
reference: data/GRCh38.fasta
clinvar_vcf: data/clinvar/clinvar_GRCh38_chr.vcf.gz
max_cpus: 30
scratch_dir: /scratch
```

Then run without any extra flags:

```bash
nextflow run . -profile docker,adaptive_local
```

---

## HG002 benchmark

```bash
bash scripts/prepare_hg002_reads.sh full
bash scripts/prepare_clinvar_resource.sh

nextflow run . \
  -profile docker,adaptive_local \
  --max_cpus 30 --deepvariant_cpus 30 --deepvariant_memory '96 GB' \
  --reads 'data/benchmark/giab_hg002/HG002_NIST_Illumina2x250_R{1,2}.fastq.gz' \
  --reference data/GRCh38.fasta \
  --workflow clinical --caller deepvariant \
  --clinvar_vcf data/clinvar/clinvar_GRCh38_chr.vcf.gz \
  --outdir results/my_benchmark

# Evaluate against GIAB truth
bash scripts/run_giab_happy.sh \
  results/my_benchmark/variants/HG002_NIST_Illumina2x250_R.deepvariant.vcf.gz \
  results/my_benchmark/happy
```

---

## Testing

### Python unit tests (pytest)

```bash
pip install pytest sqlalchemy
pytest                        # runs all 41 tests
pytest tests/ -v --tb=short   # verbose output
```

Test coverage:
- `tests/test_build_clinical_report.py` — clinical report rendering, XSS escaping, critical variant filtering
- `tests/test_build_samplesheet_from_sql.py` — SQL samplesheet generation, error handling, edge cases
- `tests/test_genealogy_build_guard.py` — GRCh38/GRCh37 build inference, compatibility guard logic

### Nextflow stub smoke tests

```bash
# All six workflow routes
for route in \
  "--workflow clinical --caller haplotypecaller --aligner minimap2" \
  "--workflow clinical --caller haplotypecaller --aligner bwamem2" \
  "--workflow clinical --caller deepvariant --aligner minimap2" \
  "--workflow clinical --caller deepvariant --aligner bwamem2" \
  "--workflow genealogy --caller haplotypecaller --aligner minimap2" \
  "--workflow genealogy --caller deepvariant --aligner minimap2"; do
  nextflow run . -profile test -stub-run $route && echo "✓ $route"
done
```

---

## Repository layout

```
nf-clingen/
├── main.nf                          # Orchestration, routing, build guard
├── nextflow.config                  # Profiles, containers, resource tuning
├── nextflow_schema.json             # Parameter schema + validation
├── params.yaml                      # Annotated parameter defaults
├── modules/local/
│   ├── alignment.nf                 # minimap2 / bwamem2
│   ├── clinical.nf                  # bcftools annotation, WeasyPrint report
│   ├── genealogy.nf                 # Eagle2 phasing + Beagle 5.4 imputation
│   ├── phenotype.nf                 # HPO-based gene filtering
│   ├── preprocess.nf                # samtools / elprep
│   ├── qc.nf                        # FastQC + Fastp
│   ├── reference.nf                 # faidx, GATK dict, bwamem2 index
│   └── variant_calling.nf           # HaplotypeCaller + DeepVariant + VQSR
├── subworkflows/local/
│   ├── alignment.nf
│   ├── calling.nf
│   ├── endpoints.nf                 # Clinical / genealogy routing
│   ├── preprocess.nf
│   └── upstream.nf
├── scripts/
│   ├── build_clinical_report.py     # HTML + TSV report generator
│   ├── build_samplesheet_from_sql.py
│   ├── prepare_beagle_ref_panel_1000g.sh  # 1000G b37→hg38 liftover pipeline
│   ├── prepare_clinvar_resource.sh
│   ├── prepare_vqsr_resources.sh
│   ├── setup_genealogy_resources.sh       # First-run resource bootstrap
│   └── run_hg002_benchmark.sh
├── tests/
│   ├── test_build_clinical_report.py
│   ├── test_build_samplesheet_from_sql.py
│   └── test_genealogy_build_guard.py
├── data/
│   ├── GRCh38.fasta / .fai / .gzi
│   ├── phenotype_gene_map.tsv
│   └── benchmark/                   # GIAB truth sets
└── .github/workflows/
    ├── ci-smoke.yml                 # Nextflow stub matrix + pytest
    └── release.yml
```

---

## Pre-flight checks

Run before any long clinical or benchmark job:

```bash
bash scripts/preflight_checks.sh
```

Checks: Java version, Nextflow version, Docker availability, reference files, ClinVar currency, disk space (>50 GB), and memory.

---

## Operational notes

- **Resume after failure**: always use `-resume`. Nextflow caches completed tasks by hash.
- **VQSR**: best for cohorts (≥30 exomes). For single samples, prefer `--caller deepvariant`.
- **Scratch storage**: set `--scratch_dir /path/to/fast/disk` to redirect work files.
- **AWS mode**: local CPU/RAM caps (`max_cpus`, `max_ram_gb`, `max_forks`) are ignored when `--execution_mode aws`.
- **Changing resources mid-run**: resource changes only apply to new tasks. Use `-resume` after config changes.

---

## Limitations

- Single-sample pipeline; multi-sample joint calling is not implemented.
- CNV/SV detection, inheritance modelling, and ACMG criteria scoring are outside scope.
- Population frequency filters and transcript curation are institution-specific responsibilities.
- **Not validated for regulated clinical deployment** — independent validation, quality assurance, and regulatory review are required.

---

## Changelog

See [`CHANGELOG.md`](CHANGELOG.md).

---

## License

MIT — see [`LICENSE`](LICENSE) for details.
