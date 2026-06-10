# nf-clingen

nf-clingen is a modular Nextflow DSL2 router for paired-end human sequencing runs that need to switch between clinical reporting and genealogy endpoints while keeping the upstream stack containerized and low on unnecessary intermediate disk I/O.

## Layout

- `main.nf`: orchestration, schema validation, runtime logging, and routing logic.
- `subworkflows/local/`: orchestration layers that map route decisions to module processes.
- `nextflow.config`: defaults, profiles, and process-level containers/resources.
- `modules/local/`: isolated process modules for QC, reference preparation, alignment, preprocessing, calling, and downstream endpoints.
- `nextflow_schema.json`: JSON schema used for parameter validation and generated help text.
- `.github/workflows/ci-smoke.yml`: matrix smoke checks running `-stub-run` against test data.

## Key behaviors

- Runs `RUN_FASTQC` and `RUN_FASTP` unless `--skip_qc true` is set.
- Builds reference indexes inside the workflow for faidx, GATK dictionary, and bwa-mem2 when required.
- Routes alignment with `--aligner minimap2|bwamem2`.
- Routes preprocessing with `--preprocess samtools|elprep`.
- Routes calling with `--caller haplotypecaller|deepvariant`.
- Routes downstream outputs with `--workflow clinical|genealogy`.
- Uses a concrete local clinical annotation stack based on `bcftools` triage with optional ClinVar overlay and WeasyPrint PDF compilation.

## Quick start

```bash
nextflow run . -profile docker
```

AWS Batch quick start (managed CPU/memory scheduling):

```bash
nextflow run . \
  -profile awsbatch \
  --execution_mode aws \
  --aws_region us-east-1 \
  --aws_batch_queue nf-clingen-batch-queue \
  --aws_workdir s3://my-nf-bucket/nf-clingen-work \
  --reads 's3://my-raw-bucket/fastq/*_R{1,2}.fastq.gz' \
  --reference s3://my-ref-bucket/GRCh38.fasta \
  --outdir s3://my-nf-bucket/results/clinical_run
```

Clinical run with DeepVariant and ClinVar overlay:

```bash
nextflow run . \
  -profile docker \
  --workflow clinical \
  --caller deepvariant \
  --clinvar_vcf /path/to/clinvar.vcf.gz
```

Auto-download latest ClinVar (auto-detects reference genome version):

```bash
# Download/update ClinVar from NCBI (auto-detects GRCh version from nextflow.config)
bash scripts/prepare_clinvar_resource.sh

# Then run pipeline with the downloaded resource
nextflow run . \
  -profile docker \
  --workflow clinical \
  --caller deepvariant \
  --clinvar_vcf data/clinvar/clinvar_GRCh38_chr.vcf.gz
```

**Important:** The script reads `nextflow.config` to detect the active reference genome (GRCh38 or GRCh37) and automatically downloads the matching ClinVar build. If you change the reference in the config, just re-run this script.

Performance-tuned run with adaptive local resource sizing:

```bash
nextflow run . \
  -profile docker,adaptive_local \
  --scratch_dir /path/to/7TB_scratch \
  --max_cpus 30 \
  --workflow clinical \
  --caller deepvariant \
  --deepvariant_cpus 30 \
  --deepvariant_memory '96 GB'
```

Notes:

- `adaptive_local` is hardware-portable and can cap local CPU usage via `--max_cpus`.
- For large shared-memory hosts, use `-profile docker,aggressive_local` to push harder on CPU and RAM per step.
- DeepVariant resources are tunable at runtime:
  - `--deepvariant_cpus` (default `16`)
  - `--deepvariant_memory` (default `'48 GB'`)
  - `--deepvariant_time` (default `'24h'`)
- For shared hosts, set both `--max_cpus` and `--deepvariant_cpus` to reserve CPU headroom for other tasks.
- `--scratch_dir` routes Nextflow working files to fast scratch (`<scratch_dir>/nf-clingen-work`).
- In AWS mode (`-profile awsbatch --execution_mode aws`), local host caps (`max_cpus`, `max_ram_gb`, `max_forks`, `scratch_dir`) are ignored.

## Configuration with YAML (Optional)

For users who prefer configuration files over CLI flags, **nf-clingen** supports an optional `params.yaml` that can act as both a runnable defaults file and a descriptive parameter reference.

### Using `params.yaml`

Instead of passing all parameters via command-line, you can configure them in a `params.yaml` file in the project root:

```bash
# Edit params.yaml to configure pipeline
nano params.yaml

# Then run Nextflow (it auto-loads params.yaml):
nextflow run . -profile docker,aggressive_local
```

**Example `params.yaml`:**

```yaml
# Workflow routing
workflow: clinical
caller: deepvariant
aligner: minimap2
preprocess: samtools

# Inputs and outputs
clinvar_vcf: data/clinvar/clinvar_GRCh38_chr.vcf.gz
outdir: results/clinical_run

# Resources
max_cpus: 62
max_forks: 15
deepvariant_cpus: 48

# Scratch storage (optional)
scratch_dir: /scratch
```

The repository `params.yaml` is intentionally annotated as a living template. It lists the supported options, valid values, and practical notes so teams can use it as documentation instead of keeping a separate cheat sheet.

**Precedence (highest to lowest):**
1. **CLI flags** (e.g., `--param value`) — overrides everything
2. **params.yaml** — overrides nextflow.config defaults
3. **nextflow.config** defaults — auto-detects host RAM/CPU

Execution mode behavior:

- `execution_mode: local` (default): uses host auto-detection and local tuning parameters.
- `execution_mode: aws`: uses AWS Batch scheduling and S3 work directory (`aws_workdir`).
- Keep `params.yaml` focused on user-facing parameters; infrastructure toggle remains via profile + mode.

A sample annotated `params.yaml` is provided in the repo; uncomment and adjust any parameters you need.

## Resource tuning

Use these knobs to tune throughput vs interactive headroom:

- `--max_cpus`: global CPU ceiling for local executor under `adaptive_local`.
- `--max_forks`: max concurrent processes (default `4`).
- `--deepvariant_cpus`: DeepVariant shard count/CPU allocation (default `16`).
- `--deepvariant_memory`: DeepVariant container memory (default `'48 GB'`).
- `--deepvariant_time`: DeepVariant walltime budget (default `'24h'`).
- `--scratch_dir`: place work files on fast storage (`<scratch_dir>/nf-clingen-work`).

These knobs are local-mode only. In AWS mode, Batch schedules resources based on each process `cpus`/`memory` directives.

## AWS Execution

nf-clingen supports AWS Batch execution with a single profile toggle.

### Required AWS resources

- S3 bucket(s) for input, reference, work, and output paths
- AWS Batch compute environment and job queue
- IAM permissions for S3, Batch, CloudWatch logs, and container registry pulls

### Minimal run

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

### SQL metadata -> samplesheet -> pipeline

You can drive inputs from SQL by creating a samplesheet and passing `--samplesheet`.

1. Build samplesheet from SQL metadata:

```bash
python scripts/build_samplesheet_from_sql.py \
  --db-url sqlite:///data/demo_samples.db \
  --query "SELECT sample_id, r1_uri AS r1, r2_uri AS r2 FROM fastq_manifest" \
  --output data/samplesheet.csv
```

2. Run pipeline with samplesheet (reads glob ignored):

```bash
nextflow run . \
  -profile awsbatch \
  --execution_mode aws \
  --aws_workdir s3://my-nf-bucket/nf-clingen-work \
  --samplesheet data/samplesheet.csv \
  --reference s3://my-ref-bucket/GRCh38.fasta \
  --outdir s3://my-nf-bucket/results/run_from_sql
```

Practical presets:

- Shared 32-core machine (leave room for other tasks):

```bash
nextflow run . \
  -profile docker,adaptive_local \
  --max_cpus 30 \
  --deepvariant_cpus 30 \
  --deepvariant_memory '96 GB' \
  --workflow clinical \
  --caller deepvariant
```

- High-throughput 64-core / high-RAM host:

```bash
nextflow run . \
  -profile docker,aggressive_local \
  --workflow clinical \
  --caller deepvariant \
  --clinvar_vcf data/clinvar/clinvar_GRCh38_chr.vcf.gz
```

- Dedicated host (use most available CPU):

```bash
nextflow run . \
  -profile docker,adaptive_local \
  --max_cpus 120 \
  --deepvariant_cpus 60 \
  --deepvariant_memory '96 GB' \
  --workflow clinical \
  --caller deepvariant
```

- Wrapper-script run with forwarded Nextflow args:

```bash
bash scripts/run_hg002_benchmark.sh full \
  --max_cpus 30 \
  --deepvariant_cpus 30 \
  --deepvariant_memory '96 GB' \
  -resume
```

Important:

- Changing `cpus`/`memory` for a task does not modify an already running task.
- To apply new resource values, relaunch with `-resume`; cached completed steps are reused.

Genealogy run with QC disabled:

```bash
nextflow run . \
  -profile docker \
  --workflow genealogy \
  --skip_qc true
```

**Clinical run with patient phenotype filtering:**

```bash
nextflow run . \
  -profile docker \
  --workflow clinical \
  --caller haplotypecaller \
  --patient_phenotype "HP:0002664"
```

The `patient_phenotype` parameter filters reported variants to genes relevant to the patient's phenotype/disease. Phenotypes use HPO (Human Phenotype Ontology) term IDs. Common examples:

- `HP:0002664` - Neoplasm (cancer)
- `HP:0001674` - Heart murmur (cardiovascular)
- `HP:0001250` - Seizures (neurological)
- `HP:0000365` - Hearing impairment

Customize the gene mapping by editing [data/phenotype_gene_map.tsv](data/phenotype_gene_map.tsv).

Validate parameters only:

```bash
nextflow run . --validate_only true
```

Prepare hg38 VQSR resources from the public legacy bundle:

```bash
bash scripts/prepare_vqsr_resources.sh
```

Run HaplotypeCaller with explicit VQSR resources:

```bash
nextflow run . \
  -profile docker \
  --workflow clinical \
  --caller haplotypecaller \
  --vqsr_snp_resource data/gatk_bundle/vqsr_hg38/hapmap_3.3.hg38.sites.vcf.gz \
  --vqsr_omni_resource data/gatk_bundle/vqsr_hg38/1000G_omni2.5.hg38.sites.vcf.gz \
  --vqsr_1kg_resource data/gatk_bundle/vqsr_hg38/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
  --vqsr_indel_resource data/gatk_bundle/vqsr_hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
```

Smoke-test the full graph with synthetic test inputs:

```bash
nextflow run . -profile test -stub-run --workflow clinical --caller haplotypecaller
```

Additional smoke routes:

```bash
nextflow run . -profile test -stub-run --workflow genealogy --caller haplotypecaller
nextflow run . -profile test -stub-run --workflow clinical --caller deepvariant
```

## HG002 benchmark run

Prepare HG002 Illumina 2x250 paired reads and run DeepVariant + hap.py against HG002 truth:

```bash
bash scripts/run_hg002_benchmark.sh mini
```

Recommended full run with adaptive profile:

```bash
export NF_PROFILES="docker,adaptive_local"
export NF_SCRATCH_DIR="/path/to/7TB_scratch"
export NF_RESERVE_CORES="4"  # optional
bash scripts/run_hg002_benchmark.sh full
```

If you want more interactive headroom on busy systems, run Nextflow directly and cap CPU usage:

```bash
nextflow run . \
  -profile docker,adaptive_local \
  --max_cpus 28 \
  --deepvariant_cpus 24 \
  --workflow clinical \
  --caller deepvariant
```

When using the HG002 benchmark helper script, extra args are forwarded to `nextflow run`, so this works:

```bash
bash scripts/run_hg002_benchmark.sh full --max_cpus 30 --deepvariant_cpus 30 --deepvariant_memory '96 GB' -resume
```

Run pre-flight checks before a long full-mode run:

```bash
bash scripts/pre_full_run_checklist.sh
```

Notes:

- `mini` downloads one HG002 read-pair chunk (faster setup, lower depth).
- `full` downloads all chunks from the GIAB index (slow, large disk footprint).
- You can run prep and benchmark as separate steps:

```bash
bash scripts/prepare_hg002_reads.sh mini
nextflow run . -profile docker \
  --reads "${PWD}/data/benchmark/giab_hg002/HG002_NIST_Illumina2x250_R{1,2}.fastq.gz" \
  --reference "${PWD}/data/GRCh38.fasta" \
  --workflow clinical \
  --caller deepvariant \
  --outdir "${PWD}/results/benchmark_hg002_dv"

GIAB_SAMPLE=HG002 bash scripts/run_giab_happy.sh \
  "${PWD}/results/benchmark_hg002_dv/variants/HG002_NIST_Illumina2x250_R.deepvariant.vcf.gz" \
  "${PWD}/results/benchmark_hg002_dv/happy"
```

## Current status (June 2026)

- End-to-end clinical workflow execution fully operational with `--caller deepvariant` and ClinVar overlay.
- Clinical reporting completes successfully (HTML + PDF + annotated VCF + clinical summary outputs).
- Single authoritative benchmark run available: **HG002 full-genome DeepVariant vs GIAB v4 truth** with 99.4% SNP recall / 99.9% precision.

Primary result locations:

- `results/benchmark_clinvar_dv/clinical/annotations/` — annotated clinical VCF
- `results/benchmark_clinvar_dv/clinical/reports/` — HTML and PDF clinical reports
- `results/benchmark_clinvar_dv/bench_eval/` — hap.py benchmark metrics (summary.csv, ROC curves, JSON metrics)
- `results/benchmark_clinvar_dv/variants/` — raw and ClinVar-overlaid VCFs

## Workflow sketch

```text
FASTQ R1/R2
  |
  v
[QC: FASTQC + FASTP]  (optional via --skip_qc)
  |
  v
[Alignment: minimap2 | bwamem2]
  |
  v
[Preprocess: samtools | elprep]
  |
  v
[Calling]
  |- DeepVariant -------------------------------> Clinical annotation/report
  |
  `- HaplotypeCaller
      |- Hard filter (baseline) --------------> Benchmark
      `- VQSR (trained recalibration model) --> Clinical annotation/report + Benchmark

Downstream endpoint routing:
  --workflow clinical  -> bcftools triage + optional ClinVar + WeasyPrint PDF
  --workflow genealogy -> genealogy stub path
```



## Performance benchmarks

**HG002 NIST Illumina 2x250 full-genome clinical workflow (June 2026)**

Evaluated with ClinVar overlay and clinical annotation stack against GIAB v4 truth set:

| Variant Type | Recall   | Precision | F1 Score |
|---|---|---|---|
| SNP | **99.40%** | **99.90%** | 0.9965 |
| INDEL | **98.13%** | **99.46%** | 0.9879 |

Output locations:
- Annotated VCF: `results/benchmark_clinvar_dv/clinical/annotations/HG002_NIST_Illumina2x250_R.clinical.annotated.vcf.gz`
- hap.py metrics: `results/benchmark_clinvar_dv/bench_eval/hg002_dv_bench.summary.csv`
- Clinical report: `results/benchmark_clinvar_dv/clinical/reports/`

To verify this performance on your system:

```bash
bash scripts/prepare_hg002_reads.sh full   # Download full HG002 dataset
bash scripts/prepare_clinvar_resource.sh   # get ClinVar GRCh38 (auto-detected)
nextflow run . \
  -profile docker,adaptive_local \
  --max_cpus 30 \
  --reads data/benchmark/giab_hg002/HG002_NIST_Illumina2x250_R{1,2}.fastq.gz \
  --reference data/GRCh38.fasta \
  --workflow clinical \
  --caller deepvariant \
  --clinvar_vcf data/clinvar/clinvar_GRCh38_chr.vcf.gz \
  --outdir results/my_benchmark_verify
```

**Note:** These metrics do not constitute clinical validation. They represent end-to-end pipeline performance on a single high-confidence sample. Clinical deployment requires regulatory validation, quality assurance, and environment-specific testing.

## Decision log

- Kept GIAB HG002 + exome-constrained evaluation as the primary clinical benchmark context.
- Implemented and retained VQSR (trained filtering) for HaplotypeCaller, not just hard filters.
- Added/kept hard-filter path only as a baseline comparator.
- Stabilized clinical report endpoint by:
  - moving report generation into `scripts/build_clinical_report.py`
  - staging script as an explicit Nextflow input
  - using `python:3.11` for task runtime compatibility
  - pinning `weasyprint==61.2` with `pydyf==0.10.0`
- Recorded comparison outputs in dedicated artifacts under `results/benchmark_vqsr/` for reproducibility.

## Important assumptions

- Input data are paired-end FASTQs matching the default `*_R{1,2}.fastq.gz` layout.
- The pipeline is boilerplate-grade and needs environment-specific hardening before regulated clinical deployment.
- **Reference genome is GRCh38** (hg38). All resource files (VQSR, ClinVar, phenotype mappings) must be for this build.
- The clinical branch assumes the source VCF already carries consequence annotations if no ClinVar resource is supplied.
- **ClinVar overlay** is optional; use `--clinvar_vcf` to add clinical significance labels (pathogenicity, phenotype). Pre-built GRCh38 copies exist in `data/clinvar/`. To refresh from NCBI, run `bash scripts/prepare_clinvar_resource.sh` (validates reference match).
- The genealogy branch is a production routing stub for Eagle2 and Beagle 5.4 rather than a fully wired reference-panel execution path.
- The `test` profile is for CI/stub smoke validation and not for biological correctness.
- GATK VQSR is data-hungry and is generally a poor fit for single exomes or very small exome cohorts; hard filtering or DeepVariant will often remain the more defensible option in that setting.

## Single-sample recommendation

- For single-sample exome/panel runs, prefer `--caller deepvariant`.
- Reserve HaplotypeCaller+VQSR for larger cohorts where recalibration modeling is well-supported.

## Container profiles

- `docker`
- `singularity`
- `apptainer`

Use `nextflow config` to inspect the resolved profile configuration for your target platform.

## Operational Runbooks

### Pre-Run Validation

Always run the preflight checklist before executing a full pipeline run, especially for clinical or long-duration benchmarks:

```bash
# Validate all prerequisites (reference files, disk space, memory, container availability)
bash scripts/preflight_checks.sh

# Validate parameters only without executing the pipeline
nextflow run . --validate_only true -profile docker
```

### Clinical Runs

For patient variant calls with clinical reporting:

```bash
# 1. Ensure ClinVar is current
bash scripts/prepare_clinvar_resource.sh

# 2. Validate input FASTQ pairing and sample IDs
bash scripts/validate_inputs.sh --fastq-pattern "data/*_R{1,2}.fastq.gz"

# 3. Run clinical pipeline with DeepVariant (recommended for single samples)
nextflow run . \
  -profile docker,adaptive_local \
  --workflow clinical \
  --caller deepvariant \
  --clinvar_vcf data/clinvar/clinvar_GRCh38_chr.vcf.gz \
  -resume

# 4. Review QC summary and clinical report
cat results/qc/multiqc_report.html  # FastQC + alignment summary
cat results/clinical_report.txt      # Variant report with ClinVar annotations
```

### Genealogy Runs

For pedigree-aware analysis:

```bash
# Run with HaplotypeCaller (recommended for cohorts)
nextflow run . \
  -profile docker,adaptive_local \
  --workflow genealogy \
  --caller haplotypecaller \
  --skip_qc false \
  -resume
```

### Benchmark Runs

For performance validation and accuracy comparison:

```bash
# Download and run HG002 benchmark
bash scripts/run_hg002_benchmark.sh full \
  -profile docker,adaptive_local \
  --max_cpus 30 \
  --deepvariant_cpus 30 \
  --deepvariant_memory '96 GB' \
  -resume

# Compare with truth using hap.py (if installed)
hap.py results/clinical/deepvariant.vcf.gz \
  <path-to-truth-vcf> \
  -f <path-to-confident-regions-bed> \
  -o results/benchmark_happycomparison
```

### Resuming After Failures

The pipeline integrates with Nextflow's `-resume` mechanis for efficient recovery:

```bash
# After fixing an error, re-run with cache preservation
nextflow run . -profile docker,adaptive_local -resume

# View execution timeline and identify failed tasks
nextflow run . -profile docker,adaptive_local -with-timeline results/timeline.html
```

### Monitoring Long Runs

For background execution:

```bash
# Run in the background with logging
nohup nextflow run . -profile docker,adaptive_local > nextflow.log 2>&1 &

# Monitor progress in a separate terminal
tail -f nextflow.log
nextflow log                              # Show all runs
nextflow log <runName> -f "{name} {exit}" # Show task statuses
```

## Validated Release Profiles

These profiles have undergone quality assurance for specific use cases:

### Profile: `docker,adaptive_local` (Default)

- **Use case**: Single-machine production runs with dynamic resource scaling.
- **When to use**: Most clinical and genealogy deployments.
- **Key settings**: CPU capped via `--max_cpus`, memory auto-scaled per process.
- **Validation**: Tested against HG002 (all callers), exome subsets, full genome runs.

**Validated command**:

```bash
nextflow run . -profile docker,adaptive_local \
  --max_cpus 30 --deepvariant_cpus 30 --deepvariant_memory '96 GB' \
  --workflow clinical --caller deepvariant
```

### Profile: `docker,aggressive_local`

- **Use case**: High-throughput analysis on memory-rich, single-tenant machines.
- **When to use**: Benchmark runs, cohort processing on dedicated hardware.
- **Key settings**: Aggressive CPU/RAM per process, minimal queuing.
- **Validation**: Tested for throughput on 64+ core systems; not recommended for shared compute.

**Validated command**:

```bash
nextflow run . -profile docker,aggressive_local \
  --workflow clinical --caller deepvariant \
  --clinvar_vcf data/clinvar/clinvar_GRCh38_chr.vcf.gz
```

### Profile: `test` (CI Only)

- **Use case**: Rapid stub-run validation in continuous integration.
- **When to use**: Pull request smoke tests, local development iteration.
- **Key settings**: Stub mode, minimal data, no actual compute.
- **Validation**: All routes (clinical/genealogy, all callers) complete in <2 min.

**Validated command**:

```bash
nextflow run . -profile test -stub-run \
  --workflow clinical --caller deepvariant
```

## Release & Version Management

### Version Bumping

Follow semantic versioning. Update the following before release:

1. [nextflow.config](nextflow.config) - update `manifest.version`
2. [CHANGELOG.md](CHANGELOG.md) - record all changes
3. [.github/workflows/release.yml](.github/workflows/release.yml) - optional: auto-tag on version bump

```bash
# Example: bumping to v0.2.0
sed -i "s/version = .*/version = '0.2.0'/" nextflow.config
git add CHANGELOG.md nextflow.config
git commit -m "chore: release v0.2.0"
git tag v0.2.0
git push origin main --tags
```

### Container Image Pinning

When releasing, pin container image digests for reproducibility:

```bash
# In nextflow.config, update images from tag to digest
# BEFORE:  container 'staphb/bcftools:latest'
# AFTER:   container 'staphb/bcftools@sha256:abc123...'

# Pull and inspect current digest
docker pull staphb/bcftools:latest
docker inspect staphb/bcftools:latest | grep -i digest
```

See [RELEASE_CHECKLIST.md](RELEASE_CHECKLIST.md) for complete pre-release validation steps.


