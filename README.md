# nf-prism

nf-prism is a modular Nextflow DSL2 router for paired-end human sequencing runs that need to switch between clinical reporting and genealogy endpoints while keeping the upstream stack containerized and low on unnecessary intermediate disk I/O.

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

Clinical run with DeepVariant and ClinVar overlay:

```bash
nextflow run . \
  -profile docker \
  --workflow clinical \
  --caller deepvariant \
  --clinvar_vcf /path/to/clinvar.vcf.gz
```

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

Smoke-test the full graph with synthetic test inputs:

```bash
nextflow run . -profile test -stub-run --workflow clinical --caller haplotypecaller
```

Additional smoke routes:

```bash
nextflow run . -profile test -stub-run --workflow genealogy --caller haplotypecaller
nextflow run . -profile test -stub-run --workflow clinical --caller deepvariant
```

## Important assumptions

- Input data are paired-end FASTQs matching the default `*_R{1,2}.fastq.gz` layout.
- The pipeline is boilerplate-grade and needs environment-specific hardening before regulated clinical deployment.
- The clinical branch assumes the source VCF already carries consequence annotations if no ClinVar resource is supplied.
- The genealogy branch is a production routing stub for Eagle2 and Beagle 5.4 rather than a fully wired reference-panel execution path.
- The `test` profile is for CI/stub smoke validation and not for biological correctness.

## Container profiles

- `docker`
- `singularity`
- `apptainer`

Use `nextflow config` to inspect the resolved profile configuration for your target platform.