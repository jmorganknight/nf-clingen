# Changelog

All notable changes to nf-clingen are documented in this file. Follow [Semantic Versioning](https://semver.org/).

## [0.2.0] - 2025-06-09

### Added

- **YAML Parameter Support**: Optional `params.yaml` configuration file for runnable defaults and parameter documentation
  - Known-key validation: unknown YAML keys trigger warnings
  - Three preset configurations (clinical, genealogy, benchmark) for common workflows
  - Full parameter catalog with supported values, examples, and operational notes
  
- **Expanded Nextflow Schema**: All 24 operational parameters now documented
  - `scratch_dir`: Route work files to fast scratch storage
  - `patient_phenotype`, `phenotype_gene_map`: Phenotype-aware variant filtering
  - `max_cpus`, `max_ram_gb`, `max_forks`: Resource controls for local execution
  - `deepvariant_cpus`, `deepvariant_memory`, `deepvariant_time`: DeepVariant tuning
  
- **Comprehensive Runbooks in README**: Operational guidance for clinical, genealogy, and benchmark workflows
  - Pre-run validation and parameter checking
  - Resume and failure recovery procedures
  - Background execution and monitoring for long-running jobs
  
- **Validated Release Profiles**: Three profiles certified for production use
  - `docker,adaptive_local` (default): Single-machine clinical runs with dynamic resource scaling
  - `docker,aggressive_local`: High-throughput analysis on memory-rich hosts
  - `test`: CI validation with stub-run mode
  
- **Enhanced CI/CD Matrix**: Expanded GitHub Actions smoke testing
  - Schema validation job (lint nextflow_schema.json)
  - 6-route matrix: clinical/genealogy × haplotypecaller/deepvariant × minimap2/bwamem2
  - Parameter validation and config validation checks
  
- **Preflight Validation Script** (`scripts/preflight_checks.sh`):
  - Checks Java, Nextflow, Docker/container runtime
  - Validates reference files, input FASTQs, ClinVar resources
  - Confirms disk space (>50 GB) and memory availability
  - Suggests ready-to-run example commands
  
- **Release Checklist** (`RELEASE_CHECKLIST.md`):
  - Pre-release code review, testing, and documentation steps
  - Version bumping guidance (semantic versioning)
  - Container image digest pinning for reproducibility
  - Post-release verification and rollback procedures

### Changed

- **README Reorganization**: Added operational sections for clarity
  - New: Operational Runbooks (pre-run, clinical, genealogy, benchmark, resume, monitoring)
  - New: Validated Release Profiles (with tested command examples)
  - New: Release & Version Management (version bumping, container pinning)

- **nextflow.config**: Added YAML parameter loader
  - Optional auto-loading of params.yaml if present
  - Known-key allowlist (24 supported parameters)
  - Unknown keys trigger stderr warnings
  - Parser supports quoted strings, type inference, inline comments

### Fixed

- Parameter validation now covers all 24 operational parameters (was 9)
- Help text (`nextflow run . --help`) now complete with all resource controls and tuning parameters

### Known Issues

- `genealogy` branch is a routing stub; Eagle2 and Beagle orchestration not fully wired
- GATK VQSR requires large cohorts for accurate model training; single exomes should use DeepVariant or hard filtering

### Migration Notes

No breaking changes. Existing command-line workflows continue unchanged.

To adopt YAML configuration:
1. Copy `params.yaml` and customize parameter values
2. Run: `nextflow run . -profile docker` (auto-loads params.yaml)
3. CLI flags override YAML values

---

## [0.1.0] - 2025-05-15

### Initial Release

- Basic nf-clingen functionality: clinical/genealogy routing with haplotypecaller/deepvariant support
- Local Docker-based execution with process modules for QC, alignment, calling, annotation
- Minimal smoke-test CI with 3 routes (clinical-hc, genealogy-hc, clinical-dv)
- Decision log documenting architecture and trade-offs

---

## Release Guidelines

### Versioning Strategy

- **MAJOR**: Breaking changes, schema/config reorganization, backward-incompatible refactoring
- **MINOR**: New features, enhanced functionality, non-breaking improvements
- **PATCH**: Bug fixes, documentation updates, minor optimizations

### Release Process

1. Complete pre-release checklist (see [RELEASE_CHECKLIST.md](RELEASE_CHECKLIST.md))
2. Update version in `nextflow.config` (manifest.version)
3. Update this file with new section at top
4. Commit: `git commit -m "chore: release vX.Y.Z"`
5. Tag: `git tag -a vX.Y.Z -m "Release vX.Y.Z - [brief summary]"`
6. Push: `git push origin main --tags`

### Testing Before Release

- [ ] All CI checks pass (`nextflow run . -profile test -stub-run` for all 6 routes)
- [ ] Full HG002 benchmark completes without regressions
- [ ] Parameter validation passes: `nextflow run . --validate_only true`
- [ ] Help rendering is correct: `nextflow run . --help true | head -30`
- [ ] Fresh clone test on a clean environment

---

**Latest Release**: [v0.2.0](https://github.com/[owner]/nf-clingen/releases/tag/v0.2.0)  
**Maintained by**: GitHub Copilot  
**Next Planned Release**: v0.3.0 (Q3 2025) - Reproducibility locks, clinical validation profiles, audit logging
