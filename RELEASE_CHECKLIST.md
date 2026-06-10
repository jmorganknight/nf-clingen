# Release Checklist for nf-clingen

Use this checklist before tagging a new release. All steps should be completed and verified.

## Pre-Release Code Review

- [ ] All pull requests merged to `main` have been reviewed and approved
- [ ] No `TODO`, `FIXME`, or `XXX` comments remain in critical paths (`main.nf`, `nextflow.config`, core modules)
- [ ] Breaking changes are documented in CHANGELOG.md with migration guidance
- [ ] Performance regressions have been investigated (compare benchmark results with previous release)

## Testing & Validation

- [ ] All CI checks pass (GitHub Actions: lint-schema, smoke matrix for all 6 routes)
  ```bash
  # Verify locally before pushing:
  nextflow run . -profile test -stub-run --workflow clinical --caller haplotypecaller --aligner minimap2
  nextflow run . -profile test -stub-run --workflow clinical --caller haplotypecaller --aligner bwamem2
  nextflow run . -profile test -stub-run --workflow clinical --caller deepvariant --aligner minimap2
  nextflow run . -profile test -stub-run --workflow clinical --caller deepvariant --aligner bwamem2
  nextflow run . -profile test -stub-run --workflow genealogy --caller haplotypecaller --aligner minimap2
  nextflow run . -profile test -stub-run --workflow genealogy --caller deepvariant --aligner minimap2
  ```

- [ ] Full HG002 benchmark completes without errors
  ```bash
  bash scripts/pre_full_run_checklist.sh
  bash scripts/run_hg002_benchmark.sh mini  # Quick validation
  # Store results for comparison
  cp results/benchmark_hg002_dv/happy/HG002_giab.summary.csv results/benchmark_hg002_dv/happy/HG002_giab.summary.$(date +%Y%m%d).csv
  ```

- [ ] Parameter validation passes
  ```bash
  nextflow run . --validate_only true -profile docker
  ```

- [ ] Help text renders correctly with all 24 parameters
  ```bash
  nextflow run . --help true -profile docker,adaptive_local 2>&1 | head -30
  ```

## Documentation Updates

- [ ] [README.md](README.md) is current with latest configuration examples
- [ ] [CHANGELOG.md](CHANGELOG.md) has entry for this release with:
  - Version number (semantic: MAJOR.MINOR.PATCH)
  - Release date
  - Summary of changes
  - Breaking changes (if any) with migration path
  - New features
  - Bug fixes
  - Known issues (if any)

- [ ] All parameter descriptions in `nextflow_schema.json` are clear and complete
- [ ] Container image versions are documented (inspect `nextflow.config` process containers)
- [ ] Runbook examples in README are tested and work as documented

## Version Bumping

- [ ] Update `manifest.version` in [nextflow.config](nextflow.config)
  ```bash
  # Example: upgrading from v0.1.0 to v0.2.0
  sed -i "s/version = '0.1.0'/version = '0.2.0'/" nextflow.config
  ```

- [ ] Update container image digests for production reproducibility
  ```bash
  # Replace all 'tag' image references with 'digest' format
  # Inspect current digests:
  docker pull staphb/bcftools:latest
  docker inspect staphb/bcftools:latest | grep -i digest
  
  # Update in nextflow.config:
  # BEFORE: container 'staphb/bcftools:latest'
  # AFTER:  container 'staphb/bcftools@sha256:xxxxxx...'
  ```

- [ ] Verify version number consistency
  ```bash
  grep -r "0.2.0" nextflow.config  # Should appear exactly once in manifest.version
  ```

## Git & GitHub

- [ ] Commit all changes with clear message
  ```bash
  git add CHANGELOG.md nextflow.config README.md params.yaml nextflow_schema.json
  git commit -m "chore: release v0.2.0 - Add feature X, fix bug Y"
  ```

- [ ] Create annotated git tag (standard practice)
  ```bash
  git tag -a v0.2.0 -m "Release v0.2.0 - [brief description]"
  ```

- [ ] Push commits and tags to origin
  ```bash
  git push origin main
  git push origin --tags
  ```

- [ ] Verify tag appears in GitHub: https://github.com/[owner]/nf-clingen/releases

## Post-Release Verification

- [ ] Clone a fresh copy and run smoke tests
  ```bash
  cd /tmp
  git clone https://github.com/[owner]/nf-clingen.git nf-clingen-test
  cd nf-clingen-test
  git checkout v0.2.0
  nextflow run . -profile test -stub-run --workflow clinical --caller deepvariant
  ```

- [ ] Update any internal documentation or deployment scripts that hardcode version numbers

- [ ] Announce release (if applicable):
  - GitHub releases page (auto-generated from git tag)
  - Documentation changelog
  - Internal team communication

## Rollback Plan

If critical issues are discovered after release:

```bash
# Create a hotfix branch from the release tag
git checkout -b hotfix/v0.2.1 v0.2.0

# Apply fix, test thoroughly, then:
git tag -a v0.2.1 -m "Hotfix: [description]"
git push origin hotfix/v0.2.1
git push origin v0.2.1
```

## Release Schedule & Frequency

- **Patch releases** (v0.2.1): Bug fixes, minor updates. As needed.
- **Minor releases** (v0.3.0): New features, non-breaking enhancements. Every 1-3 months.
- **Major releases** (v1.0.0): Large features, breaking changes. Planned milestones.

---

**Completed by**: _________________  
**Date**: _________________  
**Notes**:
