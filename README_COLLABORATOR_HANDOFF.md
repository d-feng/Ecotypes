# Collaborator Handoff

This repository now contains three related layers of work:

1. The original EcoTyper codebase at the repository root
2. A lightweight prototype ecotype scaffold in `phase1_ecotype_pipeline/`
3. STAD-specific downstream analyses in:
   - `mutation_integration/`
   - `publication0_10.1016_j.ccell.2026.01.015/`

There is also a repo-local Codex skill for this workflow:

- `.codex/skills/ecotype-analysis-workflow/`

This document is the fastest way for a collaborator to get oriented in a fresh Codex session.

## What is already in git

Committed code on branch `codex-stad-publication0-signatures` includes:

- Windows-compatible updates needed to run `EcoTyper_recovery_bulk.R`
- the integrated Phase 1 scaffold in `phase1_ecotype_pipeline/`
- the STAD mutation module in `mutation_integration/`
- the gastric publication-inspired module in `publication0_10.1016_j.ccell.2026.01.015/`

Large generated outputs and downloaded cohort data were intentionally not committed.

## Recommended starting point

Use the repository root for real EcoTyper-based work.

- Run actual recovery with:
  - `EcoTyper_recovery_bulk.R`
- Use `phase1_ecotype_pipeline/` only when:
  - you want a lightweight prototype
  - you already have a sample-by-cell-state abundance matrix
  - you want to test ecotype logic without the full upstream stack

## Environment assumptions

Tested in a Windows Codex environment with PowerShell and Rscript.

Required R packages across the added modules include:

- `data.table`
- `yaml`
- `ggplot2`
- `ComplexHeatmap`
- `RColorBrewer`
- `cluster`
- `AnnotationDbi`
- `org.Hs.eg.db`
- `GSVA`

The root EcoTyper README lists the broader package set needed by the original pipeline.

## Important path rule

All new scripts should be run from the repository root:

```powershell
cd <repo-root>
```

Then call scripts with repo-relative paths such as:

```powershell
Rscript mutation_integration/build_stad_mutation_features.R --maf_path data/tcga_stad_raw/TCGA.STAD.open_masked_somatic.maf.tsv.gz
```

Do not rely on `C:/Users/difen/...` paths on another machine.

## Minimal run map

### 1. Prototype ecotype workflow

See:

- `phase1_ecotype_pipeline/README.md`
- `phase1_ecotype_pipeline/RELATIONSHIP_TO_ECOTYPER.md`

Main driver:

- `phase1_ecotype_pipeline/R/run_phase1_pipeline.R`

### 2. Real EcoTyper gastric recovery

Prepare bulk STAD input:

```powershell
Rscript prepare_tcga_stad_ecotyper_input.R
```

Then run bulk recovery from the repository root:

```powershell
Rscript EcoTyper_recovery_bulk.R --discovery Carcinoma --matrix data/tcga_stad_prepared/bulk_stad_data.txt --annotation data/tcga_stad_prepared/bulk_stad_annotation.txt --output RecoveryOutput_STAD
```

Main expected recovery output:

- `RecoveryOutput_STAD/bulk_stad_data/Ecotypes/ecotype_assignment.txt`
- `RecoveryOutput_STAD/bulk_stad_data/Ecotypes/ecotype_abundance.txt`

### 3. STAD mutation integration

See:

- `mutation_integration/README_STAD_mutation_module.md`

Typical order:

```powershell
Rscript mutation_integration/build_stad_mutation_features.R --maf_path data/tcga_stad_raw/TCGA.STAD.open_masked_somatic.maf.tsv.gz
Rscript mutation_integration/run_stad_ecotype_mutation_association.R
Rscript mutation_integration/run_stad_within_ecotype_clustering.R
Rscript mutation_integration/make_stad_mutation_figures.R
Rscript mutation_integration/make_stad_gastric_manuscript_panels.R
```

### 4. Gastric publication-aligned scoring

See:

- `publication0_10.1016_j.ccell.2026.01.015/README.md`

The added panel script expects:

- bulk expression matrix
- CE assignments from the STAD EcoTyper recovery

Main script:

- `publication0_10.1016_j.ccell.2026.01.015/publication0_10.1016_j.ccell.2026.01.015_make_ce_ecscore_and_pathway_panels.R`

## What is not portable by default

These items are not guaranteed to exist in a fresh clone unless regenerated or copied separately:

- `data/` cohort downloads
- `RecoveryOutput_STAD/`
- `RecoveryOutput_actual/`
- `MutationIntegration_STAD/`
- `comparison_outputs/`

If a collaborator needs those exact files, they should either:

- regenerate them from the scripts in this repo, or
- receive a separate data/results bundle outside git

## Suggested handoff order for a collaborator

1. Read this file
2. Read `phase1_ecotype_pipeline/RELATIONSHIP_TO_ECOTYPER.md`
3. Use the repo-local skill at `.codex/skills/ecotype-analysis-workflow/` if working from Codex
4. Run or inspect the prototype scaffold if they need the ecotype logic first
5. Use real EcoTyper recovery for cohort work
6. Layer mutation and publication-style scoring afterward

## Notes about the gastric publication module

`publication0_10.1016_j.ccell.2026.01.015/` ports the authors' signature-building logic, not a private exact reproduction of their internal objects.

That means:

- the logic is portable
- exact original marker tables are not bundled here unless separately provided
