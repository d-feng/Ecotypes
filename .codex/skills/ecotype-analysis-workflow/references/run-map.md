# Run Map

## 1. Prototype ecotype scaffold

Use when you already have or want to simulate a sample-by-cell-state abundance matrix.

Main files:

- `phase1_ecotype_pipeline/R/run_phase1_pipeline.R`
- `phase1_ecotype_pipeline/R/pipeline_utils.R`
- `phase1_ecotype_pipeline/config/example_config.yml`

Typical use:

```powershell
Rscript phase1_ecotype_pipeline/R/generate_synthetic_data.R
Rscript phase1_ecotype_pipeline/R/run_phase1_pipeline.R --config phase1_ecotype_pipeline/config/example_config.yml
```

## 2. Real STAD EcoTyper recovery

Prepare inputs:

```powershell
Rscript prepare_tcga_stad_ecotyper_input.R
```

Run recovery:

```powershell
Rscript EcoTyper_recovery_bulk.R --discovery Carcinoma --matrix data/tcga_stad_prepared/bulk_stad_data.txt --annotation data/tcga_stad_prepared/bulk_stad_annotation.txt --output RecoveryOutput_STAD
```

Required downstream files:

- `RecoveryOutput_STAD/bulk_stad_data/Ecotypes/ecotype_assignment.txt`
- `RecoveryOutput_STAD/bulk_stad_data/Ecotypes/ecotype_abundance.txt`

## 3. STAD mutation integration

Use only after STAD recovery exists.

Typical order:

```powershell
Rscript mutation_integration/build_stad_mutation_features.R --maf_path data/tcga_stad_raw/TCGA.STAD.open_masked_somatic.maf.tsv.gz
Rscript mutation_integration/run_stad_ecotype_mutation_association.R
Rscript mutation_integration/run_stad_within_ecotype_clustering.R
Rscript mutation_integration/make_stad_mutation_figures.R
Rscript mutation_integration/make_stad_gastric_manuscript_panels.R
```

Main output folder:

- `MutationIntegration_STAD/`

## 4. publication0 gastric scoring

Use after STAD recovery and when cohort expression is available.

Main files:

- `publication0_10.1016_j.ccell.2026.01.015/publication0_10.1016_j.ccell.2026.01.015_ecotype_signature_logic.R`
- `publication0_10.1016_j.ccell.2026.01.015/publication0_10.1016_j.ccell.2026.01.015_mos_signature_logic.R`
- `publication0_10.1016_j.ccell.2026.01.015/publication0_10.1016_j.ccell.2026.01.015_make_ce_ecscore_and_pathway_panels.R`

STAD panel output folder:

- `publication0_10.1016_j.ccell.2026.01.015/output_stad_ce_panels/`

## Practical cautions

- Many result folders are generated locally and are not guaranteed to be present in a fresh clone.
- If a task fails because input data is missing, say which upstream folder or file must be regenerated.
- On Windows, keep using the patched repository version of `EcoTyper_recovery_bulk.R` and `pipeline/lib/multithreading.R`.
