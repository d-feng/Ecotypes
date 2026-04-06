# Phase 1 Ecotype-Only Pipeline

This prototype is preserved inside the main `ecotyper` repository so the original Phase 1 ecotype-only work stays together with the later EcoTyper, STAD mutation, and `publication0` modules. See `RELATIONSHIP_TO_ECOTYPER.md` for how this scaffold fits into the larger workflow.

This project implements a lightweight, R-based Phase 1 pipeline for ecotype discovery from tumor expression data, sample metadata, and a sample-by-cell-state abundance matrix.

The refactored logic follows the EcoTyper ecotype concept more closely:

- start from cell states that have already been defined or recovered
- mark the dominant state per cell type in each sample
- compute state-state Jaccard overlap across samples
- zero out non-significant state overlaps using a hypergeometric test
- hierarchically cluster states into ecotypes
- choose the ecotype number with maximum average silhouette
- assign each sample to the ecotype with the largest aggregate state abundance

This implementation is dependency-light so it can run locally. It does not reproduce EcoTyper's upstream CIBERSORTx and NMF discovery steps, but it is designed to accept their state abundance outputs directly.

## Project layout

- `config/example_config.yml`: runnable example configuration
- `R/run_phase1_pipeline.R`: main driver script
- `R/generate_synthetic_data.R`: creates a synthetic example dataset
- `data/`: input data location
- `output/`: pipeline results

## Input expectations

### Expression matrix

- tab-delimited text file
- rows are genes
- columns are samples
- first column contains gene identifiers

### Metadata table

- tab-delimited text file
- must contain `SampleID`

### Cell-state abundance matrix

- tab-delimited text file
- rows are cell states
- columns are samples
- first column contains state identifiers
- recommended state naming format: `CellType|StateLabel`

## Run the example

Generate synthetic data:

```powershell
Rscript R/generate_synthetic_data.R
```

Run the pipeline:

```powershell
Rscript R/run_phase1_pipeline.R --config config/example_config.yml
```

## Main outputs

- `output/tables/state_to_ecotype.tsv`
- `output/tables/ecotype_assignments.tsv`
- `output/tables/ecotype_abundance.tsv`
- `output/tables/ecotype_number_metrics.tsv`
- `output/tables/top_markers_by_ecotype.tsv`
- `output/tables/metadata_associations.tsv`
- `output/plots/jaccard_state_heatmap.png`
- `output/plots/ecotype_number_selection.png`
- `output/plots/state_abundance_heatmap.png`
- `output/plots/sample_pca_from_states.png`
- `output/plots/ecotype_marker_heatmap.png`

## EcoTyper-style mirrored outputs

The scaffold also writes a recovery-like folder layout so downstream code can consume more familiar file names:

- `output/<CellType>/state_abundances.txt`
- `output/<CellType>/state_assignment.txt`
- `output/Ecotypes/ecotype_assignment.txt`
- `output/Ecotypes/ecotype_abundance.txt`

These are mirrors of the pipeline's internally generated results and are intended to resemble EcoTyper recovery outputs.

## How this maps to EcoTyper

If you run EcoTyper discovery or recovery externally, the clean integration point here is the state abundance matrix:

1. Run EcoTyper upstream to define or recover cell states.
2. Export the state abundance matrix per sample.
3. Save it as `cell_state_abundances.tsv`.
4. Point `cell_state_abundances_path` to that file in the config.
5. Re-run this pipeline to discover ecotypes, assign samples, and generate figures.
