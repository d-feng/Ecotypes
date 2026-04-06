# STAD Mutation Integration Module

This module adds mutation-aware analysis on top of the existing STAD EcoTyper recovery run in:

- `RecoveryOutput_STAD/bulk_stad_data`

The workflow keeps ecotypes fixed, then adds mutation-derived biomarkers and within-ecotype mutation clustering.

Run all commands from the repository root.

## Scripts

- `mutation_integration/build_stad_mutation_features.R`
  - Reads a TCGA-STAD MAF
  - Harmonizes TCGA barcodes to the dotted sample IDs used by the STAD EcoTyper run
  - Builds gene-level, pathway-level, and burden-style mutation features
  - Merges them with dominant EcoTyper ecotypes and STAD metadata

- `mutation_integration/run_stad_ecotype_mutation_association.R`
  - Tests ecotype-vs-rest mutation enrichment
  - Uses Fisher tests for binary features and Wilcoxon tests for continuous features

- `mutation_integration/run_stad_within_ecotype_clustering.R`
  - Performs within-ecotype mutation clustering for large ecotypes using Gower distance and silhouette-selected hierarchical clustering

- `mutation_integration/make_stad_mutation_figures.R`
  - Generates a binary-feature ecotype heatmap and mutation burden boxplot

## Expected input

The main external input is a TCGA-STAD MAF with at least:

- `Tumor_Sample_Barcode`
- `Hugo_Symbol`
- `Variant_Classification`

Default expected path:

- `data/tcga_stad_raw/TCGA.STAD.open_masked_somatic.maf.tsv.gz`

## Example run order

```r
Rscript mutation_integration/build_stad_mutation_features.R --maf_path data/tcga_stad_raw/TCGA.STAD.open_masked_somatic.maf.tsv.gz
Rscript mutation_integration/run_stad_ecotype_mutation_association.R
Rscript mutation_integration/run_stad_within_ecotype_clustering.R
Rscript mutation_integration/make_stad_mutation_figures.R
Rscript mutation_integration/make_stad_gastric_manuscript_panels.R
```

## Outputs

Written under:

- `MutationIntegration_STAD`

Key files:

- `tables/stad_mutation_features.tsv`
- `tables/sample_ecotype_mutation_table.tsv`
- `tables/ecotype_mutation_associations.tsv`
- `tables/ecotype_mutation_associations_top.tsv`
- `tables/within_ecotype/*.tsv`
- `plots/ecotype_mutation_heatmap.png`
- `plots/ecotype_mutation_burden_boxplot.png`

## Current scope

Implemented now:

- recurrent gastric driver features
- pathway-level mutation biomarkers
- mutation burden features
- ecotype association testing
- within-ecotype mutation subclustering

Planned next:

- mutational signatures
- MSI/HRD/aneuploidy integration
- oncoprint generation
- pathway- and paper-specific comparison panels

## Upstream dependency

This module assumes the STAD EcoTyper recovery has already been run and that these files exist:

- `RecoveryOutput_STAD/bulk_stad_data/Ecotypes/ecotype_assignment.txt`
- `RecoveryOutput_STAD/bulk_stad_data/Ecotypes/ecotype_abundance.txt`

If they do not exist yet, start with the recovery step described in `README_COLLABORATOR_HANDOFF.md`.
