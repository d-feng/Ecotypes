## Relationship to the Main `ecotyper` Workflow

This folder preserves the earlier Phase 1 ecotype-only scaffold in the same repository as the later gastric cancer work.

### What this scaffold is

- A lightweight prototype for ecotype discovery from:
  - bulk expression
  - sample metadata
  - a sample-by-cell-state abundance matrix
- A dependency-light implementation of the ecotype logic we used before switching to the full EcoTyper recovery workflow
- Useful for local testing, method development, and understanding the ecotype abstraction without needing the whole upstream EcoTyper stack

### How it differs from the main workflow

- `phase1_ecotype_pipeline/`
  - prototype / scaffold
  - expects precomputed cell-state abundances
  - performs ecotype discovery and sample assignment locally
  - mirrors EcoTyper-style output files, but does not run CIBERSORTx or EcoTyper state discovery

- repository root workflow
  - real EcoTyper-based implementation
  - uses the upstream EcoTyper scripts such as `EcoTyper_recovery_bulk.R`
  - supports the STAD analyses added in:
    - `mutation_integration/`
    - `publication0_10.1016_j.ccell.2026.01.015/`

### Recommended usage

- Use `phase1_ecotype_pipeline/` when:
  - you want a small, inspectable prototype
  - you already have cell-state abundances
  - you want to test ecotype logic on synthetic or prepared matrices

- Use the main repository workflow when:
  - you want actual EcoTyper recovery outputs
  - you want the STAD gastric analyses
  - you want mutation integration or publication-aligned scoring panels

### Practical lineage

The Phase 1 scaffold came first. The later work in this repository is the production path that replaced it for the STAD analyses, but the scaffold is kept here because it still documents the design logic and remains useful for rapid prototyping.
