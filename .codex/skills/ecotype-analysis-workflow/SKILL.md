---
name: ecotype-analysis-workflow
description: Run or extend this repository's ecotype analysis workflow for prototype ecotype discovery, EcoTyper bulk recovery, STAD mutation integration, and publication0 gastric scoring. Use when Codex needs to execute, document, debug, adapt, or compare any of these connected modules in this repo.
---

# Ecotype Analysis Workflow

Use this skill when working inside this repository on any of the connected ecotype analysis layers.

## Quick route

Choose the smallest path that fits the task:

- Prototype ecotype logic only:
  - use `phase1_ecotype_pipeline/`
- Real cohort ecotype recovery:
  - use the repository root EcoTyper scripts
- Gastric STAD mutation work:
  - use `mutation_integration/`
- Gastric publication-aligned scoring:
  - use `publication0_10.1016_j.ccell.2026.01.015/`

## Working rules

- Run commands from the repository root unless a script explicitly requires another working directory.
- Prefer repo-relative paths over machine-specific absolute paths.
- Treat generated cohort data and output folders as reproducible artifacts, not guaranteed git-tracked inputs.
- Before extending an analysis, check whether the required upstream output already exists.

## Workflow order

Follow this dependency order unless the user asks for only one layer:

1. If the task is conceptual, lightweight, or prototype-only, start with `phase1_ecotype_pipeline/`.
2. For real cohort ecotypes, run or inspect `EcoTyper_recovery_bulk.R` outputs first.
3. For STAD mutation analyses, require existing ecotype assignments in `RecoveryOutput_STAD/bulk_stad_data/Ecotypes/`.
4. For `publication0` panels, require bulk expression plus CE assignments from the STAD recovery.

Read `references/run-map.md` when you need the concrete command order and expected files.

## Key repo anchors

- `README_COLLABORATOR_HANDOFF.md`
- `phase1_ecotype_pipeline/RELATIONSHIP_TO_ECOTYPER.md`
- `mutation_integration/README_STAD_mutation_module.md`
- `publication0_10.1016_j.ccell.2026.01.015/README.md`

## Output expectations

When asked to run or modify the workflow:

- state which layer you are operating in
- name the upstream files you depend on
- report the output folder that was created or updated
- call out if results depend on untracked local data or generated outputs
