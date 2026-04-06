# publication0_10.1016_j.ccell.2026.01.015

This module ports the **signature-building logic** from the authors' `nICT_GC` repository for:

- EC1-EC5 ecotype signatures from Figure 4
- MOS subtype signatures from Figure 8

It is designed to be usable without the authors' private `.rds` objects by rebuilding signatures from standard tabular inputs.

## Files

- `publication0_utils.R`
  - shared parsing and signature helper functions

- `publication0_10.1016_j.ccell.2026.01.015_ecotype_signature_logic.R`
  - ports Figure 4 logic
  - expects a marker table with columns:
    - `gene`
    - `cluster`
    - `p_val_adj`
  - creates EC-style signatures and optional GSVA scores on bulk expression

- `publication0_10.1016_j.ccell.2026.01.015_mos_signature_logic.R`
  - ports Figure 8 logic
  - expects:
    - an expression matrix
    - sample metadata with a group column
  - builds one-vs-rest DE signatures intersected with highest-expression genes
  - scores the resulting signatures with GSVA

## DOI

- `10.1016/j.ccell.2026.01.015`

## Example usage

```r
Rscript publication0_10.1016_j.ccell.2026.01.015/publication0_10.1016_j.ccell.2026.01.015_ecotype_signature_logic.R ^
  --marker_table_path C:/path/TumorCombine_Filter_Marker.txt ^
  --expression_path C:/path/expression.tsv ^
  --metadata_path C:/path/metadata.tsv
```

```r
Rscript publication0_10.1016_j.ccell.2026.01.015/publication0_10.1016_j.ccell.2026.01.015_mos_signature_logic.R ^
  --expression_path C:/path/expression.tsv ^
  --metadata_path C:/path/sample_metadata.tsv ^
  --group_column MOS
```

## Notes

- The authors' GitHub contains the **logic** for deriving signatures, but not always the final printed gene lists.
- Figure 4 uses marker-derived EC gene sets.
- Figure 8 derives MOS signatures as the intersection of:
  - one-vs-rest significant genes
  - genes with highest mean expression in that subtype
