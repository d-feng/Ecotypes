source("mutation_integration/mutation_utils.R")

defaults <- list(
  sample_table_path = "C:/Users/difen/Rcode/ecotyper/MutationIntegration_STAD/tables/sample_ecotype_mutation_table.tsv",
  output_dir = "C:/Users/difen/Rcode/ecotyper/MutationIntegration_STAD"
)
args <- parse_simple_args(defaults)

dir_create_safe(args$output_dir)
dir_create_safe(file.path(args$output_dir, "tables"))

dat <- fread(args$sample_table_path)
dat <- dat[!is.na(Ecotype)]

feature_cols <- setdiff(
  names(dat),
  c("ID", "Ecotype", "Tissue", "Histology", "Type", "Stage", "PrimaryDiagnosis", "Gender",
    "AgeYears", "VitalStatus", "OS_Time", "OS_Status")
)

binary_features <- feature_cols[vapply(dat[, ..feature_cols], function(x) all(x %in% c(0, 1)), logical(1))]
continuous_features <- setdiff(feature_cols, binary_features)

results <- rbindlist(lapply(sort(unique(dat$Ecotype)), function(ecotype) {
  in_idx <- dat$Ecotype == ecotype
  out_idx <- !in_idx

  binary_res <- rbindlist(lapply(binary_features, function(feature) {
    stats <- safe_fisher(dat[[feature]][in_idx], dat[[feature]][out_idx])
    data.table(
      Ecotype = ecotype,
      Feature = feature,
      FeatureType = "binary",
      PValue = stats$p_value,
      EffectSize = stats$odds_ratio,
      MeanInEcotype = mean(dat[[feature]][in_idx], na.rm = TRUE),
      MeanOutsideEcotype = mean(dat[[feature]][out_idx], na.rm = TRUE)
    )
  }))

  continuous_res <- rbindlist(lapply(continuous_features, function(feature) {
    stats <- safe_wilcox(dat[[feature]][in_idx], dat[[feature]][out_idx])
    data.table(
      Ecotype = ecotype,
      Feature = feature,
      FeatureType = "continuous",
      PValue = stats$p_value,
      EffectSize = stats$effect_size,
      MeanInEcotype = mean(dat[[feature]][in_idx], na.rm = TRUE),
      MeanOutsideEcotype = mean(dat[[feature]][out_idx], na.rm = TRUE)
    )
  }))

  rbind(binary_res, continuous_res, fill = TRUE)
}))

results[, FDR := p.adjust(PValue, method = "BH"), by = Ecotype]
results[, Direction := ifelse(EffectSize >= 0, "enriched", "depleted")]
results[, AbsEffectSize := abs(EffectSize)]
setorder(results, Ecotype, FDR, -AbsEffectSize)

top_hits <- results[FDR <= 0.1]
if (nrow(top_hits) == 0) {
  top_hits <- results[, .SD[1:min(.N, 10)], by = Ecotype]
}

fwrite(results, file.path(args$output_dir, "tables", "ecotype_mutation_associations.tsv"), sep = "\t")
fwrite(top_hits, file.path(args$output_dir, "tables", "ecotype_mutation_associations_top.tsv"), sep = "\t")
