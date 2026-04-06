source("mutation_integration/mutation_utils.R")

defaults <- list(
  maf_path = "C:/Users/difen/Rcode/ecotyper/data/tcga_stad_raw/TCGA.STAD.mutect2_somatic.maf.gz",
  metadata_path = "C:/Users/difen/Rcode/ecotyper/data/tcga_stad_prepared/bulk_stad_annotation.txt",
  ecotype_abundance_path = "C:/Users/difen/Rcode/ecotyper/RecoveryOutput_STAD/bulk_stad_data/Ecotypes/ecotype_abundance.txt",
  output_dir = "C:/Users/difen/Rcode/ecotyper/MutationIntegration_STAD"
)
args <- parse_simple_args(defaults)

dir_create_safe(args$output_dir)
dir_create_safe(file.path(args$output_dir, "tables"))

if (!file.exists(args$maf_path)) {
  stop(sprintf("MAF file not found: %s", args$maf_path))
}

metadata <- fread(args$metadata_path)
metadata[, ID := harmonize_tcga_sample_id(ID)]
metadata <- unique(metadata, by = "ID")

ecotype_table <- dominant_ecotype_table(args$ecotype_abundance_path)
sample_ids <- intersect(metadata$ID, ecotype_table$ID)

maf <- read_maf_table(args$maf_path)
maf <- maf[sample_id %in% sample_ids]

gene_list <- get_stad_driver_genes()
pathway_map <- get_stad_pathways()

binary_features <- build_binary_feature_table(maf, sample_ids, gene_list, pathway_map)
continuous_features <- build_continuous_feature_table(maf, sample_ids)

mutation_features <- merge(binary_features, continuous_features, by = "ID", all = TRUE)
sample_table <- merge(ecotype_table, metadata, by = "ID", all.x = TRUE)
sample_table <- merge(sample_table, mutation_features, by = "ID", all.x = TRUE)

feature_cols <- setdiff(names(mutation_features), "ID")
for (col in feature_cols) {
  if (is.numeric(sample_table[[col]])) {
    set(sample_table, which(is.na(sample_table[[col]])), col, 0)
  }
}

fwrite(mutation_features, file.path(args$output_dir, "tables", "stad_mutation_features.tsv"), sep = "\t")
fwrite(sample_table, file.path(args$output_dir, "tables", "sample_ecotype_mutation_table.tsv"), sep = "\t")
fwrite(as.data.table(maf), file.path(args$output_dir, "tables", "stad_maf_filtered.tsv.gz"), sep = "\t")

summary_lines <- c(
  sprintf("Samples with ecotype labels: %d", length(sample_ids)),
  sprintf("Samples with >=1 nonsilent mutation: %d", sum(mutation_features$nonsilent_mutation_count > 0)),
  sprintf("Median nonsilent mutation count: %.2f", median(mutation_features$nonsilent_mutation_count)),
  sprintf("Features generated: %d", ncol(mutation_features) - 1),
  sprintf("Source MAF: %s", args$maf_path)
)
writeLines(summary_lines, file.path(args$output_dir, "tables", "build_summary.txt"))
