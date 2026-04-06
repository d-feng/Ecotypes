source("publication0_10.1016_j.ccell.2026.01.015/publication0_utils.R")

defaults <- list(
  marker_table_path = "C:/Users/difen/Rcode/ecotyper/author_repo_raw/TumorCombine_Filter_Marker.txt",
  expression_path = "C:/Users/difen/Rcode/ecotyper/data/tcga_stad_prepared/bulk_stad_data.txt",
  metadata_path = "C:/Users/difen/Rcode/ecotyper/data/tcga_stad_prepared/bulk_stad_annotation.txt",
  p_adj_cutoff = "0.05",
  output_dir = "C:/Users/difen/Rcode/ecotyper/publication0_10.1016_j.ccell.2026.01.015/output_ecotype_logic"
)
args <- publication0_parse_args(defaults)

publication0_dir_create(args$output_dir)
publication0_dir_create(file.path(args$output_dir, "tables"))

marker_table <- fread(args$marker_table_path)
required_cols <- c("gene", "cluster", "p_val_adj")
missing_cols <- setdiff(required_cols, colnames(marker_table))
if (length(missing_cols) > 0) {
  stop(sprintf("Marker table missing required columns: %s", paste(missing_cols, collapse = ", ")))
}

p_adj_cutoff <- as.numeric(args$p_adj_cutoff)
clusters <- sort(unique(marker_table$cluster))
clusters <- clusters[grepl("^EC[1-9][0-9]*$", clusters)]
if (length(clusters) == 0) {
  stop("No EC-style clusters found in marker table.")
}

signature_list <- lapply(clusters, function(cluster_name) {
  unique(marker_table$gene[marker_table$cluster == cluster_name & marker_table$p_val_adj < p_adj_cutoff])
})
names(signature_list) <- clusters
signature_list <- signature_list[lengths(signature_list) > 0]

signature_table <- rbindlist(lapply(names(signature_list), function(cluster_name) {
  data.table(
    Signature = cluster_name,
    Gene = signature_list[[cluster_name]]
  )
}))

fwrite(signature_table, file.path(args$output_dir, "tables", "publication0_ecotype_signature_table.tsv"), sep = "\t")
saveRDS(signature_list, file.path(args$output_dir, "tables", "publication0_ecotype_signatures.rds"))

if (file.exists(args$expression_path)) {
  expr <- publication0_read_expr(args$expression_path)
  expr_log <- log2(expr + 1)
  score <- publication0_safe_gsva(expr_log, signature_list, method = "gsva")
  score_dt <- as.data.table(t(score), keep.rownames = "ID")

  if (file.exists(args$metadata_path)) {
    metadata <- fread(args$metadata_path)
    score_dt <- merge(metadata, score_dt, by = "ID", all.y = TRUE)
  }

  fwrite(score_dt, file.path(args$output_dir, "tables", "publication0_ecotype_signature_scores.tsv"), sep = "\t")
}

summary_lines <- c(
  sprintf("DOI: %s", "10.1016/j.ccell.2026.01.015"),
  "Ported logic: Figure 4 EC signature derivation",
  sprintf("Signatures recovered: %d", length(signature_list)),
  sprintf("Marker cutoff p_val_adj < %s", args$p_adj_cutoff),
  sprintf("Output directory: %s", args$output_dir)
)
writeLines(summary_lines, file.path(args$output_dir, "tables", "publication0_ecotype_logic_summary.txt"))
