source("publication0_10.1016_j.ccell.2026.01.015/publication0_utils.R")

defaults <- list(
  expression_path = "C:/Users/difen/Rcode/ecotyper/data/tcga_stad_prepared/bulk_stad_data.txt",
  metadata_path = "C:/Users/difen/Rcode/ecotyper/MutationIntegration_STAD/tables/sample_ecotype_mutation_table.tsv",
  group_column = "Ecotype",
  output_dir = "C:/Users/difen/Rcode/ecotyper/publication0_10.1016_j.ccell.2026.01.015/output_mos_logic",
  min_logfc = "1",
  max_fdr = "0.05"
)
args <- publication0_parse_args(defaults)

publication0_dir_create(args$output_dir)
publication0_dir_create(file.path(args$output_dir, "tables"))

expr <- publication0_read_expr(args$expression_path)
metadata <- fread(args$metadata_path)

if (!(args$group_column %in% colnames(metadata))) {
  stop(sprintf("Grouping column not found: %s", args$group_column))
}

sample_ids <- intersect(colnames(expr), metadata$ID)
if (length(sample_ids) == 0) {
  stop("No overlapping sample IDs between expression matrix and metadata.")
}

expr <- expr[, sample_ids, drop = FALSE]
metadata <- metadata[match(sample_ids, metadata$ID)]
groups <- metadata[[args$group_column]]

min_logfc <- as.numeric(args$min_logfc)
max_fdr <- as.numeric(args$max_fdr)

de_list <- publication0_one_vs_rest(expr = expr, groups = groups, min_logfc = min_logfc, max_fdr = max_fdr)
top_labels <- publication0_build_top_expression_labels(expr, groups)

signature_list <- lapply(names(de_list), function(group_name) {
  de_genes <- de_list[[group_name]]$gene
  top_genes <- top_labels$gene[top_labels$top_group == group_name]
  intersect(de_genes, top_genes)
})
names(signature_list) <- names(de_list)
signature_list <- signature_list[lengths(signature_list) > 0]

signature_table <- rbindlist(lapply(names(signature_list), function(group_name) {
  data.table(
    Signature = group_name,
    Gene = signature_list[[group_name]]
  )
}))

de_table <- rbindlist(de_list, use.names = TRUE, fill = TRUE)

fwrite(de_table, file.path(args$output_dir, "tables", "publication0_mos_one_vs_rest_de.tsv"), sep = "\t")
fwrite(signature_table, file.path(args$output_dir, "tables", "publication0_mos_signature_table.tsv"), sep = "\t")
saveRDS(signature_list, file.path(args$output_dir, "tables", "publication0_mos_signatures.rds"))

if (length(signature_list) > 0) {
  expr_log <- log2(expr + 1)
  score <- publication0_safe_gsva(expr_log, signature_list, method = "gsva")
  score_dt <- as.data.table(t(score), keep.rownames = "ID")
  score_dt <- merge(metadata, score_dt, by = "ID", all.y = TRUE)
  fwrite(score_dt, file.path(args$output_dir, "tables", "publication0_mos_signature_scores.tsv"), sep = "\t")
}

summary_lines <- c(
  sprintf("DOI: %s", "10.1016/j.ccell.2026.01.015"),
  "Ported logic: Figure 8 MOS signature derivation",
  sprintf("Grouping column used: %s", args$group_column),
  sprintf("Signatures recovered: %d", length(signature_list)),
  sprintf("DE filter: logFC > %s and FDR < %s", args$min_logfc, args$max_fdr),
  sprintf("Output directory: %s", args$output_dir)
)
writeLines(summary_lines, file.path(args$output_dir, "tables", "publication0_mos_logic_summary.txt"))
