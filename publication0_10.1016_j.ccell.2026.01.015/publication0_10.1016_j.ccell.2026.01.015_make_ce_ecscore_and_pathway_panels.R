suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(grid)
})

source("publication0_10.1016_j.ccell.2026.01.015/publication0_utils.R")

defaults <- list(
  expression_path = "C:/Users/difen/Rcode/ecotyper/data/tcga_stad_prepared/bulk_stad_data.txt",
  metadata_path = "C:/Users/difen/Rcode/ecotyper/MutationIntegration_STAD/tables/sample_ecotype_mutation_table.tsv",
  output_dir = "C:/Users/difen/Rcode/ecotyper/publication0_10.1016_j.ccell.2026.01.015/output_stad_ce_panels"
)
args <- publication0_parse_args(defaults)

publication0_dir_create(args$output_dir)
publication0_dir_create(file.path(args$output_dir, "tables"))
publication0_dir_create(file.path(args$output_dir, "plots"))

publication0_ec_signatures <- list(
  EC1_ImmuneActivation = c("CD3D", "CD3E", "CD8A", "CD8B", "IFNG", "GZMB", "PRF1", "NKG7", "CXCL9", "CXCL10", "TBX21", "STAT1", "ITGAE", "PDCD1"),
  EC2_TLSImmunity = c("CD19", "MS4A1", "CD79A", "CD79B", "CXCR5", "CXCL13", "TNFRSF13B", "MZB1", "JCHAIN", "BANK1", "CD27", "CD74", "HLA.DRA", "HLA.DPA1"),
  EC3_AngiogenesisTME = c("KDR", "FLT1", "FLT4", "PECAM1", "VWF", "EMCN", "ESAM", "ANGPT2", "VEGFA", "ACKR1", "RGS5", "CDH5", "SEMA3G"),
  EC4_ECMOrganization = c("COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL5A2", "COL6A1", "COL6A2", "COL6A3", "DCN", "LUM", "FAP", "POSTN", "MMP2", "TGFB1", "THY1"),
  EC5_MetabolismTME = c("APOA1", "APOE", "FABP1", "FABP4", "CPT1A", "ACADM", "ATP5F1A", "NDUFS1", "COX5A", "PPARGC1A", "TREM2", "CD68", "CD163", "LAMP3", "CCL2")
)

publication0_pathway_panel <- list(
  T_Cell_Activation = c("CD3D", "CD3E", "CD8A", "CD8B", "PRF1", "GZMB", "NKG7", "TBX21", "ITGAE"),
  IFNG_Response = c("IFNG", "STAT1", "IRF1", "CXCL9", "CXCL10", "IDO1", "GBP1", "HLA.DRA"),
  TLS_BCell_Immunity = c("CD19", "MS4A1", "CD79A", "CD79B", "CXCR5", "CXCL13", "MZB1", "JCHAIN"),
  Antigen_Presentation = c("HLA.DRA", "HLA.DPA1", "HLA.DPB1", "HLA.B", "B2M", "TAP1", "TAP2", "CIITA", "CD74"),
  Angiogenesis = c("VEGFA", "KDR", "FLT1", "ANGPT2", "PECAM1", "VWF", "ESAM", "EMCN"),
  ECM_Fibroblast = c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "FAP", "POSTN", "TGFB1", "THY1"),
  OxPhos_Metabolism = c("NDUFA1", "NDUFS1", "UQCRC1", "COX5A", "COX6C", "ATP5F1A", "ATP5F1B", "ATP5MC1", "SDHA", "SDHB"),
  Myeloid_Chemotaxis = c("CCL2", "CCL7", "CCL8", "CXCL8", "S100A8", "S100A9", "CSF1", "IL1B", "TREM2", "APOE"),
  FGFR_Signaling = c("FGFR1", "FGFR2", "FGFR3", "FGF3", "FGF4", "FGF7", "FGF9", "FGF19", "FRS2", "KLB"),
  APOA1_TREM2_Axis = c("APOA1", "APOE", "TREM2", "CD68", "CD163", "LAMP3", "MARCO", "C1QA", "C1QB")
)

expr <- publication0_read_expr(args$expression_path)
metadata <- fread(args$metadata_path)
metadata <- metadata[!is.na(Ecotype)]

sample_ids <- intersect(colnames(expr), metadata$ID)
expr <- expr[, sample_ids, drop = FALSE]
metadata <- metadata[match(sample_ids, metadata$ID)]
metadata$Ecotype <- factor(metadata$Ecotype, levels = sort(unique(metadata$Ecotype), decreasing = FALSE))

filter_signatures <- function(signature_list, expr_rownames) {
  filtered <- lapply(signature_list, function(g) {
    unique(g[g %in% expr_rownames])
  })
  filtered[lengths(filtered) > 1]
}

ec_sigs_filtered <- filter_signatures(publication0_ec_signatures, rownames(expr))
pathway_sigs_filtered <- filter_signatures(publication0_pathway_panel, rownames(expr))

expr_log <- log2(expr + 1)
ec_scores <- publication0_safe_gsva(expr_log, ec_sigs_filtered, method = "gsva")
pathway_scores <- publication0_safe_gsva(expr_log, pathway_sigs_filtered, method = "gsva")

score_to_ce_summary <- function(score_matrix, metadata, ce_levels) {
  dt <- as.data.table(t(score_matrix), keep.rownames = "ID")
  dt <- merge(metadata[, .(ID, Ecotype)], dt, by = "ID", all.x = TRUE)
  value_cols <- setdiff(colnames(dt), c("ID", "Ecotype"))
  summary_dt <- dt[, lapply(.SD, mean, na.rm = TRUE), by = Ecotype, .SDcols = value_cols]
  summary_dt$Ecotype <- factor(summary_dt$Ecotype, levels = ce_levels)
  summary_dt
}

ce_levels <- levels(metadata$Ecotype)
ec_summary <- score_to_ce_summary(ec_scores, metadata, ce_levels)
pathway_summary <- score_to_ce_summary(pathway_scores, metadata, ce_levels)

zscore_by_signature <- function(summary_dt) {
  mat <- as.matrix(summary_dt[, -1])
  rownames(mat) <- as.character(summary_dt$Ecotype)
  scaled <- scale(mat, center = TRUE, scale = TRUE)
  scaled[is.na(scaled)] <- 0
  as.data.table(scaled, keep.rownames = "Ecotype")
}

ec_summary_z <- zscore_by_signature(ec_summary)
pathway_summary_z <- zscore_by_signature(pathway_summary)

ec_long <- melt(ec_summary_z, id.vars = "Ecotype", variable.name = "Signature", value.name = "Score")
pathway_long <- melt(pathway_summary_z, id.vars = "Ecotype", variable.name = "Pathway", value.name = "Score")

ec_label_map <- c(
  EC1_ImmuneActivation = "EC1\nImmune activation",
  EC2_TLSImmunity = "EC2\nTLS immunity",
  EC3_AngiogenesisTME = "EC3\nAngiogenesis TME",
  EC4_ECMOrganization = "EC4\nECM organization",
  EC5_MetabolismTME = "EC5\nMetabolism TME"
)

pathway_label_map <- c(
  T_Cell_Activation = "T-cell\nactivation",
  IFNG_Response = "IFNG\nresponse",
  TLS_BCell_Immunity = "TLS / B-cell\nimmunity",
  Antigen_Presentation = "Antigen\npresentation",
  Angiogenesis = "Angiogenesis",
  ECM_Fibroblast = "ECM /\nFibroblast",
  OxPhos_Metabolism = "OxPhos /\nMetabolism",
  Myeloid_Chemotaxis = "Myeloid\nchemotaxis",
  FGFR_Signaling = "FGFR\nsignaling",
  APOA1_TREM2_Axis = "APOA1-\nTREM2 axis"
)

ec_long[, Signature := factor(Signature, levels = names(ec_sigs_filtered))]
pathway_long[, Pathway := factor(Pathway, levels = names(pathway_sigs_filtered))]
ec_long[, Ecotype := factor(Ecotype, levels = ce_levels)]
pathway_long[, Ecotype := factor(Ecotype, levels = ce_levels)]

fwrite(ec_summary, file.path(args$output_dir, "tables", "publication0_ce_ec_score_means.tsv"), sep = "\t")
fwrite(pathway_summary, file.path(args$output_dir, "tables", "publication0_ce_pathway_score_means.tsv"), sep = "\t")
fwrite(rbindlist(lapply(names(ec_sigs_filtered), function(nm) data.table(Set = nm, Gene = ec_sigs_filtered[[nm]]))), file.path(args$output_dir, "tables", "publication0_ce_ec_signature_genes_used.tsv"), sep = "\t")
fwrite(rbindlist(lapply(names(pathway_sigs_filtered), function(nm) data.table(Set = nm, Gene = pathway_sigs_filtered[[nm]]))), file.path(args$output_dir, "tables", "publication0_ce_pathway_genes_used.tsv"), sep = "\t")

tile_palette <- scale_fill_gradient2(low = "#2166ac", mid = "#f7f7f7", high = "#b2182b", midpoint = 0)

ec_heatmap <- ggplot(ec_long, aes(x = Signature, y = Ecotype, fill = Score)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.2f", Score)), size = 3.2) +
  tile_palette +
  scale_x_discrete(labels = ec_label_map) +
  labs(
    title = "CE x publication0 EC1-EC5 score heatmap",
    subtitle = "TCGA-STAD bulk RNA scored with publication-inspired EC signatures; values are z-scored CE means",
    x = NULL,
    y = "EcoTyper CE",
    fill = "Z-score"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(face = "bold", color = "#1f1f1f"),
    axis.text.y = element_text(face = "bold", color = "#1f1f1f"),
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 11, color = "#4d4d4d")
  )

pathway_panel <- ggplot(pathway_long, aes(x = Pathway, y = Ecotype, fill = Score)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.2f", Score)), size = 2.8) +
  tile_palette +
  scale_x_discrete(labels = pathway_label_map) +
  labs(
    title = "CE x pathway annotation panel",
    subtitle = "Publication-inspired pathway programs scored in TCGA-STAD and summarized by EcoTyper ecotype",
    x = NULL,
    y = "EcoTyper CE",
    fill = "Z-score"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(face = "bold", color = "#1f1f1f"),
    axis.text.y = element_text(face = "bold", color = "#1f1f1f"),
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 11, color = "#4d4d4d")
  )

ggsave(file.path(args$output_dir, "plots", "publication0_ce_ecscore_heatmap.png"), ec_heatmap, width = 8.8, height = 5.4, dpi = 350)
ggsave(file.path(args$output_dir, "plots", "publication0_ce_ecscore_heatmap.pdf"), ec_heatmap, width = 8.8, height = 5.4, device = cairo_pdf)

ggsave(file.path(args$output_dir, "plots", "publication0_ce_pathway_annotation_panel.png"), pathway_panel, width = 12, height = 5.6, dpi = 350)
ggsave(file.path(args$output_dir, "plots", "publication0_ce_pathway_annotation_panel.pdf"), pathway_panel, width = 12, height = 5.6, device = cairo_pdf)

png(file.path(args$output_dir, "plots", "publication0_ce_combined_panels.png"), width = 12, height = 10, units = "in", res = 320)
grid.newpage()
pushViewport(viewport(layout = grid.layout(nrow = 2, ncol = 1, heights = unit(c(1, 1.2), "null"))))
print(ec_heatmap, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(pathway_panel, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
dev.off()

pdf(file.path(args$output_dir, "plots", "publication0_ce_combined_panels.pdf"), width = 12, height = 10, useDingbats = FALSE)
grid.newpage()
pushViewport(viewport(layout = grid.layout(nrow = 2, ncol = 1, heights = unit(c(1, 1.2), "null"))))
print(ec_heatmap, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(pathway_panel, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
dev.off()

summary_lines <- c(
  "publication0 STAD CE panel generation complete",
  "Note: EC1-EC5 signatures are publication-inspired proxy signatures because the authors' exact final EC marker table was not publicly available.",
  sprintf("Samples scored: %d", length(sample_ids)),
  sprintf("CE groups summarized: %d", length(ce_levels)),
  sprintf("EC signatures used: %d", length(ec_sigs_filtered)),
  sprintf("Pathway programs used: %d", length(pathway_sigs_filtered))
)
writeLines(summary_lines, file.path(args$output_dir, "tables", "publication0_ce_panels_summary.txt"))
