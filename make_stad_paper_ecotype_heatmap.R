suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

base_dir <- "C:/Users/difen/Rcode/ecotyper/RecoveryOutput_STAD/bulk_stad_data"
out_dir <- "C:/Users/difen/Rcode/ecotyper/comparison_outputs"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

read_matrix_with_rownames <- function(path) {
  lines <- readLines(path, warn = FALSE)
  if (length(lines) < 2) {
    stop(sprintf("File has no matrix data: %s", path))
  }

  header <- strsplit(lines[[1]], "\t", fixed = TRUE)[[1]]
  body <- strsplit(lines[-1], "\t", fixed = TRUE)

  row_names <- vapply(body, function(x) gsub('"', "", x[[1]], fixed = TRUE), character(1))
  value_matrix <- do.call(rbind, lapply(body, function(x) as.numeric(x[-1])))

  colnames(value_matrix) <- gsub('"', "", header, fixed = TRUE)
  rownames(value_matrix) <- row_names

  as.data.frame(value_matrix, check.names = FALSE)
}

ecotype_abundance <- read_matrix_with_rownames(file.path(base_dir, "Ecotypes", "ecotype_abundance.txt"))

sample_to_ce <- vapply(colnames(ecotype_abundance), function(sample_id) {
  vals <- ecotype_abundance[, sample_id]
  rownames(ecotype_abundance)[which.max(vals)]
}, character(1))

celltype_dirs <- setdiff(list.dirs(base_dir, recursive = FALSE, full.names = FALSE), "Ecotypes")
celltype_means <- list()

for (celltype in celltype_dirs) {
  state_abundance <- read_matrix_with_rownames(file.path(base_dir, celltype, "state_abundances.txt"))
  sample_totals <- colSums(state_abundance, na.rm = TRUE)
  for (sample_id in names(sample_totals)) {
    ce <- sample_to_ce[[sample_id]]
    if (is.null(celltype_means[[ce]])) {
      celltype_means[[ce]] <- list()
    }
    if (is.null(celltype_means[[ce]][[celltype]])) {
      celltype_means[[ce]][[celltype]] <- c()
    }
    celltype_means[[ce]][[celltype]] <- c(celltype_means[[ce]][[celltype]], sample_totals[[sample_id]])
  }
}

ce_profiles <- lapply(celltype_means, function(ct_list) {
  prof <- vapply(ct_list, function(x) mean(x, na.rm = TRUE), numeric(1))
  prof / sum(prof)
})

paper_scores <- list(
  "EC1\nT cell activation" = c("CD8.T.cells" = 1, "CD4.T.cells" = 1, "NK.cells" = 0.5),
  "EC2\nTLS / B cell" = c("B.cells" = 1, "PCs" = 1, "Dendritic.cells" = 0.5),
  "EC3\nVascular" = c("Endothelial.cells" = 1),
  "EC4\nECM / fibroblast" = c("Fibroblasts" = 1),
  "EC5\nMacrophage suppressive" = c("Monocytes.and.Macrophages" = 1, "PMNs" = 0.5, "Mast.cells" = 0.25)
)

ce_order <- rownames(ecotype_abundance)[order(as.integer(sub("CE", "", rownames(ecotype_abundance))))]
paper_order <- names(paper_scores)

score_df <- rbindlist(lapply(paper_order, function(paper_name) {
  weights <- paper_scores[[paper_name]]
  max_possible <- sum(weights)
  rbindlist(lapply(ce_order, function(ce) {
    prof <- ce_profiles[[ce]]
    score <- 0
    for (ct in names(weights)) {
      score <- score + weights[[ct]] * ifelse(ct %in% names(prof), prof[[ct]], 0)
    }
    data.table(PaperEcotype = paper_name, EcoTyperEcotype = ce, Similarity = score / max_possible)
  }))
}))

fwrite(score_df, file.path(out_dir, "STAD_paper_vs_EcoTyper_similarity.tsv"), sep = "\t")

score_wide <- dcast(score_df, PaperEcotype ~ EcoTyperEcotype, value.var = "Similarity")
score_mat <- as.matrix(score_wide[, -1])
rownames(score_mat) <- score_wide[[1]]

row_order <- rownames(score_mat)[hclust(dist(score_mat))$order]
col_order <- colnames(score_mat)[hclust(dist(t(score_mat)))$order]

paper_labels <- c(
  "EC1\nT cell activation" = "EC1  T cell activation",
  "EC2\nTLS / B cell" = "EC2  TLS / B cell",
  "EC3\nVascular" = "EC3  Vascular",
  "EC4\nECM / fibroblast" = "EC4  ECM / fibroblast",
  "EC5\nMacrophage suppressive" = "EC5  Macrophage suppressive"
)

base_palette <- scale_fill_gradient(
  low = "#f7fbff",
  high = "#084081"
)

pub_palette <- scale_fill_gradientn(
  colours = c("#fffaf2", "#fdd49e", "#fc8d59", "#b30000")
)

score_df$PaperEcotype <- factor(score_df$PaperEcotype, levels = paper_order)
score_df$EcoTyperEcotype <- factor(score_df$EcoTyperEcotype, levels = ce_order)

plt <- ggplot(score_df, aes(x = EcoTyperEcotype, y = PaperEcotype, fill = Similarity)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.3f", Similarity)), size = 3.2) +
  base_palette +
  labs(
    title = "TCGA-STAD: Gastric Paper Ecotypes vs EcoTyper Ecotypes",
    x = "EcoTyper ecotypes",
    y = "Paper ecotypes"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold")
  )

ggsave(file.path(out_dir, "STAD_paper_vs_EcoTyper_similarity_heatmap.png"), plt, width = 11, height = 4.8, dpi = 300)

clustered_df <- copy(score_df)
clustered_df$PaperEcotype <- factor(clustered_df$PaperEcotype, levels = row_order)
clustered_df$EcoTyperEcotype <- factor(clustered_df$EcoTyperEcotype, levels = col_order)

plt_clustered <- ggplot(clustered_df, aes(x = EcoTyperEcotype, y = PaperEcotype, fill = Similarity)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.3f", Similarity)), size = 3.4, color = "black") +
  pub_palette +
  scale_y_discrete(labels = paper_labels) +
  labs(
    title = "TCGA-STAD Ecotype Similarity Crosswalk",
    subtitle = "Rows and columns reordered by hierarchical clustering",
    x = "EcoTyper ecotypes",
    y = "Gastric paper ecotypes",
    fill = "Similarity"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", color = "#222222"),
    axis.text.y = element_text(face = "bold", color = "#222222"),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(color = "#444444"),
    legend.title = element_text(face = "bold")
  )

ggsave(
  file.path(out_dir, "STAD_paper_vs_EcoTyper_similarity_heatmap_clustered.png"),
  plt_clustered,
  width = 11,
  height = 5.2,
  dpi = 300
)

plt_journal <- ggplot(clustered_df, aes(x = EcoTyperEcotype, y = PaperEcotype, fill = Similarity)) +
  geom_tile(color = "white", linewidth = 0.7) +
  pub_palette +
  scale_y_discrete(labels = paper_labels) +
  labs(
    title = "Ecotype Similarity Crosswalk",
    subtitle = "TCGA-STAD projected onto EcoTyper carcinoma ecotypes",
    x = "EcoTyper ecotypes",
    y = "Gastric paper ecotypes",
    fill = "Similarity"
  ) +
  theme_minimal(base_size = 15) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 12, color = "#1f1f1f"),
    axis.text.y = element_text(face = "bold", size = 12, color = "#1f1f1f"),
    axis.title = element_text(face = "bold", size = 13),
    plot.title = element_text(face = "bold", size = 18, color = "#111111"),
    plot.subtitle = element_text(size = 12, color = "#4d4d4d"),
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 11),
    plot.margin = margin(10, 14, 10, 10)
  )

ggsave(
  file.path(out_dir, "STAD_paper_vs_EcoTyper_similarity_heatmap_journal.png"),
  plt_journal,
  width = 10.5,
  height = 5.4,
  dpi = 400
)

ggsave(
  file.path(out_dir, "STAD_paper_vs_EcoTyper_similarity_heatmap_journal.pdf"),
  plt_journal,
  width = 10.5,
  height = 5.4,
  device = cairo_pdf
)
