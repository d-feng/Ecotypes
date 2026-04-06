suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(grid)
})

source("mutation_integration/mutation_utils.R")

defaults <- list(
  sample_table_path = "C:/Users/difen/Rcode/ecotyper/MutationIntegration_STAD/tables/sample_ecotype_mutation_table.tsv",
  output_dir = "C:/Users/difen/Rcode/ecotyper/MutationIntegration_STAD/plots"
)
args <- parse_simple_args(defaults)

dir_create_safe(args$output_dir)

sample_table <- fread(args$sample_table_path)
sample_table <- sample_table[!is.na(Ecotype)]

ecotype_levels <- sort(unique(sample_table$Ecotype))
ecotype_levels <- ecotype_levels[order(as.integer(sub("CE", "", ecotype_levels)))]
sample_table[, Ecotype := factor(Ecotype, levels = ecotype_levels)]
setorder(sample_table, Ecotype, -nonsilent_mutation_count, -mutated_gene_count, ID)
sample_table[, sample_pos := seq_len(.N)]

ecotype_colors <- c(
  CE1 = "#5B8FF9",
  CE2 = "#61DDAA",
  CE3 = "#65789B",
  CE4 = "#F6BD16",
  CE5 = "#7262FD",
  CE6 = "#78D3F8",
  CE7 = "#9661BC",
  CE8 = "#F6903D",
  CE9 = "#008685",
  CE10 = "#F08BB4"
)
ecotype_colors <- ecotype_colors[names(ecotype_colors) %in% ecotype_levels]

gene_map <- c(
  gene_TP53 = "TP53",
  gene_CDH1 = "CDH1",
  gene_ARID1A = "ARID1A",
  gene_PIK3CA = "PIK3CA",
  gene_ERBB2 = "ERBB2",
  gene_KRAS = "KRAS",
  gene_RHOA = "RHOA",
  gene_FGFR2 = "FGFR2"
)

oncoprint_long <- melt(
  sample_table[, c("ID", "Ecotype", "sample_pos", names(gene_map)), with = FALSE],
  id.vars = c("ID", "Ecotype", "sample_pos"),
  variable.name = "Feature",
  value.name = "Mutated"
)
oncoprint_long[, Gene := factor(gene_map[Feature], levels = rev(unname(gene_map)))]
oncoprint_long[, Mutated := factor(ifelse(Mutated == 1, "Altered", "Wild-type"), levels = c("Wild-type", "Altered"))]

burden_plot <- ggplot(sample_table, aes(x = sample_pos, y = log1p(nonsilent_mutation_count), fill = Ecotype)) +
  geom_col(width = 0.95, color = NA) +
  scale_fill_manual(values = ecotype_colors) +
  labs(y = "log1p\nburden", x = NULL) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    plot.margin = margin(4, 8, 0, 8)
  )

ecotype_strip <- ggplot(sample_table, aes(x = sample_pos, y = 1, fill = Ecotype)) +
  geom_tile(height = 1, width = 0.98) +
  scale_fill_manual(values = ecotype_colors) +
  labs(x = NULL, y = NULL, fill = "Ecotype") +
  theme_void(base_size = 11) +
  theme(
    legend.position = "none",
    plot.margin = margin(0, 8, 0, 8)
  )

oncoprint_plot <- ggplot(oncoprint_long, aes(x = sample_pos, y = Gene, fill = Mutated)) +
  geom_tile(color = "white", linewidth = 0.15, height = 0.92, width = 0.98) +
  scale_fill_manual(values = c("Wild-type" = "#f1efe8", "Altered" = "#b2182b")) +
  labs(
    title = "TCGA-STAD Gastric Oncoprint Grouped by EcoTyper Ecotype",
    subtitle = "Samples ordered by ecotype and nonsilent mutation burden",
    x = "Samples",
    y = NULL,
    fill = "Status"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(face = "bold", color = "#1f1f1f"),
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 11, color = "#4d4d4d"),
    legend.position = "top",
    plot.margin = margin(0, 8, 8, 8)
  )

render_composite_plot <- function(file_path, width, height) {
  if (grepl("\\.png$", file_path, ignore.case = TRUE)) {
    png(filename = file_path, width = width, height = height, res = 300, units = "in")
  } else if (grepl("\\.pdf$", file_path, ignore.case = TRUE)) {
    pdf(file = file_path, width = width, height = height, useDingbats = FALSE)
  } else {
    stop("Unsupported output format")
  }

  grid.newpage()
  pushViewport(viewport(layout = grid.layout(nrow = 3, ncol = 1, heights = unit(c(1.1, 0.28, 3.7), "null"))))
  print(burden_plot, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(ecotype_strip, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
  print(oncoprint_plot, vp = viewport(layout.pos.row = 3, layout.pos.col = 1))
  dev.off()
}

render_composite_plot(file.path(args$output_dir, "stad_ecotype_oncoprint_gastric.png"), width = 12, height = 5.8)
render_composite_plot(file.path(args$output_dir, "stad_ecotype_oncoprint_gastric.pdf"), width = 12, height = 5.8)

paper_features <- c(
  gene_TP53 = "TP53",
  gene_CDH1 = "CDH1",
  gene_ARID1A = "ARID1A",
  gene_PIK3CA = "PIK3CA",
  pathway_DNA_REPAIR = "DNA repair",
  pathway_FGFR = "FGFR pathway",
  pathway_TGF_BETA = "TGF-beta pathway",
  pathway_PI3K = "PI3K pathway"
)

panel_dt <- melt(
  sample_table[, c("Ecotype", names(paper_features)), with = FALSE],
  id.vars = "Ecotype",
  variable.name = "Feature",
  value.name = "Value"
)

panel_summary <- panel_dt[, .(
  Prevalence = mean(Value, na.rm = TRUE),
  Samples = .N
), by = .(Ecotype, Feature)]

cohort_mean <- panel_summary[, .(CohortMean = weighted.mean(Prevalence, Samples)), by = Feature]
panel_summary <- merge(panel_summary, cohort_mean, by = "Feature", all.x = TRUE)
panel_summary[, DeltaFromCohort := Prevalence - CohortMean]
panel_summary[, FeatureLabel := factor(paper_features[Feature], levels = rev(unname(paper_features)))]
panel_summary[, Ecotype := factor(Ecotype, levels = ecotype_levels)]

paper_panel <- ggplot(panel_summary, aes(x = Ecotype, y = FeatureLabel, fill = DeltaFromCohort)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%1.0f%%", 100 * Prevalence)), size = 3.6) +
  scale_fill_gradient2(
    low = "#2166ac",
    mid = "#f7f7f7",
    high = "#b2182b",
    midpoint = 0
  ) +
  labs(
    title = "Paper-Focused Gastric Mutation Biology Across EcoTyper Ecotypes",
    subtitle = "Cell labels show prevalence within each CE group; color shows deviation from the cohort-wide average",
    x = "EcoTyper ecotype",
    y = NULL,
    fill = "Delta vs\ncohort"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", color = "#1f1f1f"),
    axis.text.y = element_text(face = "bold", color = "#1f1f1f"),
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 11, color = "#4d4d4d"),
    legend.position = "right",
    plot.margin = margin(8, 8, 8, 8)
  )

ggsave(
  file.path(args$output_dir, "stad_paper_focused_comparison_panel.png"),
  paper_panel,
  width = 9.5,
  height = 5.6,
  dpi = 350
)

ggsave(
  file.path(args$output_dir, "stad_paper_focused_comparison_panel.pdf"),
  paper_panel,
  width = 9.5,
  height = 5.6,
  device = cairo_pdf
)
