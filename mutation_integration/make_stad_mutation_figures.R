source("mutation_integration/mutation_utils.R")
suppressPackageStartupMessages({
  library(ggplot2)
})

defaults <- list(
  association_path = "C:/Users/difen/Rcode/ecotyper/MutationIntegration_STAD/tables/ecotype_mutation_associations.tsv",
  sample_table_path = "C:/Users/difen/Rcode/ecotyper/MutationIntegration_STAD/tables/sample_ecotype_mutation_table.tsv",
  output_dir = "C:/Users/difen/Rcode/ecotyper/MutationIntegration_STAD"
)
args <- parse_simple_args(defaults)

dir_create_safe(file.path(args$output_dir, "plots"))

assoc <- fread(args$association_path)
sample_table <- fread(args$sample_table_path)

heatmap_df <- assoc[FeatureType == "binary"][order(FDR, -abs(EffectSize))]
heatmap_df <- heatmap_df[, .SD[1:min(.N, 8)], by = Ecotype]
heatmap_df[, Score := -log10(pmax(FDR, 1e-300)) * sign(log2(pmax(EffectSize, 1e-6)))]

plt_heatmap <- ggplot(heatmap_df, aes(x = Feature, y = Ecotype, fill = Score)) +
  geom_tile(color = "white", linewidth = 0.4) +
  scale_fill_gradient2(low = "#2166ac", mid = "#f7f7f7", high = "#b2182b", midpoint = 0) +
  labs(
    title = "STAD mutation biomarkers across EcoTyper ecotypes",
    x = "Mutation feature",
    y = "Ecotype",
    fill = "Signed\nsignal"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold")
  )

ggsave(file.path(args$output_dir, "plots", "ecotype_mutation_heatmap.png"), plt_heatmap, width = 10.5, height = 4.6, dpi = 300)

burden_plot <- ggplot(sample_table, aes(x = Ecotype, y = nonsilent_mutation_count, fill = Ecotype)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.85) +
  geom_jitter(width = 0.18, alpha = 0.35, size = 1) +
  scale_y_continuous(trans = "log1p") +
  guides(fill = "none") +
  labs(
    title = "Nonsilent mutation burden by EcoTyper ecotype",
    x = "Ecotype",
    y = "Nonsilent mutation count (log1p scale)"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(args$output_dir, "plots", "ecotype_mutation_burden_boxplot.png"), burden_plot, width = 8.5, height = 4.5, dpi = 300)
