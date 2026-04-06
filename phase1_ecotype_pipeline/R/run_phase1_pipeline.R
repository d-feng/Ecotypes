suppressPackageStartupMessages({
  library(data.table)
  library(yaml)
})

source(file.path("R", "pipeline_utils.R"))

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  config_path <- NULL
  if (length(args) >= 2) {
    for (i in seq(1, length(args) - 1, by = 2)) {
      if (args[i] == "--config") {
        config_path <- args[i + 1]
      }
    }
  }
  if (is.null(config_path)) {
    stop("Usage: Rscript R/run_phase1_pipeline.R --config config/example_config.yml")
  }
  config_path
}

main <- function() {
  config_path <- parse_args()
  config <- yaml.load_file(config_path)

  out_dir <- config$input$output_dir
  tables_dir <- ensure_dir(file.path(out_dir, "tables"))
  plots_dir <- ensure_dir(file.path(out_dir, "plots"))

  expression <- read_matrix_tsv(config$input$expression_path)
  metadata <- read_metadata_tsv(config$input$metadata_path)
  state_abundance <- read_matrix_tsv(config$input$cell_state_abundances_path)

  aligned <- align_inputs(expression, metadata, state_abundance)
  expression <- aligned$expression
  metadata <- aligned$metadata
  state_abundance <- normalize_state_abundance(aligned$state_abundance)

  expression_scaled <- preprocess_expression(expression, config$preprocessing)

  discovery <- discover_ecotypes_from_states(
    state_abundance = state_abundance,
    min_k = config$ecotypes$k_min,
    max_k = config$ecotypes$k_max,
    alpha = config$ecotypes$jaccard_significance_alpha,
    min_states_per_ecotype = config$ecotypes$min_states_per_ecotype
  )

  assignments <- assign_samples_to_ecotypes(state_abundance, discovery$state_assignments)
  sample_labels <- assignments$sample_assignments$Ecotype
  names(sample_labels) <- assignments$sample_assignments$ID

  marker_table <- derive_markers(
    expression_matrix = expression_scaled[, assignments$sample_assignments$ID, drop = FALSE],
    sample_labels = sample_labels[assignments$sample_assignments$ID],
    top_n = config$reporting$top_markers_per_ecotype
  )

  metadata_tests <- metadata_association_tests(
    metadata = metadata[match(assignments$sample_assignments$ID, metadata$SampleID), , drop = FALSE],
    sample_labels = sample_labels[assignments$sample_assignments$ID]
  )

  fwrite(as.data.table(discovery$selection$metrics), file.path(tables_dir, "ecotype_number_metrics.tsv"), sep = "\t")
  fwrite(as.data.table(discovery$state_assignments), file.path(tables_dir, "state_to_ecotype.tsv"), sep = "\t")
  fwrite(as.data.table(assignments$sample_assignments), file.path(tables_dir, "ecotype_assignments.tsv"), sep = "\t")
  fwrite(as.data.table(assignments$ecotype_abundance, keep.rownames = "SampleID"), file.path(tables_dir, "ecotype_abundance.tsv"), sep = "\t")
  fwrite(as.data.table(state_abundance, keep.rownames = "State"), file.path(tables_dir, "cell_state_abundance_normalized.tsv"), sep = "\t")
  fwrite(as.data.table(marker_table), file.path(tables_dir, "top_markers_by_ecotype.tsv"), sep = "\t")
  fwrite(as.data.table(metadata_tests), file.path(tables_dir, "metadata_associations.tsv"), sep = "\t")

  write_ecotyper_style_outputs(
    output_dir = out_dir,
    state_abundance = state_abundance,
    state_assignments = discovery$state_assignments,
    sample_assignments = assignments$sample_assignments,
    ecotype_abundance = assignments$ecotype_abundance
  )

  plot_state_pca(state_abundance, assignments$sample_assignments, file.path(plots_dir, "sample_pca_from_states.png"))
  save_similarity_heatmap(discovery$jaccard_sig, discovery$state_labels, file.path(plots_dir, "jaccard_state_heatmap.png"))
  plot_ecotype_number_curve(discovery$selection$metrics, file.path(plots_dir, "ecotype_number_selection.png"))
  save_state_ecotype_heatmap(state_abundance, discovery$state_assignments, file.path(plots_dir, "state_abundance_heatmap.png"))
  save_marker_heatmap(
    expression_matrix = expression_scaled[, assignments$sample_assignments$ID, drop = FALSE],
    marker_table = marker_table,
    sample_labels = sample_labels[assignments$sample_assignments$ID],
    out_path = file.path(plots_dir, "ecotype_marker_heatmap.png")
  )
  plot_ecotype_by_metadata(
    metadata = metadata[match(assignments$sample_assignments$ID, metadata$SampleID), , drop = FALSE],
    sample_labels = sample_labels[assignments$sample_assignments$ID],
    column_name = config$reporting$metadata_plot_column,
    out_path = file.path(plots_dir, "ecotype_by_metadata.png")
  )

  summary_lines <- c(
    sprintf("Best ecotype number: %s", discovery$selection$best_k),
    sprintf("Samples analyzed: %s", ncol(expression_scaled)),
    sprintf("States analyzed: %s", nrow(state_abundance)),
    sprintf("Genes retained after preprocessing: %s", nrow(expression_scaled)),
    sprintf("Jaccard significance alpha: %.3f", config$ecotypes$jaccard_significance_alpha)
  )
  writeLines(summary_lines, file.path(out_dir, "run_summary.txt"))
  cat(paste(summary_lines, collapse = "\n"), "\n")
}

main()
