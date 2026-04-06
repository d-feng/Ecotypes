source("mutation_integration/mutation_utils.R")
suppressPackageStartupMessages({
  library(cluster)
})

defaults <- list(
  sample_table_path = "C:/Users/difen/Rcode/ecotyper/MutationIntegration_STAD/tables/sample_ecotype_mutation_table.tsv",
  output_dir = "C:/Users/difen/Rcode/ecotyper/MutationIntegration_STAD",
  min_samples = "30",
  max_k = "4"
)
args <- parse_simple_args(defaults)

dir_create_safe(args$output_dir)
dir_create_safe(file.path(args$output_dir, "tables", "within_ecotype"))

dat <- fread(args$sample_table_path)
dat <- dat[!is.na(Ecotype)]

feature_cols <- names(dat)[startsWith(names(dat), "gene_") | startsWith(names(dat), "pathway_") |
  startsWith(names(dat), "flag_") | names(dat) %in% c("log_nonsilent_mutation_count", "log_mutated_gene_count")]

min_samples <- as.integer(args$min_samples)
max_k <- as.integer(args$max_k)
cluster_summaries <- list()

for (ecotype in sort(unique(dat$Ecotype))) {
  sub <- dat[Ecotype == ecotype]
  if (nrow(sub) < min_samples) {
    next
  }

  feature_dt <- sub[, ..feature_cols]
  variable_cols <- feature_cols[vapply(feature_dt, function(x) length(unique(x)) > 1, logical(1))]
  if (length(variable_cols) < 2) {
    next
  }

  x <- as.data.frame(feature_dt[, ..variable_cols])
  rownames(x) <- sub$ID
  binary_cols <- which(vapply(x, function(col) all(col %in% c(0, 1)), logical(1)))
  daisy_type <- if (length(binary_cols) > 0) list(asymm = binary_cols) else list()
  dist_obj <- daisy(x, metric = "gower", type = daisy_type)
  hc <- hclust(as.dist(dist_obj), method = "average")

  k_grid <- 2:min(max_k, nrow(sub) - 1)
  if (length(k_grid) == 0) {
    next
  }

  sil_scores <- vapply(k_grid, function(k) {
    clusters <- cutree(hc, k = k)
    mean(silhouette(clusters, dist_obj)[, "sil_width"])
  }, numeric(1))

  best_k <- k_grid[[which.max(sil_scores)]]
  best_clusters <- cutree(hc, k = best_k)

  cluster_table <- data.table(
    ID = names(best_clusters),
    Ecotype = ecotype,
    MutationSubcluster = paste0(ecotype, "_M", best_clusters)
  )

  fwrite(
    cluster_table,
    file.path(args$output_dir, "tables", "within_ecotype", sprintf("%s_mutation_subclusters.tsv", ecotype)),
    sep = "\t"
  )

  cluster_summaries[[ecotype]] <- data.table(
    Ecotype = ecotype,
    Samples = nrow(sub),
    VariableFeatures = length(variable_cols),
    BestK = best_k,
    MeanSilhouette = max(sil_scores)
  )
}

if (length(cluster_summaries) > 0) {
  fwrite(rbindlist(cluster_summaries), file.path(args$output_dir, "tables", "within_ecotype_clustering_summary.tsv"), sep = "\t")
}
