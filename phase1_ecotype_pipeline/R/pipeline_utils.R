suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ComplexHeatmap)
  library(RColorBrewer)
  library(yaml)
  library(cluster)
})

read_matrix_tsv <- function(path, row_key = NULL) {
  dt <- fread(path)
  if (ncol(dt) < 2) {
    stop(sprintf("Matrix file '%s' must contain an identifier column and at least one sample column.", path))
  }
  id_col <- if (is.null(row_key)) names(dt)[1] else row_key
  if (!id_col %in% names(dt)) {
    stop(sprintf("Identifier column '%s' not found in '%s'.", id_col, path))
  }
  mat_df <- as.data.frame(dt)
  rownames(mat_df) <- mat_df[[id_col]]
  mat_df[[id_col]] <- NULL
  mat <- as.matrix(mat_df)
  storage.mode(mat) <- "numeric"
  mat
}

read_metadata_tsv <- function(path) {
  meta <- fread(path)
  if (!"SampleID" %in% names(meta)) {
    stop("Metadata file must contain a 'SampleID' column.")
  }
  as.data.frame(meta)
}

ensure_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

align_inputs <- function(expression, metadata, state_abundance) {
  sample_ids <- Reduce(intersect, list(colnames(expression), metadata$SampleID, colnames(state_abundance)))
  if (length(sample_ids) < 6) {
    stop("Fewer than 6 shared samples remain after aligning expression, metadata, and cell-state abundance matrix.")
  }
  list(
    expression = expression[, sample_ids, drop = FALSE],
    metadata = metadata[match(sample_ids, metadata$SampleID), , drop = FALSE],
    state_abundance = state_abundance[, sample_ids, drop = FALSE]
  )
}

preprocess_expression <- function(expression, settings) {
  gene_means <- rowMeans(expression, na.rm = TRUE)
  expressed_counts <- rowSums(expression > 0, na.rm = TRUE)
  keep <- gene_means >= settings$min_gene_mean & expressed_counts >= settings$min_samples_expressed
  filtered <- expression[keep, , drop = FALSE]
  if (nrow(filtered) < 50) {
    stop("Expression filtering retained fewer than 50 genes. Relax preprocessing thresholds.")
  }
  if (isTRUE(settings$log_transform)) {
    filtered <- log2(filtered + 1)
  }
  if (isTRUE(settings$zscore_genes)) {
    filtered <- t(scale(t(filtered)))
    filtered[is.na(filtered)] <- 0
  }
  filtered
}

normalize_state_abundance <- function(state_abundance) {
  state_abundance[state_abundance < 0] <- 0
  sample_sums <- colSums(state_abundance, na.rm = TRUE)
  sample_sums[sample_sums == 0] <- 1
  sweep(state_abundance, 2, sample_sums, "/")
}

parse_state_metadata <- function(state_names) {
  parts <- strsplit(state_names, "\\|")
  data.frame(
    State = state_names,
    CellType = vapply(parts, function(x) if (length(x) >= 1) x[1] else "Unknown", character(1)),
    StateLabel = vapply(parts, function(x) if (length(x) >= 2) x[2] else x[1], character(1)),
    stringsAsFactors = FALSE
  )
}

format_state_label <- function(x) {
  digits <- gsub("[^0-9]", "", x)
  ifelse(
    nzchar(digits),
    sprintf("S%02d", as.integer(digits)),
    x
  )
}

format_ecotype_label <- function(x, prefix = "CE") {
  digits <- gsub("[^0-9]", "", x)
  ifelse(
    nzchar(digits),
    sprintf("%s%d", prefix, as.integer(digits)),
    x
  )
}

split_state_abundance_by_cell_type <- function(state_abundance) {
  state_meta <- parse_state_metadata(rownames(state_abundance))
  splits <- split(state_meta$State, state_meta$CellType)
  lapply(splits, function(states) state_abundance[states, , drop = FALSE])
}

select_dominant_states <- function(state_abundance) {
  state_meta <- parse_state_metadata(rownames(state_abundance))
  binary <- matrix(
    0,
    nrow = nrow(state_abundance),
    ncol = ncol(state_abundance),
    dimnames = dimnames(state_abundance)
  )
  for (cell_type in unique(state_meta$CellType)) {
    idx <- which(state_meta$CellType == cell_type)
    block <- state_abundance[idx, , drop = FALSE]
    winners <- apply(block, 2, function(x) idx[which.max(x)])
    for (sample_idx in seq_len(ncol(state_abundance))) {
      binary[winners[sample_idx], sample_idx] <- 1
    }
  }
  list(binary = binary, state_meta = state_meta)
}

build_state_assignment_table <- function(state_abundance) {
  state_meta <- parse_state_metadata(rownames(state_abundance))
  rows <- list()
  idx <- 1L
  for (cell_type in unique(state_meta$CellType)) {
    state_ids <- state_meta$State[state_meta$CellType == cell_type]
    block <- state_abundance[state_ids, , drop = FALSE]
    winners <- apply(block, 2, function(x) state_ids[which.max(x)])
    winner_scores <- apply(block, 2, max)
    rows[[idx]] <- data.frame(
      ID = colnames(block),
      CellType = cell_type,
      State = format_state_label(state_meta$StateLabel[match(winners, state_meta$State)]),
      Abundance = winner_scores,
      stringsAsFactors = FALSE
    )
    idx <- idx + 1L
  }
  do.call(rbind, rows)
}

compute_jaccard_matrix <- function(binary_state_matrix) {
  n_states <- nrow(binary_state_matrix)
  out <- matrix(0, nrow = n_states, ncol = n_states, dimnames = list(rownames(binary_state_matrix), rownames(binary_state_matrix)))
  for (i in seq_len(n_states)) {
    xi <- binary_state_matrix[i, ] > 0
    for (j in i:n_states) {
      xj <- binary_state_matrix[j, ] > 0
      union_count <- sum(xi | xj)
      val <- if (union_count == 0) 0 else sum(xi & xj) / union_count
      out[i, j] <- val
      out[j, i] <- val
    }
  }
  diag(out) <- 1
  out
}

compute_hypergeom_pvals <- function(binary_state_matrix) {
  n_states <- nrow(binary_state_matrix)
  n_samples <- ncol(binary_state_matrix)
  pmat <- matrix(1, nrow = n_states, ncol = n_states, dimnames = list(rownames(binary_state_matrix), rownames(binary_state_matrix)))
  state_counts <- rowSums(binary_state_matrix > 0)
  for (i in seq_len(n_states)) {
    xi <- binary_state_matrix[i, ] > 0
    for (j in i:n_states) {
      xj <- binary_state_matrix[j, ] > 0
      overlap <- sum(xi & xj)
      pval <- phyper(overlap - 1, state_counts[i], n_samples - state_counts[i], state_counts[j], lower.tail = FALSE)
      pmat[i, j] <- pval
      pmat[j, i] <- pval
    }
  }
  diag(pmat) <- 0
  pmat
}

filter_jaccard_by_significance <- function(jaccard, pvals, alpha) {
  filtered <- jaccard
  filtered[pvals > alpha] <- 0
  diag(filtered) <- 1
  filtered
}

choose_ecotype_number <- function(similarity_matrix, min_k, max_k, min_states_per_ecotype) {
  dist_mat <- as.dist(1 - similarity_matrix)
  hc <- hclust(dist_mat, method = "average")
  metrics <- list()
  best_k <- NULL
  best_score <- -Inf
  max_allowed <- min(max_k, nrow(similarity_matrix) - 1)
  if (max_allowed < min_k) {
    stop("Not enough states are available to evaluate the requested ecotype range.")
  }
  for (k in seq.int(min_k, max_allowed)) {
    labels <- cutree(hc, k = k)
    sizes <- table(labels)
    valid_sizes <- all(sizes >= min_states_per_ecotype)
    sil <- silhouette(labels, dist_mat)
    avg_sil <- mean(sil[, "sil_width"])
    metrics[[as.character(k)]] <- data.frame(
      K = k,
      AverageSilhouette = avg_sil,
      MinClusterSize = min(sizes),
      ValidMinStates = valid_sizes,
      stringsAsFactors = FALSE
    )
    score <- if (valid_sizes) avg_sil else -Inf
    if (score > best_score) {
      best_score <- score
      best_k <- k
    }
  }
  if (is.infinite(best_score)) {
    metric_df <- do.call(rbind, metrics)
    best_k <- metric_df$K[which.max(metric_df$AverageSilhouette)]
  }
  list(best_k = best_k, hclust = hc, metrics = do.call(rbind, metrics))
}

discover_ecotypes_from_states <- function(state_abundance, min_k, max_k, alpha, min_states_per_ecotype) {
  dominant <- select_dominant_states(state_abundance)
  jaccard <- compute_jaccard_matrix(dominant$binary)
  pvals <- compute_hypergeom_pvals(dominant$binary)
  jaccard_sig <- filter_jaccard_by_significance(jaccard, pvals, alpha)
  selection <- choose_ecotype_number(jaccard_sig, min_k, max_k, min_states_per_ecotype)
  state_labels <- cutree(selection$hclust, k = selection$best_k)
  state_meta <- dominant$state_meta
  state_assignments <- data.frame(
    CellType = state_meta$CellType,
    State = format_state_label(state_meta$StateLabel),
    Ecotype = format_ecotype_label(paste0("E", state_labels)),
    StateID = rownames(state_abundance),
    stringsAsFactors = FALSE
  )
  list(
    dominant_binary = dominant$binary,
    state_meta = state_meta,
    jaccard = jaccard,
    jaccard_sig = jaccard_sig,
    pvals = pvals,
    state_labels = state_labels,
    state_assignments = state_assignments,
    selection = selection
  )
}

assign_samples_to_ecotypes <- function(state_abundance, state_assignments) {
  ecotypes <- sort(unique(state_assignments$Ecotype))
  ecotype_abundance <- sapply(ecotypes, function(ecotype_name) {
    states <- state_assignments$StateID[state_assignments$Ecotype == ecotype_name]
    colSums(state_abundance[states, , drop = FALSE], na.rm = TRUE)
  })
  if (is.null(dim(ecotype_abundance))) {
    ecotype_abundance <- matrix(ecotype_abundance, ncol = 1, dimnames = list(colnames(state_abundance), ecotypes))
  }
  rownames(ecotype_abundance) <- colnames(state_abundance)
  assignment <- ecotypes[max.col(ecotype_abundance, ties.method = "first")]
  data.frame(
    ID = rownames(ecotype_abundance),
    Ecotype = assignment,
    Abundance = apply(ecotype_abundance, 1, max),
    stringsAsFactors = FALSE
  ) -> sample_assignments
  list(sample_assignments = sample_assignments, ecotype_abundance = ecotype_abundance)
}

write_ecotyper_style_outputs <- function(output_dir, state_abundance, state_assignments, sample_assignments, ecotype_abundance) {
  state_dirs <- split_state_abundance_by_cell_type(state_abundance)
  state_assignment_table <- build_state_assignment_table(state_abundance)
  state_meta <- parse_state_metadata(rownames(state_abundance))

  for (cell_type in names(state_dirs)) {
    cell_dir <- ensure_dir(file.path(output_dir, cell_type))
    cell_states <- rownames(state_dirs[[cell_type]])
    state_labels <- format_state_label(state_meta$StateLabel[match(cell_states, state_meta$State)])
    abundance_df <- data.frame(state_dirs[[cell_type]], check.names = FALSE)
    rownames(abundance_df) <- state_labels
    write.table(abundance_df, file.path(cell_dir, "state_abundances.txt"), sep = "\t", quote = FALSE, col.names = NA)

    assignment_df <- state_assignment_table[state_assignment_table$CellType == cell_type, c("ID", "State"), drop = FALSE]
    fwrite(as.data.table(assignment_df), file.path(cell_dir, "state_assignment.txt"), sep = "\t")
  }

  ecotype_dir <- ensure_dir(file.path(output_dir, "Ecotypes"))
  ecotype_assignment_df <- sample_assignments[, c("ID", "Ecotype"), drop = FALSE]
  fwrite(as.data.table(ecotype_assignment_df), file.path(ecotype_dir, "ecotype_assignment.txt"), sep = "\t")

  ecotype_abundance_df <- as.data.frame(t(ecotype_abundance), check.names = FALSE)
  rownames(ecotype_abundance_df) <- colnames(ecotype_abundance)
  write.table(ecotype_abundance_df, file.path(ecotype_dir, "ecotype_abundance.txt"), sep = "\t", quote = FALSE, col.names = NA)
  fwrite(as.data.table(state_assignments[, c("CellType", "State", "Ecotype")]), file.path(ecotype_dir, "ecotypes.txt"), sep = "\t")
}

derive_markers <- function(expression_matrix, sample_labels, top_n) {
  label_levels <- sort(unique(sample_labels))
  marker_rows <- list()
  for (lbl in label_levels) {
    in_group <- expression_matrix[, sample_labels == lbl, drop = FALSE]
    out_group <- expression_matrix[, sample_labels != lbl, drop = FALSE]
    effect <- rowMeans(in_group, na.rm = TRUE) - rowMeans(out_group, na.rm = TRUE)
    pvals <- vapply(seq_len(nrow(expression_matrix)), function(i) {
      tryCatch(t.test(in_group[i, ], out_group[i, ])$p.value, error = function(e) 1)
    }, numeric(1))
    padj <- p.adjust(pvals, method = "fdr")
    res <- data.frame(
      Gene = rownames(expression_matrix),
      Ecotype = lbl,
      MeanIn = rowMeans(in_group, na.rm = TRUE),
      MeanOut = rowMeans(out_group, na.rm = TRUE),
      EffectSize = effect,
      PValue = pvals,
      FDR = padj,
      stringsAsFactors = FALSE
    )
    res <- res[order(res$EffectSize, -log10(pmax(res$FDR, 1e-300)), decreasing = TRUE), ]
    marker_rows[[lbl]] <- head(res, top_n)
  }
  do.call(rbind, marker_rows)
}

metadata_association_tests <- function(metadata, sample_labels) {
  metadata$Ecotype <- factor(sample_labels)
  test_rows <- list()
  columns <- setdiff(names(metadata), c("SampleID", "Ecotype"))
  idx <- 1L
  for (column in columns) {
    values <- metadata[[column]]
    if (is.numeric(values)) {
      fit <- kruskal.test(values ~ metadata$Ecotype)
      test_rows[[idx]] <- data.frame(
        Variable = column,
        VariableType = "numeric",
        Test = "kruskal",
        Statistic = unname(fit$statistic),
        PValue = fit$p.value,
        stringsAsFactors = FALSE
      )
    } else {
      tbl <- table(values, metadata$Ecotype)
      if (all(dim(tbl) > 1)) {
        fit <- suppressWarnings(chisq.test(tbl))
        test_rows[[idx]] <- data.frame(
          Variable = column,
          VariableType = "categorical",
          Test = "chisq",
          Statistic = unname(fit$statistic),
          PValue = fit$p.value,
          stringsAsFactors = FALSE
        )
      }
    }
    idx <- idx + 1L
  }
  if (length(test_rows) == 0) {
    return(data.frame())
  }
  out <- do.call(rbind, test_rows)
  out$FDR <- p.adjust(out$PValue, method = "fdr")
  out[order(out$PValue), , drop = FALSE]
}

plot_state_pca <- function(state_abundance, state_assignments, out_path) {
  pca <- prcomp(t(state_abundance), center = TRUE, scale. = TRUE)
  df <- data.frame(
    SampleID = rownames(pca$x),
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    Ecotype = state_assignments$Ecotype[match(rownames(pca$x), state_assignments$ID)],
    stringsAsFactors = FALSE
  )
  plt <- ggplot(df, aes(x = PC1, y = PC2, color = Ecotype)) +
    geom_point(size = 3, alpha = 0.9) +
    theme_bw(base_size = 12) +
    labs(title = "Sample PCA from Cell-State Abundances", x = "PC1", y = "PC2")
  ggsave(out_path, plt, width = 8, height = 5.5, dpi = 300)
}

save_similarity_heatmap <- function(similarity_matrix, state_labels, out_path) {
  ecotype_factor <- factor(paste0("E", state_labels), levels = sort(unique(paste0("E", state_labels))))
  palette <- structure(
    brewer.pal(max(3, length(levels(ecotype_factor))), "Set2")[seq_along(levels(ecotype_factor))],
    names = levels(ecotype_factor)
  )
  ha <- HeatmapAnnotation(Ecotype = ecotype_factor, col = list(Ecotype = palette))
  ht <- Heatmap(
    similarity_matrix,
    name = "Jaccard",
    top_annotation = ha,
    show_row_names = FALSE,
    show_column_names = FALSE,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    col = colorRampPalette(c("#f7fbff", "#6baed6", "#08306b"))(100)
  )
  png(out_path, width = 1800, height = 1600, res = 200)
  draw(ht)
  dev.off()
}

plot_ecotype_number_curve <- function(metrics, out_path) {
  plt <- ggplot(metrics, aes(x = K, y = AverageSilhouette)) +
    geom_line(color = "#2c7fb8", linewidth = 1) +
    geom_point(aes(shape = ValidMinStates, color = ValidMinStates), size = 3) +
    scale_color_manual(values = c("FALSE" = "#d95f0e", "TRUE" = "#1b9e77")) +
    theme_bw(base_size = 12) +
    labs(title = "Ecotype Number Selection", x = "Number of ecotypes", y = "Average silhouette")
  ggsave(out_path, plt, width = 7, height = 4.5, dpi = 300)
}

save_state_ecotype_heatmap <- function(state_abundance, state_assignments, out_path) {
  order_states <- order(state_assignments$Ecotype, state_assignments$State)
  ordered_state_ids <- state_assignments$StateID[order_states]
  mat <- state_abundance[ordered_state_ids, , drop = FALSE]
  rownames(mat) <- state_assignments$State[order_states]
  row_split <- factor(state_assignments$Ecotype[order_states], levels = sort(unique(state_assignments$Ecotype)))
  ht <- Heatmap(
    mat,
    name = "Abundance",
    row_split = row_split,
    show_row_names = TRUE,
    show_column_names = FALSE,
    cluster_rows = FALSE,
    cluster_columns = TRUE,
    col = colorRampPalette(c("#f7fcf5", "#74c476", "#00441b"))(100)
  )
  png(out_path, width = 1800, height = 1600, res = 200)
  draw(ht)
  dev.off()
}

save_marker_heatmap <- function(expression_matrix, marker_table, sample_labels, out_path) {
  marker_genes <- unique(marker_table$Gene)
  marker_mat <- expression_matrix[marker_genes, , drop = FALSE]
  column_split <- factor(sample_labels, levels = sort(unique(sample_labels)))
  ht <- Heatmap(
    marker_mat,
    name = "ScaledExpr",
    show_row_names = FALSE,
    show_column_names = FALSE,
    column_split = column_split,
    cluster_columns = TRUE,
    cluster_rows = TRUE,
    col = colorRampPalette(c("#2166ac", "#f7f7f7", "#b2182b"))(100)
  )
  png(out_path, width = 1800, height = 1800, res = 200)
  draw(ht)
  dev.off()
}

plot_ecotype_by_metadata <- function(metadata, sample_labels, column_name, out_path) {
  if (is.null(column_name) || !column_name %in% names(metadata)) {
    return(invisible(NULL))
  }
  df <- data.frame(
    SampleID = metadata$SampleID,
    Ecotype = factor(sample_labels),
    Group = metadata[[column_name]]
  )
  plt <- ggplot(df, aes(x = Group, fill = Ecotype)) +
    geom_bar(position = "fill") +
    theme_bw(base_size = 12) +
    labs(title = paste("Ecotype composition by", column_name), y = "Fraction", x = column_name) +
    scale_y_continuous(labels = function(x) paste0(round(x * 100), "%"))
  ggsave(out_path, plt, width = 8, height = 5, dpi = 300)
}
