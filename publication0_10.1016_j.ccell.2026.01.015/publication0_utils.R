suppressPackageStartupMessages({
  library(data.table)
})

publication0_parse_args <- function(defaults = list()) {
  args <- commandArgs(trailingOnly = TRUE)
  parsed <- defaults
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) {
      stop(sprintf("Unexpected positional argument: %s", key))
    }
    key <- substring(key, 3)
    if (i == length(args) || startsWith(args[[i + 1]], "--")) {
      parsed[[key]] <- TRUE
      i <- i + 1
    } else {
      parsed[[key]] <- args[[i + 1]]
      i <- i + 2
    }
  }
  parsed
}

publication0_dir_create <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

publication0_read_expr <- function(path) {
  expr <- fread(path, data.table = FALSE)
  rownames(expr) <- expr[[1]]
  expr[[1]] <- NULL
  expr <- as.matrix(expr)
  mode(expr) <- "numeric"
  expr
}

publication0_safe_gsva <- function(expr, signature_list, method = "gsva") {
  if (!requireNamespace("GSVA", quietly = TRUE)) {
    stop("Package 'GSVA' is required for signature scoring.")
  }
  if (length(signature_list) == 0) {
    stop("No signatures available for GSVA scoring.")
  }

  legacy_try <- try(
    GSVA::gsva(expr = expr, gset.idx.list = signature_list, method = method),
    silent = TRUE
  )
  if (!inherits(legacy_try, "try-error")) {
    return(legacy_try)
  }

  if (exists("gsvaParam", where = asNamespace("GSVA"), mode = "function")) {
    param <- GSVA::gsvaParam(exprData = expr, geneSets = signature_list)
    return(GSVA::gsva(param))
  }

  stop("Unable to run GSVA with the installed package version.")
}

publication0_pick_group_column <- function(metadata) {
  candidates <- c("MOS", "Group", "Subtype", "Ecotype", "Cluster")
  present <- intersect(candidates, colnames(metadata))
  if (length(present) == 0) {
    stop("No grouping column found. Expected one of: MOS, Group, Subtype, Ecotype, Cluster")
  }
  present[[1]]
}

publication0_one_vs_rest <- function(expr, groups, min_logfc = log2(2), max_fdr = 0.05) {
  stopifnot(ncol(expr) == length(groups))
  group_levels <- unique(groups)
  res <- list()

  for (g in group_levels) {
    in_group <- groups == g
    out_group <- groups != g
    pvals <- apply(expr, 1, function(x) {
      suppressWarnings(wilcox.test(x[in_group], x[out_group], exact = FALSE)$p.value)
    })
    mean_in <- rowMeans(expr[, in_group, drop = FALSE], na.rm = TRUE)
    mean_out <- rowMeans(expr[, out_group, drop = FALSE], na.rm = TRUE)
    logfc <- log2(mean_in + 1) - log2(mean_out + 1)
    tab <- data.table(
      gene = rownames(expr),
      group = g,
      mean_in = mean_in,
      mean_out = mean_out,
      logFC = logfc,
      p_value = pvals,
      FDR = p.adjust(pvals, method = "fdr")
    )
    res[[g]] <- tab[logFC > min_logfc & FDR < max_fdr][order(-logFC, FDR)]
  }

  res
}

publication0_build_top_expression_labels <- function(expr, groups) {
  group_levels <- unique(groups)
  means <- sapply(group_levels, function(g) rowMeans(expr[, groups == g, drop = FALSE], na.rm = TRUE))
  if (is.null(dim(means))) {
    means <- matrix(means, ncol = 1)
    colnames(means) <- group_levels
    rownames(means) <- rownames(expr)
  }
  top_group <- colnames(means)[max.col(means, ties.method = "first")]
  data.table(gene = rownames(expr), top_group = top_group)
}
