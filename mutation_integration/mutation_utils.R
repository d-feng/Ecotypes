suppressPackageStartupMessages({
  library(data.table)
})

parse_simple_args <- function(defaults = list()) {
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

dir_create_safe <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

harmonize_tcga_sample_id <- function(x) {
  substr(gsub("-", ".", x, fixed = TRUE), 1, 16)
}

read_ecotyper_matrix <- function(path) {
  lines <- readLines(path, warn = FALSE)
  if (length(lines) < 2) {
    stop(sprintf("File has no data rows: %s", path))
  }

  header <- strsplit(lines[[1]], "\t", fixed = TRUE)[[1]]
  body <- strsplit(lines[-1], "\t", fixed = TRUE)

  row_names <- vapply(body, function(x) gsub('"', "", x[[1]], fixed = TRUE), character(1))
  value_matrix <- do.call(rbind, lapply(body, function(x) as.numeric(x[-1])))

  colnames(value_matrix) <- gsub('"', "", header, fixed = TRUE)
  rownames(value_matrix) <- row_names

  as.data.frame(value_matrix, check.names = FALSE)
}

dominant_ecotype_table <- function(ecotype_abundance_path) {
  abundance <- read_ecotyper_matrix(ecotype_abundance_path)
  sample_ids <- colnames(abundance)
  ecotypes <- rownames(abundance)

  assignment <- data.table(
    ID = sample_ids,
    Ecotype = vapply(sample_ids, function(sample_id) {
      vals <- abundance[, sample_id]
      ecotypes[[which.max(vals)]]
    }, character(1))
  )

  for (ecotype in ecotypes) {
    assignment[[paste0("abundance_", ecotype)]] <- as.numeric(abundance[ecotype, sample_ids])
  }

  assignment
}

get_stad_driver_genes <- function() {
  c(
    "TP53", "ARID1A", "PIK3CA", "CDH1", "KRAS", "ERBB2", "SMAD4", "RHOA",
    "KMT2D", "RNF43", "APC", "FAT4", "PTEN", "FBXW7", "CTNNB1", "FGFR2",
    "FGFR3", "MET", "POLE", "MLH1", "MSH2", "MSH6", "PMS2", "JAK1",
    "JAK2", "B2M", "HLA.A", "HLA.B", "HLA.C"
  )
}

get_stad_pathways <- function() {
  list(
    RTK_RAS = c("ERBB2", "EGFR", "MET", "KRAS", "NRAS", "BRAF", "FGFR2", "FGFR3"),
    PI3K = c("PIK3CA", "PIK3R1", "PTEN", "AKT1", "AKT2", "AKT3"),
    TP53 = c("TP53", "MDM2", "MDM4", "ATM", "CHEK2"),
    TGF_BETA = c("SMAD4", "TGFBR1", "TGFBR2", "ACVR2A"),
    WNT = c("APC", "CTNNB1", "RNF43", "FAT4"),
    CHROMATIN = c("ARID1A", "KMT2D", "KMT2C", "CREBBP", "EP300"),
    DNA_REPAIR = c("POLE", "MLH1", "MSH2", "MSH6", "PMS2", "ATM", "BRCA1", "BRCA2"),
    FGFR = c("FGFR1", "FGFR2", "FGFR3", "FGF3", "FGF4", "FGF19"),
    ANTIGEN_PRESENTATION = c("B2M", "HLA.A", "HLA.B", "HLA.C", "TAP1", "TAP2"),
    IMMUNE_SIGNALING = c("JAK1", "JAK2", "STAT3", "STAT1", "IFNGR1", "IFNGR2")
  )
}

get_nonsilent_classes <- function() {
  c(
    "Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site",
    "Nonsense_Mutation", "Nonstop_Mutation", "In_Frame_Del", "In_Frame_Ins",
    "Missense_Mutation", "Splice_Region", "Start_Codon_Del", "Start_Codon_Ins",
    "Start_Codon_SNP", "Stop_Codon_Del", "Stop_Codon_Ins"
  )
}

read_maf_table <- function(path) {
  maf <- fread(path, fill = TRUE, data.table = TRUE, showProgress = FALSE)
  required <- c("Tumor_Sample_Barcode", "Hugo_Symbol", "Variant_Classification")
  missing_cols <- setdiff(required, names(maf))
  if (length(missing_cols) > 0) {
    stop(sprintf("MAF is missing required columns: %s", paste(missing_cols, collapse = ", ")))
  }
  maf[, sample_id := harmonize_tcga_sample_id(Tumor_Sample_Barcode)]
  maf[, Hugo_Symbol := gsub("-", ".", Hugo_Symbol, fixed = TRUE)]
  maf
}

build_binary_feature_table <- function(maf, sample_ids, gene_list, pathway_map) {
  nonsilent_classes <- get_nonsilent_classes()
  maf_nonsilent <- maf[Variant_Classification %in% nonsilent_classes]

  out <- data.table(ID = sample_ids)

  gene_hits <- unique(maf_nonsilent[Hugo_Symbol %in% gene_list, .(sample_id, Hugo_Symbol)])
  for (gene in gene_list) {
    hit_samples <- gene_hits[Hugo_Symbol == gene, unique(sample_id)]
    out[[paste0("gene_", gene)]] <- as.integer(sample_ids %in% hit_samples)
  }

  for (pathway_name in names(pathway_map)) {
    hit_samples <- unique(maf_nonsilent[Hugo_Symbol %in% pathway_map[[pathway_name]], sample_id])
    out[[paste0("pathway_", pathway_name)]] <- as.integer(sample_ids %in% hit_samples)
  }

  lof_classes <- c("Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Nonsense_Mutation")
  msi_genes <- pathway_map$DNA_REPAIR
  antigen_genes <- pathway_map$ANTIGEN_PRESENTATION

  out[, flag_lof_high := as.integer(ID %in% unique(maf[Variant_Classification %in% lof_classes, sample_id]))]
  out[, flag_dna_repair_hit := as.integer(ID %in% unique(maf_nonsilent[Hugo_Symbol %in% msi_genes, sample_id]))]
  out[, flag_antigen_presentation_hit := as.integer(ID %in% unique(maf_nonsilent[Hugo_Symbol %in% antigen_genes, sample_id]))]

  out
}

build_continuous_feature_table <- function(maf, sample_ids) {
  nonsilent_classes <- get_nonsilent_classes()
  maf_nonsilent <- maf[Variant_Classification %in% nonsilent_classes]

  counts <- maf_nonsilent[, .(
    nonsilent_mutation_count = .N,
    mutated_gene_count = uniqueN(Hugo_Symbol)
  ), by = sample_id]

  missense_counts <- maf_nonsilent[, .(
    missense_mutation_count = sum(Variant_Classification == "Missense_Mutation")
  ), by = sample_id]

  out <- data.table(ID = sample_ids)
  out <- merge(out, counts, by.x = "ID", by.y = "sample_id", all.x = TRUE)
  out <- merge(out, missense_counts, by.x = "ID", by.y = "sample_id", all.x = TRUE)

  numeric_cols <- setdiff(names(out), "ID")
  for (col in numeric_cols) {
    set(out, which(is.na(out[[col]])), col, 0)
  }

  out[, log_nonsilent_mutation_count := log1p(nonsilent_mutation_count)]
  out[, log_mutated_gene_count := log1p(mutated_gene_count)]
  out
}

safe_fisher <- function(in_group, out_group) {
  tbl <- matrix(
    c(sum(in_group), length(in_group) - sum(in_group), sum(out_group), length(out_group) - sum(out_group)),
    nrow = 2,
    byrow = TRUE
  )
  test <- suppressWarnings(fisher.test(tbl))
  list(
    p_value = unname(test$p.value),
    odds_ratio = ifelse(is.null(test$estimate), NA_real_, unname(test$estimate))
  )
}

safe_wilcox <- function(x_in, x_out) {
  if (length(unique(c(x_in, x_out))) <= 1) {
    return(list(p_value = 1, effect_size = 0))
  }
  test <- suppressWarnings(wilcox.test(x_in, x_out, exact = FALSE))
  list(
    p_value = unname(test$p.value),
    effect_size = median(x_in, na.rm = TRUE) - median(x_out, na.rm = TRUE)
  )
}
