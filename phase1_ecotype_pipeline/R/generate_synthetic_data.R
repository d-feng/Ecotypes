suppressPackageStartupMessages({
  library(data.table)
})

set.seed(11)

root_dir <- normalizePath(file.path(getwd()), winslash = "/", mustWork = TRUE)
data_dir <- file.path(root_dir, "data")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

n_genes <- 700
n_samples <- 72
sample_ids <- paste0("Sample_", sprintf("%03d", seq_len(n_samples)))
truth_ecotypes <- rep(c("E1", "E2", "E3"), each = n_samples / 3)

state_names <- c(
  "Malignant|S1", "Malignant|S2", "Malignant|S3",
  "Tcell|S1", "Tcell|S2", "Tcell|S3",
  "Myeloid|S1", "Myeloid|S2", "Myeloid|S3",
  "Fibroblast|S1", "Fibroblast|S2", "Fibroblast|S3"
)

state_abundance <- matrix(runif(length(state_names) * n_samples, min = 0.01, max = 0.08), nrow = length(state_names), ncol = n_samples)
rownames(state_abundance) <- state_names
colnames(state_abundance) <- sample_ids

ecotype_programs <- list(
  E1 = c("Malignant|S1", "Tcell|S2", "Myeloid|S3", "Fibroblast|S1"),
  E2 = c("Malignant|S2", "Tcell|S1", "Myeloid|S1", "Fibroblast|S3"),
  E3 = c("Malignant|S3", "Tcell|S3", "Myeloid|S2", "Fibroblast|S2")
)

for (ecotype_name in names(ecotype_programs)) {
  sample_idx <- which(truth_ecotypes == ecotype_name)
  states <- ecotype_programs[[ecotype_name]]
  state_abundance[states, sample_idx] <- state_abundance[states, sample_idx] + 0.35
}

state_abundance <- sweep(state_abundance, 2, colSums(state_abundance), "/")

gene_ids <- paste0("Gene", sprintf("%04d", seq_len(n_genes)))
expression <- matrix(rpois(n_genes * n_samples, lambda = 20), nrow = n_genes, ncol = n_samples)
rownames(expression) <- gene_ids
colnames(expression) <- sample_ids

state_to_gene_block <- split(seq_len(480), rep(state_names, each = 40))
for (state in names(state_to_gene_block)) {
  genes <- state_to_gene_block[[state]]
  expression[genes, ] <- expression[genes, ] + matrix(
    round(120 * rep(state_abundance[state, ], each = length(genes))),
    nrow = length(genes),
    byrow = FALSE
  )
}

tumor_type <- rep(c("LUAD", "BRCA", "HNSC"), length.out = n_samples)
stage <- sample(c("I", "II", "III"), n_samples, replace = TRUE, prob = c(0.35, 0.35, 0.30))
purity <- round(runif(n_samples, min = 0.45, max = 0.90), 2)

metadata <- data.table(
  SampleID = sample_ids,
  TruthEcotype = truth_ecotypes,
  TumorType = tumor_type,
  Stage = stage,
  Purity = purity
)

expression_dt <- data.table(Gene = rownames(expression), expression)
setnames(expression_dt, old = names(expression_dt)[-1], new = sample_ids)

state_dt <- data.table(State = rownames(state_abundance), state_abundance)
setnames(state_dt, old = names(state_dt)[-1], new = sample_ids)

fwrite(expression_dt, file.path(data_dir, "expression.tsv"), sep = "\t")
fwrite(metadata, file.path(data_dir, "metadata.tsv"), sep = "\t")
fwrite(state_dt, file.path(data_dir, "cell_state_abundances.tsv"), sep = "\t")

cat("Synthetic data written to:\n")
cat(file.path(data_dir, "expression.tsv"), "\n")
cat(file.path(data_dir, "metadata.tsv"), "\n")
cat(file.path(data_dir, "cell_state_abundances.tsv"), "\n")
