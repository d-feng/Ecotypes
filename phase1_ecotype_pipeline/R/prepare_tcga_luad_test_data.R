suppressPackageStartupMessages({
  library(data.table)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

project_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
raw_dir <- file.path(project_dir, "data", "tcga_luad_test")
out_dir <- file.path(project_dir, "data", "tcga_luad_prepared")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

expr_raw <- fread(file.path(raw_dir, "TCGA-LUAD.star_counts.tsv.gz"))
clin_raw <- fread(file.path(raw_dir, "TCGA-LUAD.clinical.tsv.gz"))
surv_raw <- fread(file.path(raw_dir, "TCGA-LUAD.survival.tsv.gz"))

sample_cols <- setdiff(names(expr_raw), "Ensembl_ID")
primary_samples <- unique(clin_raw[tissue_type.samples == "Tumor" & sample_type.samples == "Primary Tumor", sample])
sample_cols <- intersect(sample_cols, primary_samples)

expr_mat <- as.matrix(expr_raw[, ..sample_cols])
rownames(expr_mat) <- sub("\\..*$", "", expr_raw$Ensembl_ID)
storage.mode(expr_mat) <- "numeric"

symbol_map <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys = rownames(expr_mat),
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

keep <- !is.na(symbol_map) & symbol_map != ""
expr_mat <- expr_mat[keep, , drop = FALSE]
symbols <- symbol_map[keep]

expr_dt <- as.data.table(expr_mat)
expr_dt[, Gene := symbols]
expr_gene <- expr_dt[, lapply(.SD, mean), by = Gene]
setcolorder(expr_gene, c("Gene", sample_cols))

case_meta <- unique(clin_raw[, .(
  SampleID = sample,
  PatientID = submitter_id,
  TumorType = project_id.project,
  Stage = ajcc_pathologic_stage.diagnoses,
  Gender = gender.demographic,
  AgeYears = age_at_index.demographic,
  SmokingYears = years_smoked.exposures,
  PackYears = pack_years_smoked.exposures,
  VitalStatus = vital_status.demographic
)])
case_meta <- case_meta[SampleID %in% sample_cols]

surv_use <- unique(surv_raw[, .(SampleID = sample, OS_time = OS.time, OS = OS)])
meta <- merge(case_meta, surv_use, by = "SampleID", all.x = TRUE)
meta <- meta[match(sample_cols, SampleID)]

score_signature <- function(expr_dt, genes) {
  hits <- intersect(genes, expr_dt$Gene)
  if (length(hits) == 0) {
    return(rep(0, ncol(expr_dt) - 1))
  }
  vals <- as.matrix(expr_dt[Gene %in% hits, -1])
  if (length(hits) == 1) {
    vals <- matrix(vals, nrow = 1)
  }
  colMeans(vals)
}

state_signatures <- list(
  "Malignant|S1" = c("EPCAM", "KRT8", "KRT18", "KRT19", "MUC1"),
  "Malignant|S2" = c("MKI67", "TOP2A", "BIRC5", "UBE2C", "TYMS"),
  "Malignant|S3" = c("VIM", "ZEB1", "ITGA6", "LAMC2", "CXCL8"),
  "Tcell|S1" = c("CD3D", "CD3E", "TRBC1", "IL7R", "LTB"),
  "Tcell|S2" = c("NKG7", "GZMB", "PRF1", "CTSW", "GNLY"),
  "Tcell|S3" = c("PDCD1", "LAG3", "TIGIT", "HAVCR2", "CXCL13"),
  "Myeloid|S1" = c("LST1", "FCER1G", "TYROBP", "AIF1", "C1QB"),
  "Myeloid|S2" = c("S100A8", "S100A9", "CXCL2", "IL1B", "CTSB"),
  "Myeloid|S3" = c("HLA-DRA", "HLA-DPA1", "HLA-DPB1", "CD74", "C1QC"),
  "Fibroblast|S1" = c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM"),
  "Fibroblast|S2" = c("ACTA2", "TAGLN", "MYL9", "TPM2", "DES"),
  "Fibroblast|S3" = c("CXCL14", "APOD", "CFD", "C7", "DPT")
)

state_scores <- lapply(state_signatures, function(genes) score_signature(expr_gene, genes))
state_mat <- do.call(rbind, state_scores)
rownames(state_mat) <- names(state_signatures)
colnames(state_mat) <- sample_cols

state_mat <- t(scale(t(state_mat)))
state_mat[is.na(state_mat)] <- 0
state_mat <- state_mat - apply(state_mat, 1, min)
state_mat <- state_mat + 1e-3
state_mat <- sweep(state_mat, 2, colSums(state_mat), "/")

fwrite(expr_gene, file.path(out_dir, "expression.tsv"), sep = "\t")
fwrite(meta, file.path(out_dir, "metadata.tsv"), sep = "\t")
fwrite(as.data.table(data.frame(State = rownames(state_mat), state_mat, check.names = FALSE)), file.path(out_dir, "cell_state_abundances.tsv"), sep = "\t")

cat("Prepared TCGA LUAD test inputs written to:\n")
cat(file.path(out_dir, "expression.tsv"), "\n")
cat(file.path(out_dir, "metadata.tsv"), "\n")
cat(file.path(out_dir, "cell_state_abundances.tsv"), "\n")
cat(sprintf("Primary tumor samples retained: %d\n", length(sample_cols)))
cat(sprintf("Mapped genes retained: %d\n", nrow(expr_gene)))
