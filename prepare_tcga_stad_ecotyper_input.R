suppressPackageStartupMessages({
  library(data.table)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

project_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
raw_dir <- file.path(project_dir, "data", "tcga_stad_raw")
out_dir <- file.path(project_dir, "data", "tcga_stad_prepared")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

expr_raw <- fread(file.path(raw_dir, "TCGA-STAD.star_counts.tsv.gz"))
clin_raw <- fread(file.path(raw_dir, "TCGA-STAD.clinical.tsv.gz"))
surv_raw <- fread(file.path(raw_dir, "TCGA-STAD.survival.tsv.gz"))

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

safe_ids <- function(x) gsub("-", ".", x, fixed = TRUE)
sample_cols_safe <- safe_ids(sample_cols)
setnames(expr_gene, old = setdiff(names(expr_gene), "Gene"), new = sample_cols_safe)
setcolorder(expr_gene, c("Gene", sample_cols_safe))

case_meta <- unique(clin_raw[, .(
  SampleID = sample,
  Tissue = tissue_type.samples,
  Histology = project_id.project,
  Type = sample_type.samples,
  Stage = ajcc_pathologic_stage.diagnoses,
  PrimaryDiagnosis = primary_diagnosis.diagnoses,
  Gender = gender.demographic,
  AgeYears = age_at_index.demographic,
  VitalStatus = vital_status.demographic
)])
case_meta <- case_meta[SampleID %in% sample_cols]

surv_use <- unique(surv_raw[, .(SampleID = sample, OS_Time = OS.time, OS_Status = OS)])
anno <- merge(case_meta, surv_use, by = "SampleID", all.x = TRUE)
anno <- anno[match(sample_cols, SampleID)]
anno[, ID := safe_ids(SampleID)]
anno[, SampleID := NULL]
setcolorder(anno, c("ID", setdiff(names(anno), "ID")))

fwrite(expr_gene, file.path(out_dir, "bulk_stad_data.txt"), sep = "\t")
fwrite(anno, file.path(out_dir, "bulk_stad_annotation.txt"), sep = "\t")

cat("Prepared TCGA STAD EcoTyper inputs written to:\n")
cat(file.path(out_dir, "bulk_stad_data.txt"), "\n")
cat(file.path(out_dir, "bulk_stad_annotation.txt"), "\n")
cat(sprintf("Primary tumor samples retained: %d\n", length(sample_cols_safe)))
cat(sprintf("Mapped genes retained: %d\n", nrow(expr_gene)))
