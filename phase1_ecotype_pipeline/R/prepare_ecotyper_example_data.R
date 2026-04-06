suppressPackageStartupMessages({
  library(data.table)
})

project_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
raw_dir <- file.path(project_dir, "data", "ecotyper_example_raw")
out_dir <- file.path(project_dir, "data", "ecotyper_example_prepared")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

expr <- fread(file.path(raw_dir, "bulk_lung_data.txt"))
anno <- fread(file.path(raw_dir, "bulk_lung_annotation.txt"))
fractions <- fread(file.path(raw_dir, "bulk_fractions_example.txt"))

sample_cols <- setdiff(names(expr), "Gene")
sample_ids <- Reduce(intersect, list(sample_cols, anno$ID, fractions$Mixture))
expr_use <- expr[, c("Gene", sample_ids), with = FALSE]
anno_use <- unique(anno[ID %in% sample_ids])
fractions_use <- fractions[match(sample_ids, Mixture)]

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
  "Epithelial.cells|S1" = c("EPCAM", "KRT8", "KRT18", "KRT19", "MUC1"),
  "Epithelial.cells|S2" = c("MKI67", "TOP2A", "BIRC5", "UBE2C", "TYMS"),
  "Epithelial.cells|S3" = c("VIM", "LAMC2", "ITGA6", "CXCL8", "SPP1"),
  "Fibroblasts|S1" = c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM"),
  "Fibroblasts|S2" = c("ACTA2", "TAGLN", "MYL9", "TPM2", "DES"),
  "Fibroblasts|S3" = c("CXCL14", "APOD", "CFD", "C7", "DPT"),
  "Endothelial.cells|S1" = c("KDR", "ESAM", "EMCN", "PECAM1", "VWF"),
  "Endothelial.cells|S2" = c("SELE", "ACKR1", "PLVAP", "RAMP2", "RGCC"),
  "Endothelial.cells|S3" = c("CXCL12", "CD34", "FLT1", "ENG", "EPAS1"),
  "B.cells|S1" = c("MS4A1", "CD79A", "CD79B", "HLA-DRA", "CD74"),
  "B.cells|S2" = c("MZB1", "SDC1", "JCHAIN", "XBP1", "DERL3"),
  "B.cells|S3" = c("BANK1", "CD22", "FCRL1", "PAX5", "HVCN1"),
  "CD4.T.cells|S1" = c("IL7R", "LTB", "MALAT1", "LTB", "MALAT1"),
  "CD4.T.cells|S2" = c("PDCD1", "LAG3", "TIGIT", "CXCL13", "HAVCR2"),
  "CD4.T.cells|S3" = c("FOXP3", "IL2RA", "CTLA4", "TNFRSF4", "TIGIT"),
  "CD8.T.cells|S1" = c("NKG7", "PRF1", "GZMB", "CTSW", "GNLY"),
  "CD8.T.cells|S2" = c("PDCD1", "LAG3", "TIGIT", "HAVCR2", "TOX"),
  "CD8.T.cells|S3" = c("IL7R", "LTB", "CXCR3", "CCL5", "KLRD1"),
  "Dendritic.cells|S1" = c("FCER1A", "CLEC10A", "HLA-DRA", "CD1C", "CST3"),
  "Dendritic.cells|S2" = c("LILRA4", "GZMB", "IRF7", "TCF4", "JCHAIN"),
  "Dendritic.cells|S3" = c("CCR7", "CCL19", "LAMP3", "IDO1", "FSCN1"),
  "Mast.cells|S1" = c("TPSAB1", "TPSB2", "CPA3", "KIT", "HPGDS"),
  "Mast.cells|S2" = c("HDC", "GATA2", "MS4A2", "CLC", "RGS13"),
  "Mast.cells|S3" = c("IL1RL1", "AREG", "RGS1", "CCL2", "CCL4"),
  "Monocytes.and.Macrophages|S1" = c("LST1", "FCER1G", "TYROBP", "AIF1", "C1QB"),
  "Monocytes.and.Macrophages|S2" = c("S100A8", "S100A9", "IL1B", "CXCL2", "CTSB"),
  "Monocytes.and.Macrophages|S3" = c("HLA-DRA", "CD74", "C1QC", "APOE", "MSR1"),
  "NK.cells|S1" = c("NKG7", "GNLY", "PRF1", "FCGR3A", "CTSW"),
  "NK.cells|S2" = c("XCL1", "XCL2", "KLRD1", "TRAC", "CCL5"),
  "NK.cells|S3" = c("KLRC1", "TYROBP", "IFITM1", "SPON2", "FGFBP2"),
  "PCs|S1" = c("MZB1", "JCHAIN", "SDC1", "XBP1", "TNFRSF17"),
  "PCs|S2" = c("IGHG1", "IGHG3", "IGKC", "IGLC2", "CD27"),
  "PCs|S3" = c("DERL3", "SSR4", "SEC11C", "FKBP11", "TXNDC5"),
  "PMNs|S1" = c("CXCR2", "FCGR3B", "CSF3R", "S100A8", "S100A9"),
  "PMNs|S2" = c("MMP8", "MMP9", "RETN", "LCN2", "HP"),
  "PMNs|S3" = c("IFITM2", "IFITM3", "OASL", "ISG15", "RSAD2")
)

fraction_map <- c(
  "Epithelial.cells" = "Epithelial.cells",
  "Fibroblasts" = "Fibroblasts",
  "Endothelial.cells" = "Endothelial.cells",
  "B.cells" = "B.cells",
  "CD4.T.cells" = "CD4.T.cells",
  "CD8.T.cells" = "CD8.T.cells",
  "Dendritic.cells" = "Dendritic.cells",
  "Mast.cells" = "Mast.cells",
  "Monocytes.and.Macrophages" = "Monocytes.and.Macrophages",
  "NK.cells" = "NK.cells",
  "PCs" = "PCs",
  "PMNs" = "PMNs"
)

state_scores <- lapply(names(state_signatures), function(state_name) {
  cell_type <- sub("\\|.*$", "", state_name)
  signature_score <- score_signature(expr_use, state_signatures[[state_name]])
  signature_score <- as.numeric(scale(signature_score))
  signature_score[is.na(signature_score)] <- 0
  fraction_values <- fractions_use[[fraction_map[[cell_type]]]]
  pmax(signature_score + 2, 0) * pmax(fraction_values, 1e-4)
})

state_mat <- do.call(rbind, state_scores)
rownames(state_mat) <- names(state_signatures)
colnames(state_mat) <- sample_ids
state_mat <- state_mat + 1e-3
state_mat <- sweep(state_mat, 2, colSums(state_mat), "/")

metadata <- data.table(
  SampleID = anno_use$ID[match(sample_ids, anno_use$ID)],
  Tissue = anno_use$Tissue[match(sample_ids, anno_use$ID)],
  Histology = anno_use$Histology[match(sample_ids, anno_use$ID)],
  Type = anno_use$Type[match(sample_ids, anno_use$ID)],
  OS_Time = anno_use$OS_Time[match(sample_ids, anno_use$ID)],
  OS_Status = anno_use$OS_Status[match(sample_ids, anno_use$ID)]
)

fwrite(expr_use, file.path(out_dir, "expression.tsv"), sep = "\t")
fwrite(metadata, file.path(out_dir, "metadata.tsv"), sep = "\t")
fwrite(as.data.table(data.frame(State = rownames(state_mat), state_mat, check.names = FALSE)), file.path(out_dir, "cell_state_abundances.tsv"), sep = "\t")

cat("Prepared EcoTyper example inputs written to:\n")
cat(file.path(out_dir, "expression.tsv"), "\n")
cat(file.path(out_dir, "metadata.tsv"), "\n")
cat(file.path(out_dir, "cell_state_abundances.tsv"), "\n")
cat(sprintf("Samples retained: %d\n", length(sample_ids)))
cat(sprintf("Genes retained: %d\n", nrow(expr_use)))
