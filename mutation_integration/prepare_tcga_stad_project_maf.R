source("mutation_integration/mutation_utils.R")
suppressPackageStartupMessages({
  library(data.table)
})

defaults <- list(
  maf_dir = "C:/Users/difen/Rcode/ecotyper/data/tcga_stad_raw/open_maf_files",
  output_path = "C:/Users/difen/Rcode/ecotyper/data/tcga_stad_raw/TCGA.STAD.open_masked_somatic.maf.tsv.gz"
)
args <- parse_simple_args(defaults)

dir_create_safe(dirname(args$output_path))

maf_files <- list.files(args$maf_dir, pattern = "\\.maf\\.gz$", full.names = TRUE, recursive = TRUE)
if (length(maf_files) == 0) {
  stop(sprintf("No .maf.gz files found under %s", args$maf_dir))
}

maf_list <- lapply(maf_files, function(path) {
  dt <- fread(path, skip = "Hugo_Symbol", fill = TRUE, data.table = TRUE, showProgress = FALSE)
  dt[, source_file := basename(path)]
  dt
})

maf <- rbindlist(maf_list, fill = TRUE, use.names = TRUE)
fwrite(maf, args$output_path, sep = "\t")

summary_lines <- c(
  sprintf("Input MAF files: %d", length(maf_files)),
  sprintf("Combined rows: %d", nrow(maf)),
  sprintf("Combined unique samples: %d", uniqueN(harmonize_tcga_sample_id(maf$Tumor_Sample_Barcode))),
  sprintf("Output: %s", args$output_path)
)
writeLines(summary_lines, sub("\\.tsv\\.gz$", "_summary.txt", args$output_path))
