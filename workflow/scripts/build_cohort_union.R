## Merge all per-contrast union MAF TSVs (from consolidate_mafs.R) into a
## single cohort-wide maftools MAF list + merged MAF object, cached as
## MAFlist.rds.gz / UnionMAF.rds.gz. Generalizes GENOMICS/oncoplots_from_maf.R
## (which did the same for one hard-coded cell-line project) to any cohort,
## with an optional sample-subset regex.

suppressPackageStartupMessages({
  library(data.table)
  library(maftools)
})

union_tsvs   <- unlist(snakemake@input[["union_tsvs"]])
contrasts    <- unlist(snakemake@params[["contrasts"]])
subset_regex <- snakemake@params[["sample_subset_regex"]]
out_maflist  <- snakemake@output[["maf_list_rds"]]
out_union    <- snakemake@output[["union_rds"]]

stopifnot(length(union_tsvs) == length(contrasts))

if (!is.null(subset_regex) && nzchar(subset_regex)) {
  keep <- grepl(subset_regex, contrasts)
  union_tsvs <- union_tsvs[keep]
  contrasts  <- contrasts[keep]
}

read_contrast_maf <- function(path, contrast) {
  if (!file.exists(path) || file.info(path)$size == 0) return(NULL)
  dt <- tryCatch(fread(path), error = function(e) NULL)
  if (is.null(dt) || nrow(dt) == 0) return(NULL)
  m <- tryCatch(read.maf(maf = dt, verbose = FALSE, vc_nonSyn = maftools::vc.nonSyn), error = function(e) {
    message("read.maf failed for contrast ", contrast, ": ", e$message)
    NULL
  })
  if (is.null(m)) return(NULL)
  tsb <- unique(m@data$Tumor_Sample_Barcode)
  m@clinical.data <- data.table(Tumor_Sample_Barcode = tsb, Contrast = contrast)
  m
}

maf_list <- Map(read_contrast_maf, union_tsvs, contrasts)
maf_list <- Filter(Negate(is.null), maf_list)

dir.create(dirname(out_maflist), showWarnings = FALSE, recursive = TRUE)

if (length(maf_list) == 0) {
  message("No non-empty contrast MAFs found; writing placeholder RDS files")
  saveRDS(list(), gzfile(out_maflist))
  saveRDS(NULL, gzfile(out_union))
  quit(save = "no", status = 0)
}

maf_union <- if (length(maf_list) == 1) maf_list[[1]] else merge_mafs(maf = maf_list, verbose = FALSE)
all_clin <- rbindlist(lapply(maf_list, function(m) m@clinical.data), use.names = TRUE, fill = TRUE)
all_clin <- unique(all_clin, by = "Tumor_Sample_Barcode")
maf_union@clinical.data <- all_clin

saveRDS(maf_union, gzfile(out_union))
saveRDS(maf_list, gzfile(out_maflist))

message("Cohort union built from ", length(maf_list), " contrasts")
