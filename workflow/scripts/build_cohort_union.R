## Merge all per-contrast union MAF TSVs (from consolidate_mafs.R) into a
## single cohort-wide maftools MAF list + merged MAF object, cached as
## MAFlist.rds.gz / UnionMAF.rds.gz. Generalizes GENOMICS/oncoplots_from_maf.R
## (which did the same for one hard-coded cell-line project) to any cohort,
## with an optional sample-subset regex.

suppressPackageStartupMessages({
  library(data.table)
  library(maftools)
})

source(snakemake@params[["common_r"]])

union_tsvs   <- unlist(snakemake@input[["union_tsvs"]])
contrasts    <- unlist(snakemake@params[["contrasts"]])
vcset        <- snakemake@wildcards[["vcset"]]
filtering    <- snakemake@params[["filtering"]]
subset_regex <- snakemake@params[["sample_subset_regex"]]
out_maflist  <- snakemake@output[["maf_list_rds"]]
out_union    <- snakemake@output[["union_rds"]]

stopifnot(length(union_tsvs) == length(contrasts))

## The per-contrast union TSVs were already filtered to this vcset's
## Variant_Classification set by consolidate_mafs.R; re-apply the same set
## here (rather than maftools' own stringent default) so read.maf() doesn't
## silently re-drop the very classes (e.g. Splice_Region/UTR/Intron) the
## "custom" track was built to keep.
vc_nonsyn <- resolve_vc_nonsyn(vcset, filtering)

if (!is.null(subset_regex) && nzchar(subset_regex)) {
  keep <- grepl(subset_regex, contrasts)
  union_tsvs <- union_tsvs[keep]
  contrasts  <- contrasts[keep]
}

read_contrast_maf <- function(path, contrast) {
  if (!file.exists(path) || file.info(path)$size == 0) return(NULL)
  dt <- tryCatch(fread(path), error = function(e) NULL)
  if (is.null(dt) || nrow(dt) == 0) return(NULL)
  m <- tryCatch(read.maf(maf = dt, verbose = FALSE, vc_nonSyn = vc_nonsyn), error = function(e) {
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

## Pull the (tiny) per-contrast clinical rows out before merge_mafs() runs,
## so nothing extra needs to stay attached to maf_list afterwards.
all_clin <- rbindlist(lapply(maf_list, function(m) m@clinical.data), use.names = TRUE, fill = TRUE)
all_clin <- unique(all_clin, by = "Tumor_Sample_Barcode")

maf_union <- if (length(maf_list) == 1) maf_list[[1]] else merge_mafs(maf = maf_list, verbose = FALSE)
maf_union@clinical.data <- all_clin

## Write + drop the two large cohort-wide objects one at a time instead of
## keeping both fully resident while both gzip writes happen back to back -
## on a WGS "custom" (permissive) union track this pair is the single
## biggest peak-memory point in the pipeline.
saveRDS(maf_union, gzfile(out_union))
rm(maf_union); gc()

saveRDS(maf_list, gzfile(out_maflist))
rm(maf_list); gc()

message("Cohort union built from ", length(contrasts), " contrasts")
