## Build per-contrast union and consensus MAFs from the per-caller MAFs
## produced by convert_to_maf.sh.
##
## Union: all PASS variants from any caller, for sensitivity-first cohort
##        plots (oncoplot / mutational burden).
## Consensus: variants confirmed by >= consensus_min_callers callers on
##        (Chromosome, Start_Position, End_Position, Tumor_Sample_Barcode,
##        Hugo_Symbol) - used as the higher-confidence input to driver
##        candidate ranking. Adapted from the multi-caller overlap logic in
##        GENOMICS/Intersect_MAF_COSMIC.R.
##
## Called by Snakemake's `script:` directive - uses the `snakemake` S4 object
## for input/output/params/wildcards.

suppressPackageStartupMessages({
  library(data.table)
  library(maftools)
})

source(snakemake@params[["common_r"]])

maf_paths   <- unlist(snakemake@input[["mafs"]])
callers     <- unlist(snakemake@params[["callers"]])
contrast    <- snakemake@wildcards[["contrast"]]
filtering   <- snakemake@params[["filtering"]]
out_union     <- snakemake@output[["union_maf"]]
out_consensus <- snakemake@output[["consensus_maf"]]

stopifnot(length(maf_paths) == length(callers))

read_one <- function(path, caller) {
  if (!file.exists(path) || file.info(path)$size == 0) {
    message("Skipping empty/missing MAF: ", path)
    return(NULL)
  }
  m <- tryCatch(read.maf(maf = path, verbose = FALSE), error = function(e) {
    message("read.maf failed for ", path, ": ", e$message)
    NULL
  })
  if (is.null(m)) return(NULL)
  m <- pass_filter_maf(m, caller, filtering$qual_override_by_caller)
  if (nrow(m@data) == 0) return(NULL)
  dt <- as.data.table(m@data)
  dt[, Caller := caller]
  dt[, Contrast := contrast]
  dt
}

dt_list <- Map(read_one, maf_paths, callers)
dt_list <- Filter(Negate(is.null), dt_list)

dir.create(dirname(out_union), showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(out_consensus), showWarnings = FALSE, recursive = TRUE)

if (length(dt_list) == 0) {
  message("No variants remain for contrast ", contrast, " after PASS filtering; writing empty outputs")
  fwrite_gz(data.table(), out_union, sep = "\t")
  fwrite_gz(data.table(), out_consensus, sep = "\t")
  quit(save = "no", status = 0)
}

union_dt <- rbindlist(dt_list, use.names = TRUE, fill = TRUE)
union_dt <- apply_quality_filters(union_dt, filtering)
fwrite_gz(union_dt, out_union, sep = "\t", quote = FALSE)

key_cols <- c("Chromosome", "Start_Position", "End_Position", "Tumor_Sample_Barcode", "Hugo_Symbol")
key_cols <- intersect(key_cols, names(union_dt))

if (length(key_cols) == length(c("Chromosome", "Start_Position", "End_Position", "Tumor_Sample_Barcode", "Hugo_Symbol"))) {
  caller_counts <- union_dt[, .(n_callers = uniqueN(Caller)), by = key_cols]
  min_callers <- filtering$consensus_min_callers %||% 2L
  consensus_keys <- caller_counts[n_callers >= min_callers]
  consensus_dt <- union_dt[consensus_keys, on = key_cols, mult = "first"]
} else {
  message("Cannot build consensus (missing key columns); falling back to union for contrast ", contrast)
  consensus_dt <- union_dt
}

fwrite_gz(consensus_dt, out_consensus, sep = "\t", quote = FALSE)

message("Contrast ", contrast, ": union=", nrow(union_dt), " consensus=", nrow(consensus_dt))
