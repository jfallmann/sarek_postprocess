## Shared helpers for the sarek_postprocess R scripts.
## Adapted from GENOMICS/MAF_Analysis.R and GENOMICS/Intersect_MAF_COSMIC.R,
## generalized to be config-driven instead of hard-coded per cell line.

suppressPackageStartupMessages({
  library(data.table)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

#' Write a data.table to a gzip-compressed .tsv.gz file.
#'
#' data.table's fwrite() cannot write to connections, so this writes a plain
#' temp file next to the target, gzips it into place, and removes the temp.
#' fread() transparently reads .gz, so downstream consumers need no changes.
fwrite_gz <- function(dt, path, ...) {
  tmp <- tempfile(
    pattern = paste0(basename(sub("\\.gz$", "", path)), "."),
    tmpdir = dirname(path)
  )
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
  fwrite(dt, tmp, ...)
  if (requireNamespace("R.utils", quietly = TRUE)) {
    R.utils::gzip(tmp, destname = path, remove = TRUE)
  } else {
    system2("gzip", c("-f", shQuote(tmp)))
    stopifnot(file.exists(paste0(tmp, ".gz")))
    if (file.exists(path)) unlink(path)
    file.rename(paste0(tmp, ".gz"), path)
  }
  invisible(path)
}

#' PASS-filter a maftools MAF object, honouring a per-caller override list
#' (e.g. Strelka indels have no reliable FILTER/QUAL and should be skipped).
pass_filter_maf <- function(maf_obj, caller, qual_override_by_caller = list()) {
  skip <- isFALSE(qual_override_by_caller[[caller]])
  if (skip || !"FILTER" %in% colnames(maf_obj@data)) {
    return(maf_obj)
  }
  tryCatch(
    maftools::subsetMaf(maf_obj, query = "FILTER == 'PASS' | FILTER == '.'"),
    error = function(e) {
      message("PASS filter failed for caller ", caller, ": ", e$message)
      maf_obj
    }
  )
}

#' Apply gnomAD AF / depth / alt-count / VAF / gene-blacklist filters to a
#' MAF data.table (m@data), returning the filtered data.table.
apply_quality_filters <- function(dt, filtering_cfg) {
  keep <- rep(TRUE, nrow(dt))

  if (!is.null(filtering_cfg$exclude_gene_regex) && nzchar(filtering_cfg$exclude_gene_regex) &&
      "Hugo_Symbol" %in% names(dt)) {
    keep <- keep & !grepl(filtering_cfg$exclude_gene_regex, dt$Hugo_Symbol, perl = TRUE)
  }

  if (!is.null(filtering_cfg$min_depth) && "t_depth" %in% names(dt)) {
    keep <- keep & (is.na(dt$t_depth) | dt$t_depth >= filtering_cfg$min_depth)
  }

  if (!is.null(filtering_cfg$min_alt_count) && "t_alt_count" %in% names(dt)) {
    keep <- keep & (is.na(dt$t_alt_count) | dt$t_alt_count >= filtering_cfg$min_alt_count)
  }

  if (!is.null(filtering_cfg$min_vaf) && all(c("t_alt_count", "t_depth") %in% names(dt))) {
    vaf <- ifelse(!is.na(dt$t_depth) & dt$t_depth > 0, dt$t_alt_count / dt$t_depth, NA_real_)
    keep <- keep & (is.na(vaf) | vaf >= filtering_cfg$min_vaf)
  }

  if (!is.null(filtering_cfg$max_gnomad_af)) {
    af_col <- intersect(c("gnomADe_AF", "gnomADg_AF"), names(dt))
    if (length(af_col) > 0) {
      af <- suppressWarnings(as.numeric(dt[[af_col[1]]]))
      keep <- keep & (is.na(af) | af <= filtering_cfg$max_gnomad_af)
    }
  }

  dt[keep]
}

#' Load an optional "known baseline mutations" TSV (e.g. COSMIC Cell Lines
#' Project export for the parental line) and return a data.table keyed on
#' Chromosome/Start_Position/End_Position/Gene, or NULL if not configured.
load_baseline_mutations <- function(tsv_path) {
  if (is.null(tsv_path) || !nzchar(tsv_path)) return(NULL)
  if (!file.exists(tsv_path)) {
    message("Baseline mutations TSV not found at ", tsv_path, "; skipping baseline exclusion")
    return(NULL)
  }
  dt <- fread(tsv_path, sep = "\t", header = TRUE)
  if ("Position" %in% names(dt) && !all(c("Start_Position", "End_Position") %in% names(dt))) {
    dt[, c("Chromosome", "Start_Position", "End_Position") := {
      pos <- gsub("chr", "", Position, ignore.case = TRUE)
      parts <- tstrsplit(pos, "[:.]+", fixed = FALSE)
      list(parts[[1]], as.integer(parts[[2]]), as.integer(parts[[3]]))
    }]
  }
  if (!"Gene" %in% names(dt) && "Hugo_Symbol" %in% names(dt)) setnames(dt, "Hugo_Symbol", "Gene")
  need <- c("Chromosome", "Start_Position", "End_Position", "Gene")
  if (!all(need %in% names(dt))) {
    message("Baseline mutations TSV missing required columns (", paste(need, collapse = ", "), "); skipping")
    return(NULL)
  }
  dt[, Chromosome := gsub("^chr", "", Chromosome, ignore.case = TRUE)]
  dt
}

#' Remove variants already present in the baseline mutation catalog.
exclude_baseline_mutations <- function(dt, baseline_dt) {
  if (is.null(baseline_dt) || nrow(dt) == 0) return(dt)
  if (!all(c("Hugo_Symbol", "Chromosome", "Start_Position", "End_Position") %in% names(dt))) return(dt)
  dtc <- copy(dt)
  dtc[, Chromosome_clean := gsub("^chr", "", Chromosome, ignore.case = TRUE)]
  hit <- dtc[baseline_dt, on = .(Chromosome_clean = Chromosome, Start_Position, End_Position, Hugo_Symbol = Gene),
             nomatch = 0L, which = TRUE]
  if (length(hit) == 0) return(dt)
  dt[-hit]
}

#' Read the gene panel CSV (single "Gene" column) used for driver/resistance
#' prioritization.
read_gene_panel <- function(csv_path) {
  if (is.null(csv_path) || !nzchar(csv_path) || !file.exists(csv_path)) {
    message("Gene panel CSV not found at ", csv_path, "; candidate ranking will not use a panel")
    return(character(0))
  }
  unique(fread(csv_path, header = TRUE)[[1]])
}

impactful_variant_classes <- c(
  "Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins",
  "Splice_Site", "In_Frame_Del", "In_Frame_Ins", "Nonstop_Mutation",
  "Translation_Start_Site", "Start_Codon_Del", "Start_Codon_Ins", "Start_Codon_SNP",
  "De_novo_Start_InFrame", "De_novo_Start_OutOfFrame"
)
