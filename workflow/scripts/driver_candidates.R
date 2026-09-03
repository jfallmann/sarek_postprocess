## Candidate driver/resistance gene ranking per contrast and cohort-wide.
##
## With only a handful of clones (no cohort-scale statistical power), this
## does NOT rely on a single p-value cutoff. Instead it combines, per gene:
##   - presence in the consensus MAF (>= N callers agreeing, see
##     consolidate_mafs.R) for the contrast
##   - recurrence across independent contrasts (clone lineages)
##   - membership in the curated driver/resistance gene panel
##   - VEP IMPACT / consequence severity already present in the MAF
##   - optional mafCompare/forestPlot vs the pooled "sensitive" contrasts,
##     reported as exploratory (unadjusted) evidence, not a hard filter
##   - optional corroborating CNVkit copy-number evidence (gain/loss calls
##     for the same gene in the same or other contrasts)
##
## Adapted from GENOMICS/MAF_Analysis.R filtering and
## GENOMICS/MAF_filter_high_confidence.R recurrence logic.

suppressPackageStartupMessages({
  library(data.table)
  library(maftools)
})

source(snakemake@params[["common_r"]])

consensus_paths <- unlist(snakemake@input[["consensus_mafs"]])
contrasts       <- unlist(snakemake@params[["contrasts"]])
sensitive_contrasts <- unlist(snakemake@params[["sensitive_contrasts"]])
filtering       <- snakemake@params[["filtering"]]
gene_panel_csv  <- snakemake@params[["gene_panel_csv"]]
baseline_tsv    <- snakemake@params[["baseline_mutations_tsv"]]
cnv_calls_tsv   <- snakemake@params[["cnv_calls_tsv"]]
out_per_contrast_dir <- snakemake@params[["out_per_contrast_dir"]]
out_cohort_ranked     <- snakemake@output[["cohort_ranked"]]

dir.create(out_per_contrast_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(out_cohort_ranked), showWarnings = FALSE, recursive = TRUE)

gene_panel <- read_gene_panel(gene_panel_csv)
baseline_dt <- load_baseline_mutations(baseline_tsv)

read_consensus <- function(path, contrast) {
  if (!file.exists(path) || file.info(path)$size == 0) return(NULL)
  dt <- tryCatch(fread(path), error = function(e) NULL)
  if (is.null(dt) || nrow(dt) == 0) return(NULL)
  dt <- exclude_baseline_mutations(dt, baseline_dt)
  if (nrow(dt) == 0) return(NULL)
  if ("Variant_Classification" %in% names(dt)) {
    dt <- dt[Variant_Classification %in% impactful_variant_classes]
  }
  if (nrow(dt) == 0) return(NULL)
  dt[, Contrast := contrast]
  dt
}

dt_list <- Map(read_consensus, consensus_paths, contrasts)
dt_list <- Filter(Negate(is.null), dt_list)

if (length(dt_list) == 0) {
  message("No consensus variants remain across any contrast; writing empty cohort ranking")
  fwrite(data.table(), out_cohort_ranked, sep = "\t")
  quit(save = "no", status = 0)
}

all_dt <- rbindlist(dt_list, use.names = TRUE, fill = TRUE)

## Per-contrast candidate tables -------------------------------------------
for (contrast in unique(all_dt$Contrast)) {
  sub <- all_dt[Contrast == contrast]
  sub[, In_Gene_Panel := Hugo_Symbol %in% gene_panel]
  fwrite(sub[order(-In_Gene_Panel, Hugo_Symbol)],
         file.path(out_per_contrast_dir, paste0(contrast, ".candidate_drivers.tsv")),
         sep = "\t", quote = FALSE)
}

## Cohort-wide recurrence ranking -------------------------------------------
recur <- all_dt[, .(
  n_contrasts = uniqueN(Contrast),
  contrasts   = paste(sort(unique(Contrast)), collapse = ";"),
  n_variants  = .N
), by = Hugo_Symbol]
recur[, In_Gene_Panel := Hugo_Symbol %in% gene_panel]
recur[, Score := n_contrasts + 0.5 * In_Gene_Panel]

## Optional corroborating CNVkit evidence ------------------------------------
if (!is.null(cnv_calls_tsv) && nzchar(cnv_calls_tsv) && file.exists(cnv_calls_tsv) &&
    file.info(cnv_calls_tsv)$size > 0) {
  cnv_calls <- fread(cnv_calls_tsv)
  if (all(c("Hugo_Symbol", "Call", "Contrast") %in% names(cnv_calls))) {
    cnv_summary <- cnv_calls[Call != "neutral", .(
      CNV_Gain_Contrasts = paste(sort(unique(Contrast[Call == "gain"])), collapse = ";"),
      CNV_Loss_Contrasts = paste(sort(unique(Contrast[Call == "loss"])), collapse = ";"),
      n_CNV_gain = uniqueN(Contrast[Call == "gain"]),
      n_CNV_loss = uniqueN(Contrast[Call == "loss"])
    ), by = Hugo_Symbol]
    recur <- merge(recur, cnv_summary, by = "Hugo_Symbol", all.x = TRUE)
    for (col in c("n_CNV_gain", "n_CNV_loss")) recur[is.na(get(col)), (col) := 0L]
    for (col in c("CNV_Gain_Contrasts", "CNV_Loss_Contrasts")) recur[is.na(get(col)), (col) := ""]
    recur[, Score := Score + 0.25 * (n_CNV_gain > 0 | n_CNV_loss > 0)]
    message("Merged CNVkit evidence for ", nrow(cnv_summary), " genes with a gain/loss call")
  } else {
    message("CNV calls TSV missing expected columns; skipping CNV evidence merge")
  }
} else {
  message("No CNV calls TSV configured/found; skipping CNV evidence merge")
}

setorder(recur, -Score, -n_contrasts, Hugo_Symbol)

fwrite(recur, out_cohort_ranked, sep = "\t", quote = FALSE)
message("Cohort driver ranking (", nrow(recur), " genes) written to ", out_cohort_ranked)

## Optional exploratory mafCompare vs pooled sensitive contrasts -------------
if (length(sensitive_contrasts) > 0) {
  resistant_contrasts <- setdiff(unique(all_dt$Contrast), sensitive_contrasts)
  sens_dt <- all_dt[Contrast %in% sensitive_contrasts]
  for (rc in resistant_contrasts) {
    res_dt <- all_dt[Contrast == rc]
    if (nrow(res_dt) == 0 || nrow(sens_dt) == 0) next
    m_res  <- tryCatch(read.maf(maf = res_dt, verbose = FALSE), error = function(e) NULL)
    m_sens <- tryCatch(read.maf(maf = sens_dt, verbose = FALSE), error = function(e) NULL)
    if (is.null(m_res) || is.null(m_sens)) next
    cmp <- tryCatch(
      mafCompare(m1 = m_res, m2 = m_sens, m1Name = rc, m2Name = "Sensitive_pool", minMut = 1),
      error = function(e) {
        message("mafCompare failed for ", rc, ": ", e$message)
        NULL
      }
    )
    if (is.null(cmp)) next
    fwrite(cmp$results, file.path(out_per_contrast_dir, paste0(rc, ".vs_sensitive_pool.mafCompare.tsv")),
           sep = "\t", quote = FALSE)
    pdf_path <- file.path(out_per_contrast_dir, paste0(rc, ".vs_sensitive_pool.forestplot.pdf"))
    tryCatch({
      pdf(pdf_path, width = 8, height = 8)
      print(forestPlot(mafCompareRes = cmp, pVal = 1, color = c("maroon", "royalblue")))
      dev.off()
    }, error = function(e) {
      message("forestPlot failed for ", rc, ": ", e$message)
      if (dev.cur() > 1) dev.off()
    })
  }
} else {
  message("No sensitive-pool contrasts configured; skipping exploratory mafCompare step")
}
