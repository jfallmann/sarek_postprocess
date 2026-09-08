## Cohort-wide oncoplots and summary dashboard from the cached union MAF.
## Generalizes GENOMICS/oncoplots_from_maf.R plotting section (the caching
## section moved to build_cohort_union.R so it can be a separate, reusable
## Snakemake target).

suppressPackageStartupMessages({
  library(maftools)
})

source(snakemake@params[["common_r"]])

maf_union_rds  <- snakemake@input[["union_rds"]]
out_dir        <- snakemake@params[["out_dir"]]
gene_panel_csv <- snakemake@params[["gene_panel_csv"]]
done_marker    <- snakemake@output[["done"]]

## tryCatch wrapper around a single plotting call, mirroring
## GENOMICS/MAF_Analysis.R's safe_plot() - one gene/plot failing (e.g. no
## amino-acid annotation for lollipopPlot) must not abort the whole rule.
safe_plot <- function(expr) {
  tryCatch(expr, error = function(e) message("Plot failed: ", e$message))
}

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

maf_union <- readRDS(maf_union_rds)

if (is.null(maf_union) || nrow(maf_union@data) == 0) {
  message("No cohort-wide variants available; skipping oncoplots")
  file.create(file.path(out_dir, "NO_VARIANTS.txt"))
  file.create(done_marker)
  quit(save = "no", status = 0)
}

maf_pass_all <- tryCatch(
  maftools::subsetMaf(maf_union, query = "FILTER == 'PASS' | FILTER == '.'"),
  error = function(e) maf_union
)

gene_panel <- read_gene_panel(gene_panel_csv)

## Contrasts fall into two, background-incompatible groups that must never be
## plotted together: bare tumor-only calls ("vs_reference") and calls made
## against the matched Bulk_sensitive contrast ("vs_contrast", i.e. any
## Contrast containing "_vs_"). Mixing them makes the cohort overviews
## meaningless (different variant-calling background/sensitivity), and with
## every contrast combined the oncoplots/heatmaps become unreadably wide.
all_clin <- maf_pass_all@clinical.data
groups <- list(
  vs_reference = unique(all_clin$Tumor_Sample_Barcode[!grepl("_vs_", all_clin$Contrast)]),
  vs_contrast  = unique(all_clin$Tumor_Sample_Barcode[grepl("_vs_", all_clin$Contrast)])
)

## One self-contained set of cohort plots per group. Ported from
## GENOMICS/MAF_Analysis.R's oncoplot/summary/interactions blocks, now
## parameterised over the tsb subset and with sample-name labels forced on
## and plot dimensions scaled to the number of samples so long sample names
## are not clipped.
run_group_oncoplots <- function(maf_pass, group_out_dir, group_label) {
  dir.create(group_out_dir, showWarnings = FALSE, recursive = TRUE)

  n_genes   <- length(unique(maf_pass@gene.summary$Hugo_Symbol))
  n_samples <- length(unique(maf_pass@data$Tumor_Sample_Barcode))
  max_tsb_len <- max(nchar(unique(maf_pass@data$Tumor_Sample_Barcode)), 1)

  if (n_genes >= 2) {
    oc_width <- max(21, n_samples * 0.6 + max_tsb_len * 0.15)
    pdf(file.path(group_out_dir, "SummaryOncoplot.pdf"), width = oc_width, height = 21)
    safe_plot(print(oncoplot(maf = maf_pass, top = min(20, n_genes),
                   clinicalFeatures = "Contrast", sortByAnnotation = TRUE,
                   showTumorSampleBarcodes = TRUE, removeNonMutated = TRUE)))
    dev.off()
  } else {
    message("[", group_label, "] Not enough mutated genes (", n_genes, ") to draw cohort oncoplot")
  }

  pdf(file.path(group_out_dir, "MAFSummary.pdf"), width = 21, height = 21)
  safe_plot(plotmafSummary(maf = maf_pass, addStat = "median", dashboard = TRUE))
  dev.off()

  ## Gene-panel-restricted oncoplot, ported from GENOMICS/MAF_Analysis.R's
  ## "04_oncoplot_goi" (goi = curated driver/resistance gene list).
  present_panel <- intersect(gene_panel, maf_pass@gene.summary$Hugo_Symbol)
  if (length(present_panel) >= 2) {
    panel_width <- max(14, n_samples * 0.5 + max_tsb_len * 0.15)
    pdf(file.path(group_out_dir, "DriverPanelOncoplot.pdf"), width = panel_width, height = 10)
    safe_plot(print(oncoplot(maf = maf_pass, genes = present_panel,
                              clinicalFeatures = "Contrast", sortByAnnotation = TRUE,
                              showTumorSampleBarcodes = TRUE, removeNonMutated = FALSE)))
    dev.off()
  } else {
    message("[", group_label, "] Skipping driver-panel oncoplot (", length(present_panel),
            " panel genes present, need >=2)")
  }

  ## Lollipop plots per gene-panel gene present in this group, ported from
  ## GENOMICS/MAF_Analysis.R (AA-column fallback: HGVSp_Short ->
  ## Protein_position -> Amino_acids).
  lollipop_dir <- file.path(group_out_dir, "lollipop")
  dir.create(lollipop_dir, showWarnings = FALSE, recursive = TRUE)
  for (g in present_panel) {
    gene_df <- maf_pass@data[Hugo_Symbol == g]
    aa_col <- NULL
    for (candidate in c("HGVSp_Short", "Protein_position", "Amino_acids")) {
      if (candidate %in% names(gene_df) && any(!is.na(gene_df[[candidate]]) & gene_df[[candidate]] != "")) {
        aa_col <- candidate
        break
      }
    }
    if (is.null(aa_col)) {
      message("[", group_label, "] Skipping lollipop for ", g, " (no amino-acid/position annotations)")
      next
    }
    pdf(file.path(lollipop_dir, paste0(g, ".lollipop.pdf")), width = 9, height = 6)
    safe_plot(print(lollipopPlot(maf = maf_pass, gene = g, AACol = aa_col,
                                  labelPos = "all", showDomainLabel = FALSE)))
    dev.off()
  }

  ## Mutual exclusivity / co-occurrence, ported from GENOMICS/MAF_Analysis.R's
  ## "05_interactions" block.
  if (n_genes >= 2 && n_samples >= 2) {
    pdf(file.path(group_out_dir, "SomaticInteractions.pdf"), width = 10, height = 9)
    safe_plot(print(somaticInteractions(maf = maf_pass, top = 25, pvalue = c(0.05, 0.1))))
    dev.off()
  } else {
    message("[", group_label, "] Skipping somaticInteractions (need >=2 genes and >=2 samples)")
  }
}

for (group_label in names(groups)) {
  tsb_group <- groups[[group_label]]
  group_out_dir <- file.path(out_dir, group_label)
  if (length(tsb_group) == 0) {
    message("No samples in group '", group_label, "'; skipping")
    dir.create(group_out_dir, showWarnings = FALSE, recursive = TRUE)
    file.create(file.path(group_out_dir, "NO_SAMPLES.txt"))
    next
  }
  maf_group <- tryCatch(
    maftools::subsetMaf(maf_pass_all, tsb = tsb_group),
    error = function(e) NULL
  )
  if (is.null(maf_group) || nrow(maf_group@data) == 0) {
    message("No variants in group '", group_label, "'; skipping")
    dir.create(group_out_dir, showWarnings = FALSE, recursive = TRUE)
    file.create(file.path(group_out_dir, "NO_VARIANTS.txt"))
    next
  }
  run_group_oncoplots(maf_group, group_out_dir, group_label)
}

file.create(done_marker)
message("Cohort oncoplots written to ", out_dir)
