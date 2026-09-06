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

maf_pass <- tryCatch(
  maftools::subsetMaf(maf_union, query = "FILTER == 'PASS' | FILTER == '.'"),
  error = function(e) maf_union
)

n_genes <- length(unique(maf_pass@gene.summary$Hugo_Symbol))

if (n_genes >= 2) {
  pdf(file.path(out_dir, "SummaryOncoplot.pdf"), width = 21, height = 21)
  print(oncoplot(maf = maf_pass, top = min(20, n_genes),
                 clinicalFeatures = "Contrast", sortByAnnotation = TRUE,
                 removeNonMutated = TRUE))
  dev.off()
} else {
  message("Not enough mutated genes (", n_genes, ") to draw cohort oncoplot")
}

pdf(file.path(out_dir, "MAFSummary.pdf"), width = 21, height = 21)
plotmafSummary(maf = maf_pass, addStat = "median", dashboard = TRUE)
dev.off()

## Gene-panel-restricted oncoplot, ported from GENOMICS/MAF_Analysis.R's
## "04_oncoplot_goi" (goi = curated driver/resistance gene list).
gene_panel <- read_gene_panel(gene_panel_csv)
present_panel <- intersect(gene_panel, maf_pass@gene.summary$Hugo_Symbol)
if (length(present_panel) >= 2) {
  pdf(file.path(out_dir, "DriverPanelOncoplot.pdf"), width = 14, height = 10)
  safe_plot(print(oncoplot(maf = maf_pass, genes = present_panel,
                            clinicalFeatures = "Contrast", sortByAnnotation = TRUE,
                            showTumorSampleBarcodes = TRUE, removeNonMutated = FALSE)))
  dev.off()
} else {
  message("Skipping driver-panel oncoplot (", length(present_panel), " panel genes present, need >=2)")
}

## Lollipop plots per gene-panel gene present in the cohort, ported from
## GENOMICS/MAF_Analysis.R (AA-column fallback: HGVSp_Short -> Protein_position
## -> Amino_acids).
lollipop_dir <- file.path(out_dir, "lollipop")
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
    message("Skipping lollipop for ", g, " (no amino-acid/position annotations)")
    next
  }
  pdf(file.path(lollipop_dir, paste0(g, ".lollipop.pdf")), width = 9, height = 6)
  safe_plot(print(lollipopPlot(maf = maf_pass, gene = g, AACol = aa_col,
                                labelPos = "all", showDomainLabel = FALSE)))
  dev.off()
}

## Mutual exclusivity / co-occurrence, ported from GENOMICS/MAF_Analysis.R's
## "05_interactions" block.
n_samples <- length(unique(maf_pass@data$Tumor_Sample_Barcode))
if (n_genes >= 2 && n_samples >= 2) {
  pdf(file.path(out_dir, "SomaticInteractions.pdf"), width = 10, height = 9)
  safe_plot(print(somaticInteractions(maf = maf_pass, top = 25, pvalue = c(0.05, 0.1))))
  dev.off()
} else {
  message("Skipping somaticInteractions (need >=2 genes and >=2 samples)")
}

file.create(done_marker)
message("Cohort oncoplots written to ", out_dir)
