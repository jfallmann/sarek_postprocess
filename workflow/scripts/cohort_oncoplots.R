## Cohort-wide oncoplots and summary dashboard from the cached union MAF.
## Generalizes GENOMICS/oncoplots_from_maf.R plotting section (the caching
## section moved to build_cohort_union.R so it can be a separate, reusable
## Snakemake target).

suppressPackageStartupMessages({
  library(maftools)
})

maf_union_rds <- snakemake@input[["union_rds"]]
out_dir       <- snakemake@params[["out_dir"]]
done_marker   <- snakemake@output[["done"]]

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

file.create(done_marker)
message("Cohort oncoplots written to ", out_dir)
