## Lowest-priority actionability cross-reference: join the cohort-ranked
## driver candidate list against optional static, offline gene/variant
## tables (Cancer Gene Census, OncoKB cancer gene list, CIViC variant
## summary). Every table is optional and skipped with a message if its path
## is empty or missing - this is intentionally a simple join, not a new
## reporting framework.

suppressPackageStartupMessages({
  library(data.table)
})

source(snakemake@params[["common_r"]])

cohort_ranked_tsv <- snakemake@input[["cohort_ranked"]]
cgc_tsv    <- snakemake@params[["cancer_gene_census_tsv"]]
oncokb_tsv <- snakemake@params[["oncokb_gene_list_tsv"]]
civic_tsv  <- snakemake@params[["civic_variant_summary_tsv"]]
out_tsv    <- snakemake@output[["report"]]

dir.create(dirname(out_tsv), showWarnings = FALSE, recursive = TRUE)

ranked <- fread(cohort_ranked_tsv)
if (nrow(ranked) == 0) {
  message("Cohort ranked driver list is empty; writing empty actionability report")
  fwrite_gz(ranked, out_tsv, sep = "\t")
  quit(save = "no", status = 0)
}

load_gene_set <- function(path, gene_col_candidates) {
  if (is.null(path) || !nzchar(path) || !file.exists(path)) {
    message("Skipping optional table (not configured or missing): ", path)
    return(character(0))
  }
  dt <- fread(path)
  gene_col <- intersect(gene_col_candidates, names(dt))
  if (length(gene_col) == 0) {
    message("No recognizable gene column in ", path, "; skipping")
    return(character(0))
  }
  unique(dt[[gene_col[1]]])
}

cgc_genes    <- load_gene_set(cgc_tsv,    c("Gene Symbol", "Hugo_Symbol", "Gene"))
oncokb_genes <- load_gene_set(oncokb_tsv, c("Hugo Symbol", "Hugo_Symbol", "Gene"))
civic_genes  <- load_gene_set(civic_tsv,  c("gene", "Gene", "Hugo_Symbol"))

ranked[, In_CGC    := Hugo_Symbol %in% cgc_genes]
ranked[, In_OncoKB := Hugo_Symbol %in% oncokb_genes]
ranked[, In_CIViC   := Hugo_Symbol %in% civic_genes]
ranked[, Actionability_Sources := In_CGC + In_OncoKB + In_CIViC]
setorder(ranked, -Actionability_Sources, -Score)

fwrite_gz(ranked, out_tsv, sep = "\t", quote = FALSE)
message("Actionability report (", nrow(ranked), " genes) written to ", out_tsv)
