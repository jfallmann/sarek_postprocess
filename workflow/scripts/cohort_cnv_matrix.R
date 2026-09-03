## Build a cohort-wide gene x contrast CNVkit log2-ratio matrix and heatmap,
## restricted to a gene panel (defaults to the same driver/resistance panel
## used for SNV/indel ranking) to keep the heatmap readable.

suppressPackageStartupMessages({
  library(data.table)
  library(ComplexHeatmap)
  library(circlize)
})

source(snakemake@params[["common_r"]])

gene_tsvs   <- unlist(snakemake@input[["gene_tsvs"]])
contrasts   <- unlist(snakemake@params[["contrasts"]])
gene_panel_csv <- snakemake@params[["gene_panel_csv"]]
min_gain    <- snakemake@params[["min_log2_gain"]]
min_loss    <- snakemake@params[["min_log2_loss"]]
out_matrix  <- snakemake@output[["matrix_tsv"]]
out_calls   <- snakemake@output[["calls_tsv"]]
out_heatmap <- snakemake@output[["heatmap_pdf"]]

dir.create(dirname(out_matrix), showWarnings = FALSE, recursive = TRUE)

read_one <- function(path) {
  if (!file.exists(path) || file.info(path)$size == 0) return(NULL)
  dt <- tryCatch(fread(path), error = function(e) NULL)
  if (is.null(dt) || nrow(dt) == 0) return(NULL)
  dt
}

dt_list <- lapply(gene_tsvs, read_one)
dt_list <- Filter(Negate(is.null), dt_list)

if (length(dt_list) == 0) {
  message("No CNVkit gene summaries available; writing empty CNV matrix/calls")
  fwrite(data.table(), out_matrix, sep = "\t")
  fwrite(data.table(), out_calls, sep = "\t")
  pdf(out_heatmap); plot.new(); text(0.5, 0.5, "No CNV data available"); dev.off()
  quit(save = "no", status = 0)
}

all_dt <- rbindlist(dt_list, use.names = TRUE, fill = TRUE)

## Calls table (used by driver_candidates.R) - one row per gene x contrast --
all_dt[, Call := fifelse(log2_weighted >= min_gain, "gain",
                   fifelse(log2_weighted <= min_loss, "loss", "neutral"))]
fwrite(all_dt, out_calls, sep = "\t", quote = FALSE)

## Wide matrix + heatmap, restricted to a gene panel -------------------------
panel <- read_gene_panel(if (nzchar(gene_panel_csv)) gene_panel_csv else snakemake@params[["default_gene_panel_csv"]])

heat_dt <- all_dt[Hugo_Symbol %in% panel]
if (nrow(heat_dt) == 0) {
  message("No gene-panel genes found in CNVkit output; writing full-gene matrix instead")
  heat_dt <- all_dt
}

mat <- dcast(heat_dt, Hugo_Symbol ~ Contrast, value.var = "log2_weighted", fun.aggregate = mean)
gene_names <- mat$Hugo_Symbol
mat[, Hugo_Symbol := NULL]
mat_m <- as.matrix(mat)
rownames(mat_m) <- gene_names

fwrite(data.table(Hugo_Symbol = gene_names, mat), out_matrix, sep = "\t", quote = FALSE)

if (nrow(mat_m) >= 1 && ncol(mat_m) >= 1) {
  col_fun <- colorRamp2(c(min(mat_m, na.rm = TRUE), 0, max(mat_m, na.rm = TRUE)),
                         c("blue", "white", "red"))
  pdf(out_heatmap, width = max(8, ncol(mat_m) * 0.6), height = max(6, nrow(mat_m) * 0.25))
  print(Heatmap(mat_m, name = "log2 ratio", col = col_fun,
                na_col = "grey90", cluster_rows = TRUE, cluster_columns = TRUE,
                row_names_gp = grid::gpar(fontsize = 8),
                column_names_gp = grid::gpar(fontsize = 8)))
  dev.off()
} else {
  pdf(out_heatmap); plot.new(); text(0.5, 0.5, "No CNV data available for heatmap"); dev.off()
}

message("Cohort CNV matrix (", nrow(mat_m), " genes x ", ncol(mat_m), " contrasts) written to ", out_matrix)
