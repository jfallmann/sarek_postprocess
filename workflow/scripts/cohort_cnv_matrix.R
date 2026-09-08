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
out_heatmap_vs_reference <- snakemake@output[["heatmap_vs_reference_pdf"]]
out_heatmap_vs_contrast  <- snakemake@output[["heatmap_vs_contrast_pdf"]]

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
  fwrite_gz(data.table(), out_matrix, sep = "\t")
  fwrite_gz(data.table(), out_calls, sep = "\t")
  for (p in c(out_heatmap_vs_reference, out_heatmap_vs_contrast)) {
    pdf(p); plot.new(); text(0.5, 0.5, "No CNV data available"); dev.off()
  }
  quit(save = "no", status = 0)
}

all_dt <- rbindlist(dt_list, use.names = TRUE, fill = TRUE)

## Calls table (used by driver_candidates.R) - one row per gene x contrast,
## kept cohort-wide (all contrast types) since driver evidence merging
## wants every CNV hit regardless of contrast type.
all_dt[, Call := fifelse(log2_weighted >= min_gain, "gain",
                   fifelse(log2_weighted <= min_loss, "loss", "neutral"))]
fwrite_gz(all_dt, out_calls, sep = "\t", quote = FALSE)

panel <- read_gene_panel(if (nzchar(gene_panel_csv)) gene_panel_csv else snakemake@params[["default_gene_panel_csv"]])

## One heatmap per contrast group instead of one cohort-wide plot:
## "vs_reference" (bare tumor-only contrasts, called against a generic/
## pooled reference) and "vs_contrast" (matched tumor/normal
## "*_vs_Bulk_sensitive" contrasts) have different copy-number backgrounds,
## are not meaningfully comparable side by side, and mixing them was also
## what made the previous single heatmap too large to open/read.
plot_group_heatmap <- function(dt, out_path, group_label) {
  if (nrow(dt) == 0) {
    pdf(out_path); plot.new(); text(0.5, 0.5, paste0("No CNV data available for group '", group_label, "'")); dev.off()
    return(invisible())
  }
  heat_dt <- dt[Hugo_Symbol %in% panel]
  if (nrow(heat_dt) == 0) {
    message("No gene-panel genes found in CNVkit output for group '", group_label, "'; using full-gene matrix instead")
    heat_dt <- dt
  }
  mat <- dcast(heat_dt, Hugo_Symbol ~ Contrast, value.var = "log2_weighted", fun.aggregate = mean)
  gene_names <- mat$Hugo_Symbol
  mat[, Hugo_Symbol := NULL]
  mat_m <- as.matrix(mat)
  rownames(mat_m) <- gene_names
  if (nrow(mat_m) < 1 || ncol(mat_m) < 1) {
    pdf(out_path); plot.new(); text(0.5, 0.5, paste0("No CNV data available for heatmap (group '", group_label, "')")); dev.off()
    return(invisible())
  }
  col_fun <- colorRamp2(c(min(mat_m, na.rm = TRUE), 0, max(mat_m, na.rm = TRUE)),
                         c("blue", "white", "red"))
  ## Reserve extra width/height for the row/column name labels themselves
  ## (longest gene symbol, longest contrast name) - a fixed per-cell size
  ## alone left long contrast names (e.g. "BC139_resistant_vs_Bulk_sensitive")
  ## clipped at the page edge, since ComplexHeatmap does not expand the PDF
  ## canvas to fit rotated column labels on its own.
  longest_col_name <- max(nchar(colnames(mat_m)))
  longest_row_name <- max(nchar(rownames(mat_m)))
  pdf(out_path,
      width = max(8, ncol(mat_m) * 0.6) + longest_row_name * 0.08,
      height = max(6, nrow(mat_m) * 0.25) + longest_col_name * 0.12)
  print(Heatmap(mat_m, name = "log2 ratio", col = col_fun,
                na_col = "grey90", cluster_rows = TRUE, cluster_columns = TRUE,
                column_title = paste0("CNV heatmap (", group_label, ")"),
                row_names_gp = grid::gpar(fontsize = 8),
                column_names_gp = grid::gpar(fontsize = 8),
                column_names_rot = 45,
                row_names_max_width = unit(longest_row_name * 0.09, "inches"),
                column_names_max_height = unit(longest_col_name * 0.12, "inches")))
  dev.off()
  message("CNV heatmap (group=", group_label, "; ", nrow(mat_m), " genes x ", ncol(mat_m), " contrasts) written to ", out_path)
}

is_vs_contrast <- grepl("_vs_", all_dt$Contrast)
plot_group_heatmap(all_dt[!is_vs_contrast], out_heatmap_vs_reference, "vs_reference")
plot_group_heatmap(all_dt[is_vs_contrast], out_heatmap_vs_contrast, "vs_contrast")

full_mat <- dcast(all_dt, Hugo_Symbol ~ Contrast, value.var = "log2_weighted", fun.aggregate = mean)
fwrite_gz(full_mat, out_matrix, sep = "\t", quote = FALSE)

message("Cohort CNV matrix (", uniqueN(all_dt$Hugo_Symbol), " genes x ", uniqueN(all_dt$Contrast), " contrasts) written to ", out_matrix)
