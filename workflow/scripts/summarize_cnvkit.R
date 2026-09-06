## Gene-level summary of a single CNVkit .cns/.call.cns file for one contrast.
##
## CNVkit segment files (.cns from `cnvkit.py segment`, .call.cns from
## `cnvkit.py call`) carry chromosome/start/end/gene/log2[/cn/depth/probes/
## weight] columns, where "gene" can be a comma-separated list of genes
## covered by that segment (or "-"/"Background" for off-target bins). This
## expands multi-gene segments to one row per gene and collapses to one
## row per gene via a probe-count-weighted mean log2 ratio (falling back to
## an unweighted mean if no weight/probes column is present).

suppressPackageStartupMessages({
  library(data.table)
})

source(snakemake@params[["common_r"]])

cns_path   <- snakemake@input[["cns"]]
contrast   <- snakemake@wildcards[["contrast"]]
out_tsv    <- snakemake@output[["gene_tsv"]]
min_weight <- snakemake@params[["min_weight"]] %||% 0.5
min_probes <- snakemake@params[["min_probes"]] %||% 0

dir.create(dirname(out_tsv), showWarnings = FALSE, recursive = TRUE)

if (!file.exists(cns_path) || file.info(cns_path)$size == 0) {
  message("CNVkit segment file missing/empty for ", contrast, ": ", cns_path)
  fwrite_gz(data.table(), out_tsv, sep = "\t")
  quit(save = "no", status = 0)
}

cns <- fread(cns_path)

if (!"gene" %in% names(cns) || nrow(cns) == 0) {
  message("No 'gene' column or no rows in ", cns_path, "; cannot summarize CNVkit output for ", contrast)
  fwrite_gz(data.table(), out_tsv, sep = "\t")
  quit(save = "no", status = 0)
}

weight_col <- intersect(c("weight", "probes"), names(cns))
cns[, `:=`(.weight = if (length(weight_col) > 0) get(weight_col[1]) else 1)]

## Drop low-confidence segments before collapsing to gene level (mirrors
## GENOMICS/CNV/CNV_from_WES.R's weight > 0.5 quality gate), so background/
## antitarget-adjacent bins can't pull log2_weighted toward a false
## gain/loss call. Only applied when the corresponding column exists.
n_before <- nrow(cns)
if ("weight" %in% names(cns)) cns <- cns[is.na(weight) | weight >= min_weight]
if ("probes" %in% names(cns) && min_probes > 0) cns <- cns[is.na(probes) | probes >= min_probes]
if (nrow(cns) < n_before) {
  message(nrow(cns), "/", n_before, " segments retained for ", contrast,
          " after weight/probes quality filter (min_weight=", min_weight, ", min_probes=", min_probes, ")")
}

if (nrow(cns) == 0) {
  message("No segments remain after quality filtering for ", contrast)
  fwrite_gz(data.table(), out_tsv, sep = "\t")
  quit(save = "no", status = 0)
}

expanded <- cns[, .(Hugo_Symbol = trimws(unlist(strsplit(gene, "[,;]")))), by = seq_len(nrow(cns))]
expanded <- cbind(expanded, cns[expanded$seq_len, .(log2, .weight, chromosome, start, end,
                                                     cn = if ("cn" %in% names(cns)) cn else NA_integer_)])
expanded <- expanded[!Hugo_Symbol %in% c("", "-", "Background", "Antitarget")]

if (nrow(expanded) == 0) {
  message("No gene-annotated segments remain for ", contrast)
  fwrite_gz(data.table(), out_tsv, sep = "\t")
  quit(save = "no", status = 0)
}

gene_dt <- expanded[, .(
  log2_mean     = mean(log2, na.rm = TRUE),
  log2_weighted = weighted.mean(log2, w = .weight, na.rm = TRUE),
  cn_mode       = if (all(is.na(cn))) NA_integer_ else as.integer(names(sort(table(cn), decreasing = TRUE))[1]),
  n_segments    = .N
), by = Hugo_Symbol]

gene_dt[, Contrast := contrast]
setorder(gene_dt, Hugo_Symbol)

fwrite_gz(gene_dt, out_tsv, sep = "\t", quote = FALSE)
message("CNVkit gene summary for ", contrast, ": ", nrow(gene_dt), " genes -> ", out_tsv)
