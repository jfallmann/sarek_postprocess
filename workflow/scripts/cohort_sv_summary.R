## Cohort-wide summary of Manta structural variants: burden per contrast by
## SV type, gene-level hit calls (for corroborating driver_candidates.R,
## mirroring how CNVkit gain/loss evidence is merged in there already), and
## a cohort recurrence ranking restricted/annotated against the driver gene
## panel. There is no equivalent step in the WES-era GENOMICS scripts (SVs
## were not analyzed there); this is new for WGS.
##
## Per-contrast input: results/sv/<contrast>.manta_sv.tsv.gz, written by
## summarize_manta_sv.R (one row per Manta record: Chromosome, Start, End,
## SV_Type, Filter, Gene_Annotation - the raw ";"-joined ANN/CSQ strings).

suppressPackageStartupMessages({
  library(data.table)
})

source(snakemake@params[["common_r"]])

sv_tsvs         <- unlist(snakemake@input[["sv_tsvs"]])
contrasts       <- unlist(snakemake@params[["contrasts"]])
pass_only       <- isTRUE(snakemake@params[["pass_only"]])
gene_field_idx  <- as.integer(snakemake@params[["gene_field_index"]])
gene_panel_csv  <- snakemake@params[["gene_panel_csv"]]
large_sv_min_bp <- as.numeric(snakemake@params[["large_sv_min_bp"]] %||% 1e6)
out_burden      <- snakemake@output[["burden_tsv"]]
out_gene_calls  <- snakemake@output[["gene_calls_tsv"]]
out_recurrent   <- snakemake@output[["recurrent_tsv"]]
out_barplot     <- snakemake@output[["burden_barplot"]]
out_large_sv    <- snakemake@output[["large_sv_tsv"]]
out_large_plot  <- snakemake@output[["large_sv_barplot"]]

dir.create(dirname(out_burden), showWarnings = FALSE, recursive = TRUE)

gene_panel <- read_gene_panel(gene_panel_csv)

read_one <- function(path, contrast) {
  if (!file.exists(path) || file.info(path)$size == 0) return(NULL)
  dt <- tryCatch(fread(path), error = function(e) NULL)
  if (is.null(dt) || nrow(dt) == 0) return(NULL)
  dt[, Contrast := contrast]
  dt
}

dt_list <- Map(read_one, sv_tsvs, contrasts)
dt_list <- Filter(Negate(is.null), dt_list)

if (length(dt_list) == 0) {
  message("No Manta SV records across any contrast; writing empty cohort SV outputs")
  fwrite_gz(data.table(), out_burden, sep = "\t")
  fwrite_gz(data.table(), out_gene_calls, sep = "\t")
  fwrite_gz(data.table(), out_recurrent, sep = "\t")
  fwrite_gz(data.table(), out_large_sv, sep = "\t")
  pdf(out_barplot); plot.new(); text(0.5, 0.5, "No SV data available"); dev.off()
  pdf(out_large_plot); plot.new(); text(0.5, 0.5, "No SV data available"); dev.off()
  quit(save = "no", status = 0)
}

all_dt <- rbindlist(dt_list, use.names = TRUE, fill = TRUE)
rm(dt_list); gc()

all_dt[, Is_Pass := Filter %in% c("PASS", ".")]
all_dt[is.na(SV_Type) | SV_Type == "", SV_Type := "UNKNOWN"]
if (!"Is_Primary_Mate" %in% names(all_dt)) all_dt[, Is_Primary_Mate := TRUE]  # older per-contrast TSVs without this column
all_dt[is.na(Is_Primary_Mate), Is_Primary_Mate := TRUE]

## Burden/gene counts below use only Is_Primary_Mate rows: Manta represents
## each BND (translocation) event as two mated VCF records (one per
## breakend, see summarize_manta_sv.R), so counting every record would
## double BND burden relative to DEL/DUP/INV/INS, which are single-record
## events.
count_dt <- all_dt[Is_Primary_Mate == TRUE]

## Burden: total Manta calls per contrast x SV type, PASS-only and all -------
burden <- count_dt[, .(
  n_total = .N,
  n_pass  = sum(Is_Pass)
), by = .(Contrast, SV_Type)]
setorder(burden, Contrast, SV_Type)
fwrite_gz(burden, out_burden, sep = "\t", quote = FALSE)

pdf(out_barplot, width = max(8, uniqueN(burden$Contrast) * 0.6), height = 6)
tryCatch({
  plot_dt <- if (pass_only) burden[n_pass > 0, .(Contrast, SV_Type, n = n_pass)] else burden[, .(Contrast, SV_Type, n = n_total)]
  if (nrow(plot_dt) == 0) {
    plot.new(); text(0.5, 0.5, "No SV calls remain to plot")
  } else {
    mat <- dcast(plot_dt, SV_Type ~ Contrast, value.var = "n", fun.aggregate = sum, fill = 0)
    sv_types <- mat$SV_Type
    mat[, SV_Type := NULL]
    mat_m <- t(as.matrix(mat))
    colnames(mat_m) <- sv_types
    par(mar = c(8, 4, 2, 8), xpd = TRUE)
    cols <- rainbow(ncol(mat_m))
    bp <- barplot(t(mat_m), beside = FALSE, col = cols, las = 2,
                  main = paste0("Manta SV burden per contrast (", if (pass_only) "PASS-only" else "all", ")"),
                  ylab = "SV count", cex.names = 0.7)
    legend("topright", inset = c(-0.18, 0), legend = sv_types, fill = cols, bty = "n", cex = 0.7)
  }
}, error = function(e) {
  plot.new(); text(0.5, 0.5, paste("Plot failed:", e$message))
})
dev.off()

## Large-scale event detection: DEL/DUP/INV segments >= large_sv_min_bp,
## i.e. Manta's breakpoint-based counterpart to a CNVkit-style
## amplification/deletion call. This is what lets you tell, straight from
## the cohort SV overview, WHICH contrasts carry large-scale rearrangements
## instead of only many small indel-like SVs - previously the burden table
## and barplot only counted events, with no size information at all. BND is
## excluded (no span/SV_Length - see summarize_manta_sv.R) and so is INS
## (Manta's inserted-sequence length is not indicative of a genomic deletion
## /amplification region). Cross-reference with the CNVkit cohort heatmap
## (results/cnv/cohort/cnv_heatmap.pdf) for independent, depth-based
## confirmation of the same regions/genes.
if ("SV_Length" %in% names(count_dt)) {
  large_dt <- count_dt[SV_Type %in% c("DEL", "DUP", "INV") & !is.na(SV_Length) & SV_Length >= large_sv_min_bp]
  if (pass_only) large_dt <- large_dt[Is_Pass == TRUE]
} else {
  message("SV_Length column not present (per-contrast TSV predates this feature; rerun summarize_manta_sv); skipping large-SV detection")
  large_dt <- count_dt[0]
  large_dt[, SV_Length := numeric(0)]
}

large_cols <- intersect(c("Contrast", "Chromosome", "Start", "End", "SV_Type", "SV_Length", "Filter", "Gene_Annotation"), names(large_dt))
large_out <- large_dt[, ..large_cols]
if (nrow(large_out) > 0) setorder(large_out, -SV_Length)
fwrite_gz(large_out, out_large_sv, sep = "\t", quote = FALSE)

pdf(out_large_plot, width = max(8, uniqueN(count_dt$Contrast) * 0.5), height = 6)
tryCatch({
  large_counts <- large_dt[, .N, by = Contrast]
  ## Include every contrast (even those with zero large SVs) so the absence
  ## of large-scale events is visible, not just omitted.
  all_contrasts_dt <- data.table(Contrast = sort(unique(count_dt$Contrast)))
  large_counts <- merge(all_contrasts_dt, large_counts, by = "Contrast", all.x = TRUE)
  large_counts[is.na(N), N := 0L]
  setorder(large_counts, -N, Contrast)
  if (all(large_counts$N == 0)) {
    plot.new()
    text(0.5, 0.5, paste0("No SV >= ", format(large_sv_min_bp, big.mark = ","), " bp (DEL/DUP/INV) in any contrast"))
  } else {
    par(mar = c(8, 4, 2, 2))
    barplot(large_counts$N, names.arg = large_counts$Contrast, las = 2,
            main = paste0("Large-scale SV (DEL/DUP/INV >= ", format(large_sv_min_bp, big.mark = ","), " bp) per contrast",
                           if (pass_only) ", PASS-only" else ""),
            ylab = "Count", cex.names = 0.7, col = "firebrick")
  }
}, error = function(e) {
  plot.new(); text(0.5, 0.5, paste("Plot failed:", e$message))
})
dev.off()

## Gene-level hit calls: explode Gene_Annotation into one row per gene ------
## symbol per SV record, restricted to Is_Pass by default (config
## sv_analysis.pass_only). Mirrors the CNVkit gain/loss "calls" table
## consumed by driver_candidates.R.
sv_for_genes <- if (pass_only) count_dt[Is_Pass == TRUE] else count_dt

extract_genes <- function(ann_string) {
  if (is.na(ann_string) || !nzchar(ann_string)) return(character(0))
  records <- unlist(strsplit(ann_string, ";", fixed = TRUE))
  genes <- vapply(records, function(rec) {
    fields <- strsplit(rec, "|", fixed = TRUE)[[1]]
    if (length(fields) >= gene_field_idx) trimws(fields[gene_field_idx]) else NA_character_
  }, character(1))
  genes <- genes[!is.na(genes) & nzchar(genes)]
  unique(genes)
}

if (nrow(sv_for_genes) > 0 && "Gene_Annotation" %in% names(sv_for_genes)) {
  gene_hits <- sv_for_genes[, .(Hugo_Symbol = extract_genes(Gene_Annotation)),
                             by = .(Contrast, SV_Type, Chromosome, Start, End, Filter)]
  gene_hits <- gene_hits[!is.na(Hugo_Symbol) & nzchar(Hugo_Symbol)]
} else {
  gene_hits <- data.table(Contrast = character(0), SV_Type = character(0), Chromosome = character(0),
                           Start = integer(0), End = integer(0), Filter = character(0), Hugo_Symbol = character(0))
}

gene_calls <- gene_hits[, .(
  n_sv = .N,
  SV_Types = paste(sort(unique(SV_Type)), collapse = ";")
), by = .(Hugo_Symbol, Contrast)]
gene_calls[, In_Gene_Panel := Hugo_Symbol %in% gene_panel]
fwrite_gz(gene_calls, out_gene_calls, sep = "\t", quote = FALSE)

## Cohort-wide gene recurrence ranking, gene-panel-aware like driver_candidates.R
if (nrow(gene_calls) > 0) {
  recur <- gene_calls[, .(
    n_contrasts = uniqueN(Contrast),
    contrasts   = paste(sort(unique(Contrast)), collapse = ";"),
    n_sv_total  = sum(n_sv),
    SV_Types    = paste(sort(unique(unlist(strsplit(SV_Types, ";")))), collapse = ";")
  ), by = Hugo_Symbol]
  recur[, In_Gene_Panel := Hugo_Symbol %in% gene_panel]
  recur[, Score := n_contrasts + 0.5 * In_Gene_Panel]
  setorder(recur, -Score, -n_contrasts, Hugo_Symbol)
} else {
  recur <- data.table()
}
fwrite_gz(recur, out_recurrent, sep = "\t", quote = FALSE)

message("Cohort SV summary: ", nrow(count_dt), " SV events (BND mate pairs deduplicated) across ", uniqueN(all_dt$Contrast),
        " contrasts -> ", nrow(gene_calls), " gene x contrast hits, ", nrow(recur), " recurrent genes")
