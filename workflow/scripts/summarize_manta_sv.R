## Summarize Manta structural variants per contrast into a flat TSV
## (BND/DEL/DUP/INV per gene, from the SnpEff/VEP ANN/CSQ field if present).
## Kept separate from the SNV/indel MAF track since Manta VCFs are not
## MAF-compatible.

suppressPackageStartupMessages({
  library(data.table)
  library(VariantAnnotation)
})

source(snakemake@params[["common_r"]])

vcf_path <- snakemake@input[["vcf"]]
contrast <- snakemake@wildcards[["contrast"]]
out_tsv  <- snakemake@output[["sv_tsv"]]

dir.create(dirname(out_tsv), showWarnings = FALSE, recursive = TRUE)

if (!file.exists(vcf_path) || file.info(vcf_path)$size == 0) {
  message("Manta VCF missing/empty for ", contrast)
  fwrite_gz(data.table(), out_tsv, sep = "\t")
  quit(save = "no", status = 0)
}

vcf <- tryCatch(readVcf(vcf_path, genome = "GRCh38"), error = function(e) {
  message("readVcf failed for ", vcf_path, ": ", e$message)
  NULL
})

if (is.null(vcf) || nrow(vcf) == 0) {
  message("No SVs in ", vcf_path)
  fwrite_gz(data.table(), out_tsv, sep = "\t")
  quit(save = "no", status = 0)
}

info_df <- as.data.table(info(vcf))
gr <- rowRanges(vcf)

## Concatenate BOTH annotators' output when a VCF was run through SnpEff
## (ANN) and VEP (CSQ), rather than the previous else-if which silently
## dropped CSQ whenever ANN was also present. Both formats put the gene
## symbol at the same pipe-delimited field position (see
## sv_analysis.gene_field_index in config.yaml), so downstream gene
## extraction doesn't need to know which annotator a given sub-record came
## from.
list_to_joined <- function(col) {
  vapply(col, function(x) paste(unique(x), collapse = ";"), character(1))
}
ann_txt <- if ("ANN" %in% names(info_df)) list_to_joined(info_df$ANN) else rep("", nrow(info_df))
csq_txt <- if ("CSQ" %in% names(info_df)) list_to_joined(info_df$CSQ) else rep("", nrow(info_df))
gene_annotation <- ifelse(nzchar(ann_txt) & nzchar(csq_txt), paste(ann_txt, csq_txt, sep = ";"),
                    ifelse(nzchar(ann_txt), ann_txt, ifelse(nzchar(csq_txt), csq_txt, NA_character_)))

## Manta represents each BND (translocation/breakend) event as a MATE PAIR:
## two separate VCF records, one per breakend, cross-referenced by
## MATEID/ID. Counting both mates as separate SV events doubles BND burden
## relative to DEL/DUP/INV/INS (which are single-record events) - flag one
## mate per pair as primary here (kept ID < MateID, an arbitrary but
## deterministic and stable tie-break) so cohort_sv_summary.R can count/rank
## by unique event instead of by VCF record. Both records are still kept in
## this per-contrast file for manual inspection of both breakpoints.
ids <- names(rowRanges(vcf))
if (is.null(ids)) ids <- as.character(VariantAnnotation::ID(vcf))
mate_ids <- if ("MATEID" %in% names(info_df)) {
  vapply(info_df$MATEID, function(x) if (length(x) > 0) as.character(x[[1]]) else NA_character_, character(1))
} else {
  rep(NA_character_, nrow(info_df))
}
is_primary <- ifelse(is.na(mate_ids), TRUE, ids < mate_ids)

start_pos <- GenomicRanges::start(gr)
end_pos   <- if ("END" %in% names(info_df)) info_df$END else start_pos
sv_type   <- if ("SVTYPE" %in% names(info_df)) info_df$SVTYPE else NA_character_

## Event size: prefer Manta's own SVLEN (signed - negative for deletions in
## some Manta versions, hence abs()); fall back to END-START for DEL/DUP/INV
## when SVLEN is absent. Left NA for BND (a breakpoint between two loci, not
## a span - "size" is meaningless/undefined for a translocation) and for
## INS/other types with neither SVLEN nor a usable END. This is what lets
## cohort_sv_summary.R distinguish "many small indel-like SVs" from "one
## chromosome-arm-scale deletion/duplication" per contrast.
svlen_raw <- if ("SVLEN" %in% names(info_df)) {
  vapply(info_df$SVLEN, function(x) if (length(x) > 0) as.numeric(x[[1]]) else NA_real_, numeric(1))
} else {
  rep(NA_real_, nrow(info_df))
}
sv_length <- abs(svlen_raw)
needs_fallback <- is.na(sv_length) & sv_type %in% c("DEL", "DUP", "INV")
sv_length[needs_fallback] <- abs(end_pos[needs_fallback] - start_pos[needs_fallback])

sv_dt <- data.table(
  Contrast   = contrast,
  Chromosome = as.character(GenomicRanges::seqnames(gr)),
  Start      = start_pos,
  End        = end_pos,
  SV_Type    = sv_type,
  SV_Length  = sv_length,
  Filter     = as.character(VariantAnnotation::filt(vcf)),
  ID         = ids,
  Mate_ID    = mate_ids,
  Is_Primary_Mate = is_primary,
  Gene_Annotation = gene_annotation
)

fwrite_gz(sv_dt, out_tsv, sep = "\t", quote = FALSE)
message("Manta SV summary for ", contrast, ": ", nrow(sv_dt), " records -> ", out_tsv)
