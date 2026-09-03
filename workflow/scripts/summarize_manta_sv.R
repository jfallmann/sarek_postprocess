## Summarize Manta structural variants per contrast into a flat TSV
## (BND/DEL/DUP/INV per gene, from the SnpEff/VEP ANN/CSQ field if present).
## Kept separate from the SNV/indel MAF track since Manta VCFs are not
## MAF-compatible.

suppressPackageStartupMessages({
  library(data.table)
  library(VariantAnnotation)
})

vcf_path <- snakemake@input[["vcf"]]
contrast <- snakemake@wildcards[["contrast"]]
out_tsv  <- snakemake@output[["sv_tsv"]]

dir.create(dirname(out_tsv), showWarnings = FALSE, recursive = TRUE)

if (!file.exists(vcf_path) || file.info(vcf_path)$size == 0) {
  message("Manta VCF missing/empty for ", contrast)
  fwrite(data.table(), out_tsv, sep = "\t")
  quit(save = "no", status = 0)
}

vcf <- tryCatch(readVcf(vcf_path, genome = "GRCh38"), error = function(e) {
  message("readVcf failed for ", vcf_path, ": ", e$message)
  NULL
})

if (is.null(vcf) || nrow(vcf) == 0) {
  message("No SVs in ", vcf_path)
  fwrite(data.table(), out_tsv, sep = "\t")
  quit(save = "no", status = 0)
}

info_df <- as.data.table(info(vcf))
gr <- rowRanges(vcf)

sv_dt <- data.table(
  Contrast   = contrast,
  Chromosome = as.character(GenomicRanges::seqnames(gr)),
  Start      = GenomicRanges::start(gr),
  End        = if ("END" %in% names(info_df)) info_df$END else GenomicRanges::start(gr),
  SV_Type    = if ("SVTYPE" %in% names(info_df)) info_df$SVTYPE else NA_character_,
  Filter     = as.character(VariantAnnotation::filt(vcf)),
  Gene_Annotation = if ("ANN" %in% names(info_df)) {
    vapply(info_df$ANN, function(x) paste(unique(x), collapse = ";"), character(1))
  } else if ("CSQ" %in% names(info_df)) {
    vapply(info_df$CSQ, function(x) paste(unique(x), collapse = ";"), character(1))
  } else {
    NA_character_
  }
)

fwrite(sv_dt, out_tsv, sep = "\t", quote = FALSE)
message("Manta SV summary for ", contrast, ": ", nrow(sv_dt), " records -> ", out_tsv)
