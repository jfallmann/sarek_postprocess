## Mutational signature fitting per contrast + cohort-wide heatmap.
##
## Uses MutationalPatterns (Bioconductor) fitting PASS Mutect2 SNVs against
## the COSMIC reference signature set - NOT de novo NMF extraction, since a
## handful of clones/contrasts do not give enough mutations/samples for
## stable de novo signature discovery.

suppressPackageStartupMessages({
  library(data.table)
  library(MutationalPatterns)
  library(BSgenome)
})

source(snakemake@params[["common_r"]])

vcf_paths   <- unlist(snakemake@input[["mutect2_vcfs"]])
contrasts   <- unlist(snakemake@params[["contrasts"]])
genome_pkg  <- snakemake@params[["genome_package"]]
ref_set     <- snakemake@params[["reference_set"]]
out_dir     <- snakemake@params[["out_dir"]]
done_marker <- snakemake@output[["done"]]

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

if (!requireNamespace(genome_pkg, quietly = TRUE)) {
  stop("Genome package ", genome_pkg, " is not installed. Install it via BiocManager on the ",
       "execution environment used for this rule before running on the cluster.")
}
library(genome_pkg, character.only = TRUE)

vcf_paths <- vcf_paths[file.exists(vcf_paths)]
contrasts <- contrasts[seq_along(vcf_paths)]

if (length(vcf_paths) == 0) {
  message("No Mutect2 VCFs available for signature analysis; skipping")
  file.create(file.path(out_dir, "NO_VCFS.txt"))
  file.create(done_marker)
  quit(save = "no", status = 0)
}

grl <- tryCatch(
  read_vcfs_as_granges(vcf_paths, sample_names = contrasts, genome = genome_pkg, type = "snv"),
  error = function(e) {
    message("read_vcfs_as_granges failed: ", e$message)
    NULL
  }
)

if (is.null(grl) || length(grl) == 0) {
  message("No usable SNVs across contrasts; skipping signature fitting")
  file.create(file.path(out_dir, "NO_VARIANTS.txt"))
  file.create(done_marker)
  quit(save = "no", status = 0)
}

mut_mat <- mut_matrix(vcf_list = grl, ref_genome = genome_pkg)

signatures <- tryCatch(
  get_known_signatures(muttype = "snv", source = "COSMIC", genome = "GRCh38"),
  error = function(e) {
    message("get_known_signatures failed (check MutationalPatterns/BSgenome versions): ", e$message)
    NULL
  }
)

if (is.null(signatures)) {
  message("Reference signature set unavailable; writing raw mutation matrix only")
  fwrite_gz(as.data.table(mut_mat, keep.rownames = "context"), file.path(out_dir, "mutation_matrix.tsv.gz"), sep = "\t")
  file.create(done_marker)
  quit(save = "no", status = 0)
}

fit_res <- fit_to_signatures(mut_mat, signatures)
contribution <- fit_res$contribution

fwrite_gz(as.data.table(contribution, keep.rownames = "signature"),
       file.path(out_dir, "signature_contributions.tsv.gz"), sep = "\t")

## Both plots put one bar/column per contrast; with the default fixed-size
## PDFs long contrast names (e.g. "*_vs_Bulk_sensitive") get clipped off the
## x-axis. Scale the canvas with the number/length of contrasts instead, and
## angle the heatmap's x-axis labels so they fit without truncation.
n_contrasts <- ncol(contribution)
max_contrast_len <- max(nchar(colnames(contribution)), 1)
heatmap_width  <- max(10, n_contrasts * 0.4 + max_contrast_len * 0.12)
barplot_height <- max(6, n_contrasts * 0.3 + max_contrast_len * 0.05)

pdf(file.path(out_dir, "signature_contribution_heatmap.pdf"), width = heatmap_width, height = 8)
heatmap_plot <- plot_contribution_heatmap(contribution, cluster_samples = FALSE)
print(heatmap_plot + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)))
dev.off()

## coord_flip = TRUE puts contrasts on the (now horizontal) axis, so it is the
## plot height, not width, that must grow with the number/length of contrasts.
pdf(file.path(out_dir, "signature_contribution_barplot.pdf"), width = 10, height = barplot_height)
print(plot_contribution(contribution, coord_flip = TRUE, mode = "relative"))
dev.off()

file.create(done_marker)
message("Mutational signature analysis written to ", out_dir)
