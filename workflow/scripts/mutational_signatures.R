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

## Known proposed aetiologies for the COSMIC SBS reference signatures
## (https://cancer.sanger.ac.uk/signatures/sbs/), so a signature ID like
## "SBS4" in the outputs can be read as "tobacco smoking" without having to
## look it up separately. Restricted to the widely-established/common
## signatures; anything not listed here (rarer/less-characterised SBS IDs)
## is reported with "Unknown / not in built-in lookup - see COSMIC" rather
## than guessed at.
SBS_AETIOLOGY <- c(
  SBS1  = "Clock-like (spontaneous deamination of 5-methylcytosine)",
  SBS2  = "APOBEC cytidine deaminase activity",
  SBS3  = "Homologous recombination deficiency (BRCA1/BRCA2)",
  SBS4  = "Tobacco smoking",
  SBS5  = "Clock-like (unknown mechanism)",
  SBS6  = "Defective DNA mismatch repair",
  SBS7a = "Ultraviolet light exposure",
  SBS7b = "Ultraviolet light exposure",
  SBS7c = "Ultraviolet light exposure",
  SBS7d = "Ultraviolet light exposure",
  SBS8  = "Unknown (possibly late-replicating DNA damage)",
  SBS9  = "Polymerase eta somatic hypermutation",
  SBS10a = "Defective DNA polymerase epsilon proofreading (POLE)",
  SBS10b = "Defective DNA polymerase epsilon proofreading (POLE)",
  SBS11 = "Temozolomide treatment",
  SBS13 = "APOBEC cytidine deaminase activity",
  SBS14 = "Defective DNA mismatch repair + POLE proofreading",
  SBS15 = "Defective DNA mismatch repair",
  SBS17a = "Unknown (associated with gastric/oesophageal cancer)",
  SBS17b = "Unknown (associated with gastric/oesophageal cancer)",
  SBS18 = "Reactive oxygen species damage",
  SBS20 = "Defective DNA mismatch repair + POLE proofreading",
  SBS21 = "Defective DNA mismatch repair",
  SBS22 = "Aristolochic acid exposure",
  SBS24 = "Aflatoxin exposure",
  SBS25 = "Chemotherapy treatment",
  SBS26 = "Defective DNA mismatch repair",
  SBS29 = "Tobacco chewing",
  SBS30 = "Defective base excision repair (NTHL1)",
  SBS31 = "Platinum chemotherapy treatment",
  SBS32 = "Azathioprine treatment",
  SBS35 = "Platinum chemotherapy treatment",
  SBS36 = "Defective base excision repair (MUTYH)",
  SBS40 = "Unknown (correlated with age in some cancers)",
  SBS44 = "Defective DNA mismatch repair"
)

sig_dt <- as.data.table(contribution, keep.rownames = "signature")
fwrite_gz(sig_dt, file.path(out_dir, "signature_contributions.tsv.gz"), sep = "\t")

## Companion mapping file, restricted to the signatures actually fit in this
## cohort (not the full ~80-signature COSMIC catalogue), so it can be
## cross-referenced directly against signature_contributions.tsv.gz /
## the two plots below without hunting through the full reference set.
sig_meaning_dt <- data.table(
  Signature = sig_dt$signature,
  Proposed_Aetiology = ifelse(sig_dt$signature %in% names(SBS_AETIOLOGY),
                               SBS_AETIOLOGY[sig_dt$signature],
                               "Unknown / not in built-in lookup - see COSMIC (https://cancer.sanger.ac.uk/signatures/sbs/)")
)
fwrite_gz(sig_meaning_dt, file.path(out_dir, "signature_meanings.tsv.gz"), sep = "\t")

## Both plots put one bar/column per contrast; with the default fixed-size
## PDFs long contrast names (e.g. "*_vs_Bulk_sensitive") get clipped off the
## x-axis. Scale the canvas with the number/length of contrasts instead, and
## angle the heatmap's x-axis labels so they fit without truncation. The
## previous version of this scaling (0.12 in/char, no extra plot margin)
## still let the last, longest rotated label run past the right edge of the
## canvas on real cohorts (e.g. 15 contrasts, 37-char names) - bump the
## per-char factor and add an explicit right-hand plot margin sized to the
## longest label so ggplot has room to actually draw it instead of clipping
## it at the device boundary.
n_contrasts <- ncol(contribution)
max_contrast_len <- max(nchar(colnames(contribution)), 1)
heatmap_width  <- max(10, n_contrasts * 0.4 + max_contrast_len * 0.25)
barplot_height <- max(6, n_contrasts * 0.3 + max_contrast_len * 0.05)

pdf(file.path(out_dir, "signature_contribution_heatmap.pdf"), width = heatmap_width, height = 8)
heatmap_plot <- plot_contribution_heatmap(contribution, cluster_samples = FALSE)
print(heatmap_plot + ggplot2::theme(
  axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
  plot.margin = ggplot2::margin(t = 5, r = max_contrast_len * 3, b = 5, l = 5)
))
dev.off()

## coord_flip = TRUE puts contrasts on the (now horizontal) axis, so it is the
## plot height, not width, that must grow with the number/length of contrasts.
pdf(file.path(out_dir, "signature_contribution_barplot.pdf"), width = 10, height = barplot_height)
print(plot_contribution(contribution, coord_flip = TRUE, mode = "relative"))
dev.off()

file.create(done_marker)
message("Mutational signature analysis written to ", out_dir)
