#!/usr/bin/env Rscript
## Standalone runner for driver_candidates.R, for use OUTSIDE the Snakemake
## pipeline - e.g. ranking driver candidates across arbitrary sample/group
## consensus MAFs that were never Sarek "contrasts" at all (a normal-control
## population vs. a treatment-sensitive population vs. a treatment-resistant
## population, or any other future dataset).
##
## driver_candidates.R itself is NOT modified - it is sourced unchanged, only
## fed a minimal mock of the `snakemake` S4 object Snakemake's `script:`
## directive normally injects (input$consensus_mafs, params$..., etc.), so
## this stays in sync with the pipeline's own driver logic automatically.
##
## INPUT REQUIREMENT: one "consensus" MAF/TSV per group, in the same
## maftools-MAF-column format written by consolidate_mafs.R's *.consensus.tsv.gz
## (Hugo_Symbol, Chromosome, Start_Position, End_Position, Reference_Allele,
## Tumor_Seq_Allele2, Tumor_Sample_Barcode, Variant_Classification, ...). If
## you already ran this pipeline, just point at
## results/maf_consolidated/<vcset>/<contrast>.consensus.tsv.gz. If you did
## not, you can hand it any MAF-like TSV with those columns - e.g. an
## uncalled union/consensus you built yourself from other callers - for
## reuse on non-Sarek data.
##
## Usage:
##   Rscript run_driver_candidates_standalone.R \
##     --group Normal=normal.consensus.tsv.gz \
##     --group Sensitive=sensitive.consensus.tsv.gz \
##     --group Resistant=resistant.consensus.tsv.gz \
##     --sensitive-group Normal \
##     --out cohort_ranked_drivers.tsv.gz \
##     --out-per-group-dir per_group/ \
##     [--gene-panel gene_panel.csv] \
##     [--baseline-mutations baseline.tsv] \
##     [--cnv-calls cnv_calls.tsv.gz] [--sv-calls sv_gene_calls.tsv.gz] \
##     [--consensus-min-callers 2] [--min-depth 20] [--min-alt-count 3] \
##     [--min-vaf 0.05] [--max-gnomad-af 0.001] \
##     [--exclude-gene-regex '(?i)^MUC|^TTN|^HLA|^NEB|^OBSCN|^USH2A']
##
## --group can be given multiple times, one per population/sample/contrast.
## Its consensus MAF is expected to ALREADY be caller-agreement-filtered
## (i.e. this script does not re-run consolidate_mafs.R's union/consensus
## step - only the driver ranking on top of already-consensus MAFs).
##
## More than one TSV per group (e.g. several samples/replicates that should
## be pooled into one "Resistant" population rather than ranked separately)
## is supported two ways - pick whichever reads more naturally:
##   --group Resistant=sample1.consensus.tsv.gz,sample2.consensus.tsv.gz
##   --group Resistant=sample1.consensus.tsv.gz --group Resistant=sample2.consensus.tsv.gz
## Both are read independently and then combined under the same group label
## before ranking - driver_candidates.R groups by the label you give here,
## not by path, so this needs no change to that script.
##
## --sensitive-group can be given multiple times to name the "control" pool
## for the optional exploratory mafCompare step (same role as
## config.yaml's sensitive_contrasts).

args <- commandArgs(trailingOnly = TRUE)

get_opt <- function(flag, default = NULL, multiple = FALSE) {
  idx <- which(args == flag)
  if (length(idx) == 0) return(default)
  vals <- args[idx + 1]
  if (multiple) vals else vals[1]
}

group_args <- get_opt("--group", character(0), multiple = TRUE)
if (length(group_args) == 0) {
  stop("At least one --group Label=path.tsv.gz is required")
}
parsed <- strsplit(group_args, "=", fixed = TRUE)
labels_raw <- vapply(parsed, `[`, character(1), 1)
paths_raw  <- vapply(parsed, function(x) paste(x[-1], collapse = "="), character(1))
## Split each --group's path portion on commas, so one flag can list several
## TSVs for the same label; the label is repeated for each resulting path.
paths_split  <- strsplit(paths_raw, ",", fixed = TRUE)
group_labels <- unlist(Map(function(lbl, ps) rep(lbl, length(ps)), labels_raw, paths_split))
group_paths  <- unlist(paths_split)

out_cohort_ranked   <- get_opt("--out", "cohort_ranked_drivers.tsv.gz")
out_per_contrast_dir <- get_opt("--out-per-group-dir", "per_group")
gene_panel_csv      <- get_opt("--gene-panel", "")
baseline_tsv        <- get_opt("--baseline-mutations", "")
cnv_calls_tsv       <- get_opt("--cnv-calls", "")
sv_calls_tsv        <- get_opt("--sv-calls", "")
sensitive_contrasts <- get_opt("--sensitive-group", character(0), multiple = TRUE)

## Defaults mirror config.yaml's `filtering:` block as of this writing - keep
## in sync if you change the pipeline defaults there.
filtering <- list(
  consensus_min_callers = as.integer(get_opt("--consensus-min-callers", "2")),
  min_depth             = as.numeric(get_opt("--min-depth", "20")),
  min_alt_count         = as.numeric(get_opt("--min-alt-count", "3")),
  min_vaf               = as.numeric(get_opt("--min-vaf", "0.05")),
  max_gnomad_af         = as.numeric(get_opt("--max-gnomad-af", "0.001")),
  exclude_gene_regex    = get_opt("--exclude-gene-regex", "(?i)^MUC|^TTN|^HLA|^NEB|^OBSCN|^USH2A")
)

script_dir <- dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)))
if (length(script_dir) == 0 || !nzchar(script_dir)) script_dir <- "."

setClass("MockSnakemake", representation(input = "list", output = "list", params = "list"))

snakemake <- new(
  "MockSnakemake",
  input  = list(consensus_mafs = group_paths),
  output = list(cohort_ranked = out_cohort_ranked),
  params = list(
    common_r             = file.path(script_dir, "common.R"),
    contrasts             = group_labels,
    sensitive_contrasts   = sensitive_contrasts,
    filtering             = filtering,
    gene_panel_csv        = gene_panel_csv,
    baseline_mutations_tsv = baseline_tsv,
    cnv_calls_tsv         = cnv_calls_tsv,
    sv_calls_tsv          = sv_calls_tsv,
    out_per_contrast_dir  = out_per_contrast_dir
  )
)

source(file.path(script_dir, "driver_candidates.R"))
