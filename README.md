# sarek_postprocess

Standalone Snakemake pipeline that takes the output of an nf-core/Sarek run
(Strelka2, Mutect2, Manta, VEP+SnpEff `--everything` annotation) and derives:

1. **Candidate driver/resistance genes** per contrast and cohort-wide,
   ranked by cross-caller consensus, recurrence across contrasts, and
   membership in a curated driver/resistance gene panel.
2. **Cohort-wide summary plots**: oncoplot, MAF summary dashboard, and
   mutational signature fitting (COSMIC reference signatures via
   MutationalPatterns) across all contrasts.
3. *(optional, lowest priority)* An **actionability report** cross-referencing
   candidate genes against static, offline Cancer Gene Census / OncoKB /
   CIViC exports.

Manta structural variants are summarized separately (not MAF-compatible).

## How it works

Sarek's own `--outdir` already encodes tumor-only vs tumor-normal pairing in
its folder names (`annotation/<caller>/<contrast>/...`, where `<contrast>`
is `<tumor>` or `<tumor>_vs_<normal>`). `workflow/scripts/discover_contrasts.py`
walks that structure directly, resolves tumor/normal sample IDs from the VCF
header via `bcftools`, and only reads the Sarek samplesheet for metadata
(patient/sex/status) - not to re-derive pairing. Everything downstream is
driven by that discovery step (a Snakemake checkpoint), so the same pipeline
adapts to however many samples/contrasts an outdir actually contains.

## Usage

1. Copy `config/config.yaml` and fill in `samplesheet`, `sarek_outdir`,
   `outdir`, `ref_fasta`, and the `vcf2maf` VEP cache paths (same offline
   cache Sarek used, since Sarek was run with `--everything`).
2. Adjust `config/driver_resistance_genes.csv` if you want a different gene
   panel, and `sensitive_contrasts` in the config if you want the optional
   mafCompare/forestPlot step (resistant clone vs pooled sensitive
   contrasts).
3. Dry-run first:
   ```bash
   snakemake -s workflow/Snakefile --configfile config/my_run.yaml -n
   ```
4. Run on the cluster with a Snakemake executor plugin/profile (e.g. Slurm):
   ```bash
   snakemake -s workflow/Snakefile --configfile config/my_run.yaml \
     --use-conda --executor slurm --jobs 50
   ```
   `--use-conda` builds the two environments in `workflow/envs/` (vcf2maf+VEP,
   R+maftools/MutationalPatterns) per rule.

## Layout

```
config/config.yaml                  template config (one per Sarek run/cohort)
config/driver_resistance_genes.csv  curated gene panel used for ranking
workflow/Snakefile                  entry point, checkpoint-driven DAG
workflow/rules/*.smk                one rule file per pipeline stage
workflow/scripts/                   discovery, vcf2maf, R analysis scripts
workflow/envs/                      conda envs for vcf2maf/VEP and R
test/fixtures/                      directory-layout-only smoke test fixtures
```

## Notes

- Not covered here but worth evaluating separately: `nf-core/tumourevo` for
  subclonal deconvolution/clone-tree inference, and PCGR/CPSR as a more
  complete offline actionability reporter than the static-table join used
  in `actionability_report.R`.
- VEP was run with `--everything` (no extra plugins), which already provides
  gnomAD AF and ClinVar/PHENO fields used for filtering - no VEP rerun is
  required for goals 1-2.
