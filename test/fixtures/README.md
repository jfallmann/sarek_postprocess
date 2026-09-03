# Test fixtures

`fake_sarek_output/` mimics the directory/filename layout Sarek 3.x produces
(empty placeholder `.vcf.gz` files, no real variant data) and
`fake_samplesheet.csv` mimics a Sarek samplesheet. They exist to smoke-test
`workflow/scripts/discover_contrasts.py` and a Snakemake dry-run of the
whole DAG without needing bcftools/VEP/vcf2maf/R installed - they cannot be
used to test the actual bioinformatics (VCF parsing, MAF conversion,
signature fitting), which requires the real environment. Real end-to-end
validation should be done on the cluster.

Smoke test:

```bash
python3 workflow/scripts/discover_contrasts.py \
  --sarek-outdir test/fixtures/fake_sarek_output \
  --samplesheet test/fixtures/fake_samplesheet.csv \
  --callers mutect2,strelka --include-manta --out /tmp/contrasts.tsv
```

(`bcftools` is required to resolve tumor/normal sample IDs from the VCF
header; without it, `tumor_id`/`normal_id` stay blank but directory/caller
discovery still runs.)
