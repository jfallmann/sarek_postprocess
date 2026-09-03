rule actionability_report:
    input:
        cohort_ranked=f"{OUTDIR}/drivers/cohort_ranked_drivers.tsv",
    output:
        report=f"{OUTDIR}/actionability/actionability_report.tsv",
    params:
        cancer_gene_census_tsv=config["actionability"].get("cancer_gene_census_tsv", ""),
        oncokb_gene_list_tsv=config["actionability"].get("oncokb_gene_list_tsv", ""),
        civic_variant_summary_tsv=config["actionability"].get("civic_variant_summary_tsv", ""),
    log:
        f"{OUTDIR}/logs/actionability_report.log",
    conda:
        "../envs/r_env.yaml"
    script:
        "../scripts/actionability_report.R"
