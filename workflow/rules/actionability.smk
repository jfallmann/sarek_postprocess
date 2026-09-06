rule actionability_report:
    input:
        # Actionability is based on the more inclusive "custom" driver track;
        # use the "protein_coding" track's output if you want the strict one.
        cohort_ranked=f"{OUTDIR}/drivers/custom/cohort_ranked_drivers.tsv.gz",
    output:
        report=f"{OUTDIR}/actionability/actionability_report.tsv.gz",
    params:
        common_r=f"{workflow.basedir}/scripts/common.R",
        cancer_gene_census_tsv=config["actionability"].get("cancer_gene_census_tsv", ""),
        oncokb_gene_list_tsv=config["actionability"].get("oncokb_gene_list_tsv", ""),
        civic_variant_summary_tsv=config["actionability"].get("civic_variant_summary_tsv", ""),
    log:
        f"{OUTDIR}/logs/actionability_report.log",
    resources:
        tmpdir=R_TMPDIR,
        mem_mb=2000,
    conda:
        "../envs/r_env.yaml"
    script:
        "../scripts/actionability_report.R"
