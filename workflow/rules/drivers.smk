def cohort_consensus_tsvs(wildcards):
    return [f"{OUTDIR}/maf_consolidated/{c}.consensus.tsv" for c in all_contrasts(wildcards)]


rule driver_candidates:
    input:
        consensus_mafs=cohort_consensus_tsvs,
    output:
        cohort_ranked=f"{OUTDIR}/drivers/cohort_ranked_drivers.tsv",
    params:
        common_r=f"{workflow.basedir}/scripts/common.R",
        contrasts=all_contrasts,
        sensitive_contrasts=config.get("sensitive_contrasts", []),
        filtering=config["filtering"],
        gene_panel_csv=config["gene_panel_csv"],
        baseline_mutations_tsv=config.get("baseline_mutations_tsv", ""),
        out_per_contrast_dir=f"{OUTDIR}/drivers/per_contrast",
    log:
        f"{OUTDIR}/logs/driver_candidates.log",
    conda:
        "../envs/r_env.yaml"
    script:
        "../scripts/driver_candidates.R"
