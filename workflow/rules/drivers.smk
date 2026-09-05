def cohort_consensus_tsvs(wildcards):
    return [f"{OUTDIR}/maf_consolidated/{c}.consensus.tsv.gz" for c in all_contrasts(wildcards)]


def driver_cnv_input(wildcards):
    if not config.get("cnv", {}).get("enabled", False):
        return []
    # Only wire this dependency in if CNVkit discovery actually found
    # anything - otherwise cohort_cnv_matrix would fail with no inputs.
    if len(get_cnv_df()) == 0:
        return []
    return [f"{OUTDIR}/cnv/cohort/cnv_calls.tsv.gz"]


rule driver_candidates:
    input:
        consensus_mafs=cohort_consensus_tsvs,
        cnv_calls=driver_cnv_input,
    output:
        cohort_ranked=f"{OUTDIR}/drivers/cohort_ranked_drivers.tsv.gz",
    params:
        common_r=f"{workflow.basedir}/scripts/common.R",
        contrasts=all_contrasts,
        sensitive_contrasts=config.get("sensitive_contrasts", []),
        filtering=config["filtering"],
        gene_panel_csv=config["gene_panel_csv"],
        baseline_mutations_tsv=config.get("baseline_mutations_tsv", ""),
        cnv_calls_tsv=lambda wc, input: (input.cnv_calls[0] if len(input.cnv_calls) > 0 else ""),
        out_per_contrast_dir=f"{OUTDIR}/drivers/per_contrast",
    log:
        f"{OUTDIR}/logs/driver_candidates.log",
    conda:
        "../envs/r_env.yaml"
    script:
        "../scripts/driver_candidates.R"
