def cohort_consensus_tsvs(wildcards):
    return [f"{OUTDIR}/maf_consolidated/{wildcards.vcset}/{c}.consensus.tsv.gz" for c in all_contrasts(wildcards)]


def driver_cnv_input(wildcards):
    if not config.get("cnv", {}).get("enabled", False):
        return []
    # Only wire this dependency in if CNVkit discovery actually found
    # anything - otherwise cohort_cnv_matrix would fail with no inputs.
    if len(get_cnv_df()) == 0:
        return []
    return [f"{OUTDIR}/cnv/cohort/cnv_calls.tsv.gz"]


def driver_sv_input(wildcards):
    if not config.get("sv_analysis", {}).get("enabled", False):
        return []
    # Only wire this dependency in if Manta discovery actually found
    # anything - otherwise cohort_sv_summary would have no inputs.
    if len(manta_contrasts(wildcards)) == 0:
        return []
    return [f"{OUTDIR}/sv/cohort/sv_gene_calls.tsv.gz"]


rule driver_candidates:
    input:
        consensus_mafs=cohort_consensus_tsvs,
        cnv_calls=driver_cnv_input,
        sv_calls=driver_sv_input,
    output:
        cohort_ranked=f"{OUTDIR}/drivers/{{vcset}}/cohort_ranked_drivers.tsv.gz",
    params:
        common_r=f"{workflow.basedir}/scripts/common.R",
        contrasts=all_contrasts,
        sensitive_contrasts=config.get("sensitive_contrasts", []),
        filtering=config["filtering"],
        gene_panel_csv=config["gene_panel_csv"],
        baseline_mutations_tsv=config.get("baseline_mutations_tsv", ""),
        cnv_calls_tsv=lambda wc, input: (input.cnv_calls[0] if len(input.cnv_calls) > 0 else ""),
        sv_calls_tsv=lambda wc, input: (input.sv_calls[0] if len(input.sv_calls) > 0 else ""),
        out_per_contrast_dir=f"{OUTDIR}/drivers/{{vcset}}/per_contrast",
    log:
        f"{OUTDIR}/logs/driver_candidates/{{vcset}}.log",
    resources:
        tmpdir=R_TMPDIR,
        mem_mb=8000,
    conda:
        "../envs/r_env.yaml"
    script:
        "../scripts/driver_candidates.R"
