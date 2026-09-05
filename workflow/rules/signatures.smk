def signature_caller_vcfs(wildcards):
    df = get_contrasts_df()
    caller = config["signatures"]["caller"]
    sub = df[df.caller == caller]
    return sub["vcf_path"].tolist()


def signature_contrast_names(wildcards):
    df = get_contrasts_df()
    caller = config["signatures"]["caller"]
    sub = df[df.caller == caller]
    return sub["contrast"].tolist()


rule mutational_signatures:
    input:
        mutect2_vcfs=signature_caller_vcfs,
    output:
        done=f"{OUTDIR}/signatures/cohort/.done",
    params:
        contrasts=signature_contrast_names,
        genome_package=config["signatures"]["genome_package"],
        reference_set=config["signatures"]["reference_set"],
        out_dir=f"{OUTDIR}/signatures/cohort",
        common_r=f"{workflow.basedir}/scripts/common.R",
    log:
        f"{OUTDIR}/logs/mutational_signatures.log",
    conda:
        "../envs/r_env.yaml"
    script:
        "../scripts/mutational_signatures.R"
