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


def signature_caller_vcf_for(wildcards):
    df = get_contrasts_df()
    caller = config["signatures"]["caller"]
    sub = df[(df.caller == caller) & (df.contrast == wildcards.contrast)]
    if sub.empty:
        raise WorkflowError(f"No {caller} VCF found for signature contrast={wildcards.contrast}")
    return sub["vcf_path"].iloc[0]


# WGS Mutect2 VCFs carry large numbers of non-PASS candidate/germline-leakage
# records; fitting signatures on raw calls dilutes/biases the contributions
# (every other stage of the pipeline treats PASS-only as baseline hygiene -
# see common.R::pass_filter_maf and convert_to_maf.sh's own bcftools
# pre-filter). This also substantially cuts the memory footprint of the
# downstream read_vcfs_as_granges()/mut_matrix() call for WGS-scale cohorts.
rule filter_signature_vcf_pass:
    input:
        vcf=signature_caller_vcf_for,
    output:
        vcf=f"{OUTDIR}/signatures/pass_vcfs/{{contrast}}.vcf.gz",
    log:
        f"{OUTDIR}/logs/signatures/pass_filter/{{contrast}}.log",
    resources:
        tmpdir=R_TMPDIR,
        mem_mb=1000,
    conda:
        "../envs/vcf2maf.yaml"
    shell:
        "bcftools view -f 'PASS,.' -O z -o {output.vcf} {input.vcf} > {log} 2>&1 "
        "&& bcftools index -t {output.vcf} >> {log} 2>&1"


def signature_pass_vcfs(wildcards):
    return [f"{OUTDIR}/signatures/pass_vcfs/{c}.vcf.gz" for c in signature_contrast_names(wildcards)]


rule mutational_signatures:
    input:
        mutect2_vcfs=signature_pass_vcfs,
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
    resources:
        tmpdir=R_TMPDIR,
        mem_mb=16000,
    conda:
        "../envs/r_env.yaml"
    script:
        "../scripts/mutational_signatures.R"
