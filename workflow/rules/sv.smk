def get_manta_row(wildcards):
    df = get_contrasts_df()
    sub = df[(df.contrast == wildcards.contrast) & (df.caller.str.startswith("manta"))]
    if sub.empty:
        raise WorkflowError(f"No Manta VCF discovered for contrast={wildcards.contrast}")
    # Prefer somatic calls over diploid/candidate if several are present.
    if (sub.caller == "manta_somatic").any():
        return sub[sub.caller == "manta_somatic"].iloc[0]
    return sub.iloc[0]


rule summarize_manta_sv:
    input:
        vcf=lambda wc: get_manta_row(wc)["vcf_path"],
    output:
        sv_tsv=f"{OUTDIR}/sv/{{contrast}}.manta_sv.tsv",
    log:
        f"{OUTDIR}/logs/manta/{{contrast}}.log",
    conda:
        "../envs/r_env.yaml"
    script:
        "../scripts/summarize_manta_sv.R"


def manta_contrasts(wildcards):
    df = get_contrasts_df()
    return sorted(df[df.caller.str.startswith("manta")].contrast.unique().tolist())


def all_manta_sv_tsvs(wildcards):
    return [f"{OUTDIR}/sv/{c}.manta_sv.tsv" for c in manta_contrasts(wildcards)]
