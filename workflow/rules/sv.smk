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
        sv_tsv=f"{OUTDIR}/sv/{{contrast}}.manta_sv.tsv.gz",
    params:
        common_r=f"{workflow.basedir}/scripts/common.R",
    log:
        f"{OUTDIR}/logs/manta/{{contrast}}.log",
    resources:
        tmpdir=R_TMPDIR,
        mem_mb=2000,
    conda:
        "../envs/r_env.yaml"
    script:
        "../scripts/summarize_manta_sv.R"


def manta_contrasts(wildcards):
    df = get_contrasts_df()
    return sorted(df[df.caller.str.startswith("manta")].contrast.unique().tolist())


def all_manta_sv_tsvs(wildcards):
    return [f"{OUTDIR}/sv/{c}.manta_sv.tsv.gz" for c in manta_contrasts(wildcards)]


rule cohort_sv_summary:
    input:
        sv_tsvs=all_manta_sv_tsvs,
    output:
        burden_tsv=f"{OUTDIR}/sv/cohort/sv_burden.tsv.gz",
        gene_calls_tsv=f"{OUTDIR}/sv/cohort/sv_gene_calls.tsv.gz",
        recurrent_tsv=f"{OUTDIR}/sv/cohort/sv_recurrent_genes.tsv.gz",
        burden_barplot=f"{OUTDIR}/sv/cohort/sv_burden_barplot.pdf",
        large_sv_tsv=f"{OUTDIR}/sv/cohort/sv_large_events.tsv.gz",
        large_sv_barplot=f"{OUTDIR}/sv/cohort/sv_large_events_barplot.pdf",
    params:
        common_r=f"{workflow.basedir}/scripts/common.R",
        contrasts=manta_contrasts,
        pass_only=config.get("sv_analysis", {}).get("pass_only", True),
        gene_field_index=config.get("sv_analysis", {}).get("gene_field_index", 4),
        gene_panel_csv=lambda wc: (
            config.get("sv_analysis", {}).get("gene_panel_csv", "") or config["gene_panel_csv"]
        ),
        large_sv_min_bp=config.get("sv_analysis", {}).get("large_sv_min_bp", 1000000),
    log:
        f"{OUTDIR}/logs/cohort_sv_summary.log",
    resources:
        tmpdir=R_TMPDIR,
        mem_mb=2000,
    conda:
        "../envs/r_env.yaml"
    script:
        "../scripts/cohort_sv_summary.R"
