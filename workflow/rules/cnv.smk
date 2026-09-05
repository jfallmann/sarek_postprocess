def get_cnv_row(wildcards):
    df = get_cnv_df()
    row = df[df.contrast == wildcards.contrast]
    if row.empty:
        raise WorkflowError(f"No discovered CNVkit segment file for contrast={wildcards.contrast}")
    return row.iloc[0]


rule summarize_cnvkit:
    input:
        cns=lambda wc: get_cnv_row(wc)["cns_path"],
    output:
        gene_tsv=f"{OUTDIR}/cnv/{{contrast}}.cnvkit_genes.tsv.gz",
    params:
        common_r=f"{workflow.basedir}/scripts/common.R",
    log:
        f"{OUTDIR}/logs/cnvkit/{{contrast}}.log",
    conda:
        "../envs/r_env.yaml"
    script:
        "../scripts/summarize_cnvkit.R"


def cnv_contrasts(wildcards):
    return sorted(get_cnv_df().contrast.unique().tolist())


def cnv_gene_tsvs(wildcards):
    return [f"{OUTDIR}/cnv/{c}.cnvkit_genes.tsv.gz" for c in cnv_contrasts(wildcards)]


rule cohort_cnv_matrix:
    input:
        gene_tsvs=cnv_gene_tsvs,
    output:
        matrix_tsv=f"{OUTDIR}/cnv/cohort/cnv_matrix.tsv.gz",
        calls_tsv=f"{OUTDIR}/cnv/cohort/cnv_calls.tsv.gz",
        heatmap_pdf=f"{OUTDIR}/cnv/cohort/cnv_heatmap.pdf",
    params:
        common_r=f"{workflow.basedir}/scripts/common.R",
        contrasts=cnv_contrasts,
        gene_panel_csv=config.get("cnv", {}).get("heatmap_gene_panel_csv", ""),
        default_gene_panel_csv=config["gene_panel_csv"],
        min_log2_gain=config.get("cnv", {}).get("min_log2_gain", 0.3),
        min_log2_loss=config.get("cnv", {}).get("min_log2_loss", -0.3),
    log:
        f"{OUTDIR}/logs/cohort_cnv_matrix.log",
    conda:
        "../envs/r_env.yaml"
    script:
        "../scripts/cohort_cnv_matrix.R"
