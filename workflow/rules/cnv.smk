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
        min_weight=config.get("cnv", {}).get("min_weight", 0.5),
        min_probes=config.get("cnv", {}).get("min_probes", 0),
    log:
        f"{OUTDIR}/logs/cnvkit/{{contrast}}.log",
    resources:
        tmpdir=R_TMPDIR,
        mem_mb=2000,
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
        # Two separate heatmaps instead of one cohort-wide plot: "vs_reference"
        # (bare tumor-only contrasts, called against a generic/pooled
        # reference) and "vs_contrast" (matched tumor/normal "*_vs_Bulk_sensitive"
        # contrasts) have different copy-number backgrounds and are not
        # meaningfully comparable side by side - mixing them was also what
        # made the single heatmap too large to open/read. matrix_tsv/calls_tsv
        # stay cohort-wide (driver_candidates.R needs all evidence regardless
        # of contrast type).
        heatmap_vs_reference_pdf=f"{OUTDIR}/cnv/cohort/cnv_heatmap.vs_reference.pdf",
        heatmap_vs_contrast_pdf=f"{OUTDIR}/cnv/cohort/cnv_heatmap.vs_contrast.pdf",
    params:
        common_r=f"{workflow.basedir}/scripts/common.R",
        contrasts=cnv_contrasts,
        gene_panel_csv=config.get("cnv", {}).get("heatmap_gene_panel_csv", ""),
        default_gene_panel_csv=config["gene_panel_csv"],
        min_log2_gain=config.get("cnv", {}).get("min_log2_gain", 0.3),
        min_log2_loss=config.get("cnv", {}).get("min_log2_loss", -0.3),
    log:
        f"{OUTDIR}/logs/cohort_cnv_matrix.log",
    resources:
        tmpdir=R_TMPDIR,
        mem_mb=2000,
    conda:
        "../envs/r_env.yaml"
    script:
        "../scripts/cohort_cnv_matrix.R"
