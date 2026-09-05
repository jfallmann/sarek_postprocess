def get_contrast_caller_row(wildcards):
    df = get_contrasts_df()
    row = df[(df.contrast == wildcards.contrast) & (df.caller == wildcards.caller)]
    if row.empty:
        raise WorkflowError(
            f"No discovered VCF for contrast={wildcards.contrast} caller={wildcards.caller}"
        )
    return row.iloc[0]


rule vcf2maf:
    input:
        vcf=lambda wc: get_contrast_caller_row(wc)["vcf_path"],
    output:
        maf=f"{OUTDIR}/maf/{{contrast}}.{{caller}}.maf.gz",
    params:
        # Pre-filter the VCF to FILTER=='PASS'|'.' records before handing it
        # to VEP, mirroring common.R's pass_filter_maf() so the final
        # union/consensus MAFs are unchanged - just skip the filter for any
        # caller explicitly listed as False in qual_override_by_caller (e.g.
        # Strelka indels, which have no reliable FILTER/QUAL).
        pass_filter=lambda wc: (
            "0"
            if config.get("filtering", {}).get("qual_override_by_caller", {}).get(wc.caller) is False
            else "1"
        ),
        tumor_id=lambda wc: get_contrast_caller_row(wc)["tumor_id"],
        normal_id=lambda wc: (get_contrast_caller_row(wc)["normal_id"] or "NA"),
        # literal sample-column names inside the VCF's genotype fields; only
        # differ from tumor_id/normal_id for callers with generic labels
        # (e.g. Strelka2's "TUMOR"/"NORMAL"). Fall back to tumor_id/normal_id
        # if an older contrasts.tsv without these columns is in use.
        vcf_tumor_id=lambda wc: (get_contrast_caller_row(wc).get("vcf_tumor_id") or get_contrast_caller_row(wc)["tumor_id"]),
        vcf_normal_id=lambda wc: (get_contrast_caller_row(wc).get("vcf_normal_id") or get_contrast_caller_row(wc)["normal_id"] or "NA"),
        ref_fasta=config["ref_fasta"],
        build=config["genome_build"],
        vep_path=config["vcf2maf"]["vep_path"],
        vep_data=config["vcf2maf"]["vep_data"],
        cache_version=config["vcf2maf"]["cache_version"],
        vcf2maf_pl=config["vcf2maf"]["vcf2maf_pl"],
        tmp_dir=config.get("tmp_dir", "") or "-",
    log:
        f"{OUTDIR}/logs/vcf2maf/{{contrast}}.{{caller}}.log",
    threads: config.get("threads_default", 4)
    conda:
        "../envs/vcf2maf.yaml"
    shell:
        "bash {workflow.basedir}/scripts/convert_to_maf.sh "
        "{input.vcf} {output.maf} {params.tumor_id} {params.normal_id} "
        "{params.ref_fasta} {params.build} {params.vep_path} {params.vep_data} "
        "{params.cache_version} {params.vcf2maf_pl} "
        "{params.vcf_tumor_id} {params.vcf_normal_id} {threads} {params.pass_filter} "
        "{params.tmp_dir} > {log} 2>&1"
