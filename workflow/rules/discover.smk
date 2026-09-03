checkpoint discover_contrasts:
    input:
        samplesheet=config["samplesheet"],
    output:
        tsv=f"{OUTDIR}/contrasts.tsv",
        cnv_tsv=f"{OUTDIR}/cnv_contrasts.tsv",
    params:
        sarek_outdir=config["sarek_outdir"],
        callers=",".join(config["callers"]),
        sample_subset_regex=config.get("sample_subset_regex", ""),
        mutect2_dir=config.get("tool_output_dirs", {}).get("mutect2", ""),
        strelka_dir=config.get("tool_output_dirs", {}).get("strelka", ""),
        manta_dir=config.get("tool_output_dirs", {}).get("manta", ""),
        cnvkit_dir=config.get("tool_output_dirs", {}).get("cnvkit", ""),
    log:
        f"{OUTDIR}/logs/discover_contrasts.log",
    conda:
        "../envs/vcf2maf.yaml"
    shell:
        "python3 {workflow.basedir}/scripts/discover_contrasts.py "
        "--sarek-outdir {params.sarek_outdir} "
        "--samplesheet {input.samplesheet} "
        "--callers {params.callers} --include-manta "
        "--sample-subset-regex '{params.sample_subset_regex}' "
        "--mutect2-dir '{params.mutect2_dir}' "
        "--strelka-dir '{params.strelka_dir}' "
        "--manta-dir '{params.manta_dir}' "
        "--cnvkit-dir '{params.cnvkit_dir}' "
        "--out {output.tsv} --cnv-out {output.cnv_tsv} > {log} 2>&1"
