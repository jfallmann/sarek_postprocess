checkpoint discover_contrasts:
    input:
        samplesheet=config["samplesheet"],
    output:
        tsv=f"{OUTDIR}/contrasts.tsv",
    params:
        sarek_outdir=config["sarek_outdir"],
        callers=",".join(config["callers"]),
        sample_subset_regex=config.get("sample_subset_regex", ""),
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
        "--out {output.tsv} > {log} 2>&1"
