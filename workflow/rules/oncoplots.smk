rule cohort_oncoplots:
    input:
        union_rds=f"{OUTDIR}/maf_consolidated/cohort/{{vcset}}.UnionMAF.rds.gz",
    output:
        done=f"{OUTDIR}/oncoplots/cohort/{{vcset}}/.done",
    params:
        out_dir=f"{OUTDIR}/oncoplots/cohort/{{vcset}}",
    log:
        f"{OUTDIR}/logs/cohort_oncoplots/{{vcset}}.log",
    conda:
        "../envs/r_env.yaml"
    script:
        "../scripts/cohort_oncoplots.R"
