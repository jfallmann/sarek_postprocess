rule cohort_oncoplots:
    input:
        union_rds=f"{OUTDIR}/maf_consolidated/cohort/UnionMAF.rds.gz",
    output:
        done=f"{OUTDIR}/oncoplots/cohort/.done",
    params:
        out_dir=f"{OUTDIR}/oncoplots/cohort",
    log:
        f"{OUTDIR}/logs/cohort_oncoplots.log",
    conda:
        "../envs/r_env.yaml"
    script:
        "../scripts/cohort_oncoplots.R"
