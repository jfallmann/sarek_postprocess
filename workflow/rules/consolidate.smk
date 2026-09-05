def snv_indel_callers(df):
    return df[~df.caller.str.startswith("manta")]


def contrast_caller_mafs(wildcards):
    df = snv_indel_callers(get_contrasts_df())
    sub = df[df.contrast == wildcards.contrast]
    return sub["caller"].tolist()


def contrast_maf_paths(wildcards):
    callers = contrast_caller_mafs(wildcards)
    return [f"{OUTDIR}/maf/{wildcards.contrast}.{c}.maf.gz" for c in callers]


rule consolidate_contrast:
    input:
        mafs=contrast_maf_paths,
    output:
        union_maf=f"{OUTDIR}/maf_consolidated/{{contrast}}.union.tsv.gz",
        consensus_maf=f"{OUTDIR}/maf_consolidated/{{contrast}}.consensus.tsv.gz",
    params:
        common_r=f"{workflow.basedir}/scripts/common.R",
        callers=contrast_caller_mafs,
        filtering=config["filtering"],
    log:
        f"{OUTDIR}/logs/consolidate/{{contrast}}.log",
    conda:
        "../envs/r_env.yaml"
    script:
        "../scripts/consolidate_mafs.R"


def all_contrasts(wildcards):
    df = snv_indel_callers(get_contrasts_df())
    return sorted(df.contrast.unique().tolist())


def cohort_union_tsvs(wildcards):
    return [f"{OUTDIR}/maf_consolidated/{c}.union.tsv.gz" for c in all_contrasts(wildcards)]


rule build_cohort_union:
    input:
        union_tsvs=cohort_union_tsvs,
    output:
        maf_list_rds=f"{OUTDIR}/maf_consolidated/cohort/MAFlist.rds.gz",
        union_rds=f"{OUTDIR}/maf_consolidated/cohort/UnionMAF.rds.gz",
    params:
        contrasts=all_contrasts,
        sample_subset_regex=config.get("sample_subset_regex", ""),
    log:
        f"{OUTDIR}/logs/build_cohort_union.log",
    conda:
        "../envs/r_env.yaml"
    script:
        "../scripts/build_cohort_union.R"
