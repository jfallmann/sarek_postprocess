def snv_indel_callers(df):
    return df[~df.caller.str.startswith("manta")]


def contrast_caller_mafs(wildcards):
    df = snv_indel_callers(get_contrasts_df())
    sub = df[df.contrast == wildcards.contrast]
    return sub["caller"].tolist()


def contrast_maf_paths(wildcards):
    callers = contrast_caller_mafs(wildcards)
    return [f"{OUTDIR}/maf/{wildcards.contrast}.{c}.maf.gz" for c in callers]


# Every consensus/union/oncoplot/driver output is produced in two parallel
# tracks, selected by Variant_Classification set (see
# common.R::resolve_vc_nonsyn):
#   - "custom": config-driven (filtering.vc_nonsyn_custom), less stringent,
#     includes splice-region/UTR/intron/etc. by default - the recommended
#     track for WGS.
#   - "protein_coding": maftools' own stringent built-in default
#     (missense/nonsense/frameshift/splice-site/in-frame indel only), not
#     configurable, kept as a reproducible strict baseline.
VCSETS = ["custom", "protein_coding"]


wildcard_constraints:
    vcset="|".join(VCSETS),


rule consolidate_contrast:
    input:
        mafs=contrast_maf_paths,
    output:
        union_maf=f"{OUTDIR}/maf_consolidated/{{vcset}}/{{contrast}}.union.tsv.gz",
        consensus_maf=f"{OUTDIR}/maf_consolidated/{{vcset}}/{{contrast}}.consensus.tsv.gz",
    params:
        common_r=f"{workflow.basedir}/scripts/common.R",
        callers=contrast_caller_mafs,
        filtering=config["filtering"],
    log:
        f"{OUTDIR}/logs/consolidate/{{vcset}}/{{contrast}}.log",
    conda:
        "../envs/r_env.yaml"
    script:
        "../scripts/consolidate_mafs.R"


def all_contrasts(wildcards):
    df = snv_indel_callers(get_contrasts_df())
    return sorted(df.contrast.unique().tolist())


def cohort_union_tsvs(wildcards):
    return [f"{OUTDIR}/maf_consolidated/{wildcards.vcset}/{c}.union.tsv.gz" for c in all_contrasts(wildcards)]


rule build_cohort_union:
    input:
        union_tsvs=cohort_union_tsvs,
    output:
        maf_list_rds=f"{OUTDIR}/maf_consolidated/cohort/{{vcset}}.MAFlist.rds.gz",
        union_rds=f"{OUTDIR}/maf_consolidated/cohort/{{vcset}}.UnionMAF.rds.gz",
    params:
        common_r=f"{workflow.basedir}/scripts/common.R",
        contrasts=all_contrasts,
        filtering=config["filtering"],
        sample_subset_regex=config.get("sample_subset_regex", ""),
    log:
        f"{OUTDIR}/logs/build_cohort_union/{{vcset}}.log",
    conda:
        "../envs/r_env.yaml"
    script:
        "../scripts/build_cohort_union.R"
