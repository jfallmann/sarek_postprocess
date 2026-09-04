#!/usr/bin/env bash
# Convert a single Sarek-annotated VCF to MAF via vcf2maf.pl.
#
# Adapted from the GENOMICS/vcf_to_maf.sh directory-scanning script: here the
# contrast/caller/tumor-id/normal-id have already been resolved explicitly by
# discover_contrasts.py, so this script does no filename heuristics itself -
# it only re-runs VEP inside vcf2maf.pl and writes one MAF.
#
# Usage:
#   convert_to_maf.sh <in_vcf> <out_maf> <tumor_id> <normal_id_or_NA> \
#                      <ref_fasta> <build> <vep_path> <vep_data> <cache_version> <vcf2maf_pl> \
#                      [<vcf_tumor_id_or_NA> [<vcf_normal_id_or_NA> [<vep_forks> [<pass_filter>]]]]
#
# vcf_tumor_id/vcf_normal_id are the literal sample-column names inside the
# VCF's genotype fields; they only differ from tumor_id/normal_id when a
# caller hard-codes generic labels there (e.g. Strelka2's "TUMOR"/"NORMAL"
# instead of the real sample IDs). If omitted, they default to tumor_id/
# normal_id (i.e. no relabeling, matching vcf2maf.pl's own default).
#
# vep_forks sets vcf2maf.pl's --vep-forks (parallel VEP annotation
# processes); it should match the Snakemake rule's `threads` so VEP actually
# uses all cores allocated to the job. Defaults to 4 (vcf2maf.pl's own
# default) if omitted.
#
# pass_filter (0/1, default 1): if 1, drop non-PASS/non-"." records with
# `bcftools view -f PASS,.` before running VEP, mirroring common.R's
# pass_filter_maf() so the variants that reach VEP match what the
# union/consensus MAFs keep anyway - this just cuts VEP's workload instead
# of discarding those rows only after paying for their annotation. The
# untouched original VCF is never modified, only a filtered temp copy is
# used here, so nothing upstream is lost by setting this to 1.
#
# out_maf is written gzip-compressed regardless of pass_filter.
set -euo pipefail

IN_VCF="$1"
OUT_MAF="$2"
TUMOR_ID="$3"
NORMAL_ID="$4"          # literal "NA" if tumor-only
REF_FASTA="$5"
BUILD="$6"
VEP_PATH="$7"
VEP_DATA="$8"
VEP_CACHE="$9"
V2M="${10}"
VCF_TUMOR_ID="${11:-$TUMOR_ID}"
VCF_NORMAL_ID="${12:-$NORMAL_ID}"
VEP_FORKS="${13:-4}"
PASS_FILTER="${14:-1}"
[[ -n "$VCF_TUMOR_ID" && "$VCF_TUMOR_ID" != "NA" ]] || VCF_TUMOR_ID="$TUMOR_ID"
[[ -n "$VCF_NORMAL_ID" && "$VCF_NORMAL_ID" != "NA" ]] || VCF_NORMAL_ID="$NORMAL_ID"
[[ -n "$VEP_FORKS" ]] || VEP_FORKS=4
[[ -n "$PASS_FILTER" ]] || PASS_FILTER=1

V2M_RESOLVED="$(command -v "$V2M" 2>/dev/null || true)"
[[ -n "$V2M_RESOLVED" ]] || { echo "ERROR: vcf2maf.pl not found (V2M=$V2M)"; exit 1; }
# perl does not search PATH for its script argument (only the shell does),
# so pass the fully resolved path rather than the bare/PATH-relative name.
V2M="$V2M_RESOLVED"

# vcf2maf.pl's --vep-path must be the *directory* containing the vep
# executable (it internally execs "<vep_path>/vep"), not a bare command
# name resolved via PATH.
if [[ ! -d "$VEP_PATH" ]]; then
  vep_bin="$(command -v "$VEP_PATH" 2>/dev/null || true)"
  [[ -n "$vep_bin" ]] || { echo "ERROR: vep executable not found (vep_path=$VEP_PATH)"; exit 1; }
  VEP_PATH="$(dirname "$vep_bin")"
fi
# VEP's --dir/vcf2maf's --vep-data must be the cache *root* (the directory
# that directly contains "homo_sapiens/<cache_version>_<build>/"), not that
# species/version directory itself - VEP appends "homo_sapiens/..." on top
# of whatever --dir it is given. If vep_data is configured as the full
# nested path (a common mistake when copying the path Sarek's cache actually
# lives under), auto-correct it back to the root two levels up.
vep_data_parent="$(dirname "$VEP_DATA")"
if [[ "$(basename "$vep_data_parent")" == "homo_sapiens" ]]; then
  echo "[convert_to_maf] vep_data '$VEP_DATA' looks like .../homo_sapiens/<version>_<build>; using its cache root '$(dirname "$vep_data_parent")' instead"
  VEP_DATA="$(dirname "$vep_data_parent")"
fi

[[ -f "$REF_FASTA" ]] || { echo "ERROR: REF_FASTA not found: $REF_FASTA"; exit 1; }
[[ -f "$IN_VCF" ]] || { echo "ERROR: input VCF not found: $IN_VCF"; exit 1; }
[[ -n "$TUMOR_ID" ]] || { echo "ERROR: TUMOR_ID is empty - discover_contrasts.py could not resolve it"; exit 1; }

mkdir -p "$(dirname "$OUT_MAF")"
tmpdir="$(mktemp -d)"
trap 'rm -rf "$tmpdir"' EXIT

# vcf2maf requires an uncompressed, samtools-faidx-indexed reference FASTA.
ref_fasta="$REF_FASTA"
if [[ "$REF_FASTA" == *.gz ]]; then
  if [[ -f "${REF_FASTA%.gz}" ]]; then
    # iGenomes bundles commonly ship both the .fasta.gz and an uncompressed
    # .fasta side by side - prefer that instead of re-decompressing per job.
    ref_fasta="${REF_FASTA%.gz}"
  else
    ref_fasta="$tmpdir/$(basename "${REF_FASTA%.gz}")"
    gunzip -c "$REF_FASTA" > "$ref_fasta"
  fi
fi
if [[ ! -f "${ref_fasta}.fai" ]]; then
  command -v samtools >/dev/null 2>&1 || { echo "ERROR: samtools not found, needed to index REF_FASTA"; exit 1; }
  if ! samtools faidx "$ref_fasta" 2>/dev/null; then
    # Reference directory is likely read-only (e.g. a shared iGenomes mount);
    # copy locally so we can index it instead.
    local_fasta="$tmpdir/$(basename "$ref_fasta")"
    cp "$ref_fasta" "$local_fasta"
    samtools faidx "$local_fasta"
    ref_fasta="$local_fasta"
  fi
fi
REF_FASTA="$ref_fasta"

# vcf2maf.pl requires an uncompressed VCF.
in_vcf="$IN_VCF"
if [[ "$IN_VCF" == *.vcf.gz ]]; then
  in_vcf="$tmpdir/$(basename "${IN_VCF%.gz}")"
  gunzip -c "$IN_VCF" > "$in_vcf"
fi

if [[ "$PASS_FILTER" == "1" ]]; then
  command -v bcftools >/dev/null 2>&1 || { echo "ERROR: bcftools not found, needed for PASS pre-filtering"; exit 1; }
  n_before="$(grep -vc '^#' "$in_vcf" || true)"
  filtered_vcf="$tmpdir/filtered.vcf"
  bcftools view -f "PASS,." -O v -o "$filtered_vcf" "$in_vcf"
  n_after="$(grep -vc '^#' "$filtered_vcf" || true)"
  echo "[convert_to_maf] PASS pre-filter: $n_before -> $n_after records"
  in_vcf="$filtered_vcf"
fi

# vcf2maf.pl needs an uncompressed output path; the final MAF is
# gzip-compressed afterwards.
tmp_maf="$tmpdir/output.maf"

cmd=(perl "$V2M" --input-vcf "$in_vcf" --output-maf "$tmp_maf" \
     --tumor-id "$TUMOR_ID" --vcf-tumor-id "$VCF_TUMOR_ID" \
     --ref-fasta "$REF_FASTA" --ncbi-build "$BUILD" \
     --vep-data "$VEP_DATA" --species homo_sapiens --cache-version "$VEP_CACHE" \
     --vep-path "$VEP_PATH" --vep-forks "$VEP_FORKS")

if [[ -n "$NORMAL_ID" && "$NORMAL_ID" != "NA" ]]; then
  cmd+=(--normal-id "$NORMAL_ID" --vcf-normal-id "$VCF_NORMAL_ID")
fi

echo "[convert_to_maf] ${cmd[*]}"
"${cmd[@]}"

gzip -c "$tmp_maf" > "$OUT_MAF"
