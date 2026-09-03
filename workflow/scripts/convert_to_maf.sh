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
#                      <ref_fasta> <build> <vep_path> <vep_data> <cache_version> <vcf2maf_pl>
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

command -v "$V2M" >/dev/null 2>&1 || { echo "ERROR: vcf2maf.pl not found (V2M=$V2M)"; exit 1; }
[[ -f "$REF_FASTA" ]] || { echo "ERROR: REF_FASTA not found: $REF_FASTA"; exit 1; }
[[ -f "$IN_VCF" ]] || { echo "ERROR: input VCF not found: $IN_VCF"; exit 1; }
[[ -n "$TUMOR_ID" ]] || { echo "ERROR: TUMOR_ID is empty - discover_contrasts.py could not resolve it"; exit 1; }

mkdir -p "$(dirname "$OUT_MAF")"
tmpdir="$(mktemp -d)"
trap 'rm -rf "$tmpdir"' EXIT

# vcf2maf.pl requires an uncompressed VCF.
in_vcf="$IN_VCF"
if [[ "$IN_VCF" == *.vcf.gz ]]; then
  in_vcf="$tmpdir/$(basename "${IN_VCF%.gz}")"
  gunzip -c "$IN_VCF" > "$in_vcf"
fi

cmd=(perl "$V2M" --input-vcf "$in_vcf" --output-maf "$OUT_MAF" \
     --tumor-id "$TUMOR_ID" \
     --ref-fasta "$REF_FASTA" --ncbi-build "$BUILD" \
     --vep-data "$VEP_DATA" --species homo_sapiens --cache-version "$VEP_CACHE" \
     --vep-path "$VEP_PATH")

if [[ -n "$NORMAL_ID" && "$NORMAL_ID" != "NA" ]]; then
  cmd+=(--normal-id "$NORMAL_ID")
fi

echo "[convert_to_maf] ${cmd[*]}"
"${cmd[@]}"
