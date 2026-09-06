#!/usr/bin/env bash
# One-off workaround for the Tumor_Sample_Barcode/Matched_Norm_Sample_Barcode
# mismatch fixed at the source in discover_contrasts.py::resolve_tumor_normal()
# (a Mutect2 VCF's own genotype-column names, e.g.
# "BC139_bn_BC139_naive_1", were leaking into the MAF barcode instead of the
# clean "<tumor>_vs_<normal>" contrast token Strelka already used - breaking
# the cross-caller consensus join in consolidate_mafs.R).
#
# Re-running vcf2maf/VEP is expensive; this rewrites the barcode columns of
# the ALREADY-PRODUCED per-caller MAFs (results/maf/<contrast>.<caller>.maf.gz)
# in place, deriving the correct tokens from each file's own name - matching
# what the fixed discover_contrasts.py would have produced - without
# touching any other MAF column. Once you can afford to rerun the full
# pipeline from vcf2maf onward, this script (and its output) becomes
# unnecessary; nothing here permanently changes the pipeline's *code path*,
# only these specific already-existing files.
#
# Usage:
#   fix_maf_barcodes.sh <maf_dir> [caller_suffix ...]
#
# <maf_dir>       directory containing "<contrast>.<caller>.maf.gz" files
#                  (default: results/maf, i.e. Snakemake's OUTDIR/maf)
# caller_suffix   one or more caller names as they appear in the filename
#                  (default: mutect2 strelka_snvs strelka_indels)
#
# For each matching file:
#   - contrast   = filename with ".<caller>.maf.gz" stripped
#   - tumor_tok  = contrast split on "_vs_", first half (or the whole
#                  contrast if no "_vs_", i.e. tumor-only)
#   - normal_tok = contrast split on "_vs_", second half (empty if tumor-only)
# Every data row's Tumor_Sample_Barcode is overwritten with tumor_tok, and
# Matched_Norm_Sample_Barcode (if the column exists and normal_tok is
# non-empty) with normal_tok. Column positions are looked up by header name,
# not assumed, so this is robust to any vcf2maf.pl column-order/version.
#
# Original files are backed up alongside as "<file>.bak" before being
# overwritten; delete the .baks once you've spot-checked the result.
set -euo pipefail

MAF_DIR="${1:-results/maf}"
shift || true
CALLERS=("$@")
[[ ${#CALLERS[@]} -gt 0 ]] || CALLERS=(mutect2 strelka_snvs strelka_indels)

[[ -d "$MAF_DIR" ]] || { echo "ERROR: MAF directory not found: $MAF_DIR"; exit 1; }
command -v gzip >/dev/null 2>&1 || { echo "ERROR: gzip not found"; exit 1; }

n_fixed=0
n_skipped=0

for caller in "${CALLERS[@]}"; do
  suffix=".${caller}.maf.gz"
  while IFS= read -r -d '' f; do
    base="$(basename "$f")"
    contrast="${base%"$suffix"}"
    if [[ "$contrast" == "$base" ]]; then
      continue  # suffix didn't match this file
    fi

    if [[ "$contrast" == *_vs_* ]]; then
      tumor_tok="${contrast%%_vs_*}"
      normal_tok="${contrast#*_vs_}"
    else
      tumor_tok="$contrast"
      normal_tok=""
    fi

    tmp="$(mktemp)"
    zcat "$f" > "$tmp"

    # Header is the first non-"##"/"#version" line starting with "Hugo_Symbol".
    header_line_num="$(grep -n '^Hugo_Symbol' "$tmp" | head -1 | cut -d: -f1 || true)"
    if [[ -z "$header_line_num" ]]; then
      echo "WARNING: no MAF header found in $f, skipping"
      rm -f "$tmp"
      n_skipped=$((n_skipped + 1))
      continue
    fi

    out="$(mktemp)"
    awk -F'\t' -v OFS='\t' -v hln="$header_line_num" -v ttok="$tumor_tok" -v ntok="$normal_tok" '
      NR < hln { print; next }
      NR == hln {
        for (i = 1; i <= NF; i++) {
          if ($i == "Tumor_Sample_Barcode") tcol = i
          if ($i == "Matched_Norm_Sample_Barcode") ncol = i
        }
        print
        next
      }
      {
        if (tcol) $tcol = ttok
        if (ncol && ntok != "") $ncol = ntok
        print
      }
    ' "$tmp" > "$out"

    cp "$f" "${f}.bak"
    gzip -c "$out" > "$f"
    rm -f "$tmp" "$out"

    echo "Fixed $base: Tumor_Sample_Barcode -> '$tumor_tok'$( [[ -n "$normal_tok" ]] && echo ", Matched_Norm_Sample_Barcode -> '$normal_tok'" )"
    n_fixed=$((n_fixed + 1))
  done < <(find "$MAF_DIR" -maxdepth 1 -name "*${suffix}" -print0)
done

echo "Done: $n_fixed file(s) fixed, $n_skipped skipped (no header found)."
echo "Backups written as <file>.maf.gz.bak next to each fixed file."
echo "Next: re-run Snakemake from consolidate_contrast onward, e.g.:"
echo "  snakemake -s workflow/Snakefile --configfile <your_config>.yaml --use-conda --cores 8 \\"
echo "    --resources mem_mb=55000 --forcerun consolidate_contrast build_cohort_union cohort_oncoplots driver_candidates actionability_report"
