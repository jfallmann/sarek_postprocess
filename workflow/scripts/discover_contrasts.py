#!/usr/bin/env python3
"""Discover Sarek somatic-variant contrasts from an actual Sarek --outdir.

Unlike re-deriving tumor-only/tumor-normal pairing from the samplesheet
(status column), this walks Sarek's own output directory structure, which
already encodes the pairing Sarek decided on:

    annotation/<caller>/<contrast>/<contrast>.<caller>....vcf.gz
    variant_calling/<caller>/<contrast>/...                (manta, fallback)

where <contrast> is either "<tumor>_vs_<normal>" (tumor-normal) or
"<sample>" (tumor-only / germline-only).

The samplesheet is read only to attach sample metadata (patient, sex,
status) to the discovered contrasts for later clinical annotation - it is
never used to re-derive pairing.

Any caller directory can be overridden manually (e.g. you reran Mutect2 with
different settings, or ran CNVkit yourself outside Sarek entirely) via
--<caller>-dir; each override must still contain one subdirectory per
contrast, matching Sarek's own <contrast>/<contrast>.<caller>....ext layout.

CNVkit is not part of Sarek's own --tools and is discovered separately
(--cnvkit-dir, default "<sarek-outdir>/cnvkit" if not overridden), keyed to
the same contrasts already found for the VCF-based callers, and written to
a second output file (--cnv-out).

Usage:
    discover_contrasts.py --sarek-outdir DIR --samplesheet CSV \
        --callers mutect2,strelka --out contrasts.tsv.gz --cnv-out cnv_contrasts.tsv.gz \
        [--mutect2-dir DIR] [--strelka-dir DIR] [--manta-dir DIR] [--cnvkit-dir DIR] \
        [--sample-subset-regex REGEX]
"""
import argparse
import csv
import gzip
import os
import re
import subprocess
import sys
from pathlib import Path

CALLER_FILE_PATTERNS = {
    "mutect2": [
        (re.compile(r"\.mutect2\.filtered.*\.vcf\.gz$"), "mutect2"),
        (re.compile(r"\.mutect2\.vcf\.gz$"), "mutect2"),
    ],
    "strelka": [
        (re.compile(r"\.strelka\.somatic_snvs.*\.vcf\.gz$"), "strelka_snvs"),
        (re.compile(r"\.strelka\.somatic_indels.*\.vcf\.gz$"), "strelka_indels"),
        (re.compile(r"\.strelka\.variants.*\.vcf\.gz$"), "strelka_variants"),
    ],
    "manta": [
        (re.compile(r"\.somaticSV.*\.vcf\.gz$"), "manta_somatic"),
        (re.compile(r"\.diploidSV.*\.vcf\.gz$"), "manta_diploid"),
        (re.compile(r"\.candidateSV.*\.vcf\.gz$"), "manta_candidate"),
    ],
}


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--sarek-outdir", required=True)
    p.add_argument("--samplesheet", required=True)
    p.add_argument("--callers", default="mutect2,strelka")
    p.add_argument("--include-manta", action="store_true")
    p.add_argument("--sample-subset-regex", default="")
    p.add_argument("--out", required=True)
    p.add_argument("--mutect2-dir", default="", help="Override auto-discovery of the mutect2 caller directory")
    p.add_argument("--strelka-dir", default="", help="Override auto-discovery of the strelka caller directory")
    p.add_argument("--manta-dir", default="", help="Override auto-discovery of the manta caller directory")
    p.add_argument("--cnvkit-dir", default="", help="CNVkit output directory (not part of Sarek); "
                                                      "defaults to <sarek-outdir>/cnvkit if unset")
    p.add_argument("--cnv-out", default="", help="Output TSV path for discovered CNVkit segment files; "
                                                   "skipped entirely if not given")
    return p.parse_args()


def read_samplesheet(path):
    """Return {sample_name: {patient, sex, status}} from a Sarek samplesheet."""
    meta = {}
    if not os.path.isfile(path):
        print(f"WARNING: samplesheet not found at {path}; metadata will be empty", file=sys.stderr)
        return meta
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            sample = row.get("sample")
            if not sample:
                continue
            meta[sample] = {
                "patient": row.get("patient", ""),
                "sex": row.get("sex", ""),
                "status": row.get("status", ""),
            }
    return meta


def vcf_samples(vcf_path):
    """Return the sample columns of a VCF via bcftools."""
    try:
        out = subprocess.run(
            ["bcftools", "query", "-l", str(vcf_path)],
            check=True, capture_output=True, text=True,
        )
        return [s for s in out.stdout.splitlines() if s]
    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        print(f"WARNING: could not read samples from {vcf_path}: {e}", file=sys.stderr)
        return []


# Some callers (Strelka2's somatic workflow) do not put the real sample IDs
# in the VCF genotype columns at all - they hard-code these generic labels
# instead. vcf2maf.pl still needs to select genotype columns by these exact
# literal names (--vcf-tumor-id/--vcf-normal-id), while --tumor-id/
# --normal-id (the names actually written into the output MAF) must still be
# the real sample IDs, or every Strelka-derived contrast would collide under
# the same "TUMOR"/"NORMAL" barcode once merged cohort-wide.
GENERIC_VCF_SAMPLE_LABELS = {"tumor", "normal"}


def resolve_tumor_normal(contrast, samples, sample_meta):
    """Resolve which VCF sample column is tumor vs normal.

    Returns (tumor_id, normal_id, vcf_tumor_id, vcf_normal_id) where
    tumor_id/normal_id are the real sample IDs to report in the output MAF,
    and vcf_tumor_id/vcf_normal_id are the literal column names present in
    the VCF's genotype columns (identical to tumor_id/normal_id unless the
    caller uses generic labels, e.g. Strelka2's "TUMOR"/"NORMAL").

    Preference order: exact match against contrast tokens, then substring
    match, then samplesheet 'status' (1 = tumor, 0 = normal), then
    positional fallback.
    """
    if "_vs_" in contrast:
        tumor_tok, normal_tok = contrast.split("_vs_", 1)
    else:
        tumor_tok, normal_tok = contrast, None

    if len(samples) == 1:
        return samples[0], None, samples[0], None

    if normal_tok is not None and {s.lower() for s in samples} == GENERIC_VCF_SAMPLE_LABELS:
        vcf_t = next(s for s in samples if s.lower() == "tumor")
        vcf_n = next(s for s in samples if s.lower() == "normal")
        return tumor_tok, normal_tok, vcf_t, vcf_n

    t_id = n_id = None
    for s in samples:
        if s == tumor_tok:
            t_id = s
        if normal_tok is not None and s == normal_tok:
            n_id = s
    if t_id is None or (normal_tok is not None and n_id is None):
        for s in samples:
            if t_id is None and tumor_tok in s:
                t_id = s
            if normal_tok is not None and n_id is None and normal_tok in s:
                n_id = s
    if t_id is None or (normal_tok is not None and n_id is None):
        # fall back to samplesheet status: 1 = tumor, 0 = normal
        for s in samples:
            st = sample_meta.get(s, {}).get("status")
            if st == "1" and t_id is None:
                t_id = s
            elif st == "0" and n_id is None:
                n_id = s
    if t_id is None or (normal_tok is not None and n_id is None):
        # last resort: positional (first = tumor, second = normal)
        remaining = [s for s in samples if s not in (t_id, n_id)]
        if t_id is None and remaining:
            t_id = remaining.pop(0)
        if normal_tok is not None and n_id is None and remaining:
            n_id = remaining.pop(0)
    return t_id, n_id, t_id, n_id


def find_vcfs_for_caller(outdir, caller, override_dir=""):
    """Yield (contrast, subcaller, vcf_path) for a given caller directory.

    If override_dir is set, it is used directly as the caller directory
    (must contain one subdirectory per contrast). Otherwise falls back to
    Sarek's own layout, checking annotation/ first (VEP+SnpEff annotated,
    preferred for downstream vcf2maf) then variant_calling/.
    """
    patterns = CALLER_FILE_PATTERNS[caller]

    def scan(caller_dir):
        for contrast_dir in sorted(caller_dir.iterdir()):
            if not contrast_dir.is_dir():
                continue
            contrast = contrast_dir.name
            for vcf in sorted(contrast_dir.glob("*.vcf.gz")):
                for pattern, subcaller in patterns:
                    if pattern.search(vcf.name):
                        yield contrast, subcaller, vcf
                        break

    if override_dir:
        caller_dir = Path(override_dir)
        if not caller_dir.is_dir():
            print(f"WARNING: --{caller}-dir override '{override_dir}' is not a directory", file=sys.stderr)
            return
        yield from scan(caller_dir)
        return

    for base in ("annotation", "variant_calling"):
        caller_dir = Path(outdir) / base / caller
        if not caller_dir.is_dir():
            continue
        yield from scan(caller_dir)
        # If annotation/ existed and yielded any contrast dirs at all, don't
        # also fall back to variant_calling/ for this caller (avoid dupes).
        if any(True for _ in caller_dir.iterdir()):
            break


CNS_PATTERNS = [
    re.compile(r"\.call\.cns$"),  # cnvkit.py call output (has integer CN), preferred
    re.compile(r"\.cns$"),        # segmented log2-ratio output
]
GENEMETRICS_PATTERNS = [
    re.compile(r"genemetrics.*\.(tsv|txt)$", re.IGNORECASE),
]


def _find_first(directory, patterns):
    if not directory.is_dir():
        return None
    for pattern in patterns:
        for f in sorted(directory.glob("*")):
            if f.is_file() and pattern.search(f.name):
                return f
    return None


def discover_cnvkit(cnvkit_dir, contrast_tumor_pairs):
    """Yield (contrast, tumor_id, cns_path, genemetrics_path) for each
    (contrast, tumor_id) pair already resolved from the VCF-based callers.

    CNVkit is not run by Sarek, so its output layout is whatever the user's
    own CNVkit invocation produced. Checked, in order:
      1. <cnvkit_dir>/<contrast>/*.call.cns (or *.cns), matching the same
         per-contrast subdirectory convention used elsewhere in this pipeline
      2. <cnvkit_dir>/<contrast>.call.cns / <contrast>.cns (flat layout)
      3. same two patterns keyed on the tumor sample id instead of the full
         contrast name (CNVkit is commonly run per-sample, not per-contrast)
    """
    base = Path(cnvkit_dir)
    if not base.is_dir():
        print(f"WARNING: CNVkit directory '{cnvkit_dir}' not found; skipping CNV discovery", file=sys.stderr)
        return

    for contrast, tumor_id in contrast_tumor_pairs:
        cns_path = None
        genemetrics_path = None

        for key in (contrast, tumor_id):
            if not key or cns_path is not None:
                continue
            subdir = base / key
            if subdir.is_dir():
                cns_path = _find_first(subdir, CNS_PATTERNS)
                genemetrics_path = _find_first(subdir, GENEMETRICS_PATTERNS)

        for key in (contrast, tumor_id):
            if not key or cns_path is not None:
                continue
            for suffix in (".call.cns", ".cns"):
                candidate = base / f"{key}{suffix}"
                if candidate.is_file():
                    cns_path = candidate
                    break
            if genemetrics_path is None:
                for f in sorted(base.glob(f"{key}*")):
                    if f.is_file() and any(p.search(f.name) for p in GENEMETRICS_PATTERNS):
                        genemetrics_path = f
                        break

        if cns_path is not None:
            yield contrast, tumor_id, cns_path, genemetrics_path or ""
        else:
            print(f"WARNING: no CNVkit .cns file found for contrast '{contrast}' "
                  f"(tumor_id '{tumor_id}') under {cnvkit_dir}", file=sys.stderr)


def main():
    args = parse_args()
    callers = [c.strip() for c in args.callers.split(",") if c.strip()]
    if args.include_manta and "manta" not in callers:
        callers.append("manta")

    sample_meta = read_samplesheet(args.samplesheet)
    subset_re = re.compile(args.sample_subset_regex) if args.sample_subset_regex else None

    overrides = {
        "mutect2": args.mutect2_dir,
        "strelka": args.strelka_dir,
        "manta": args.manta_dir,
    }

    rows = []
    seen = set()
    for caller in callers:
        if caller not in CALLER_FILE_PATTERNS:
            print(f"WARNING: unknown caller '{caller}', skipping", file=sys.stderr)
            continue
        for contrast, subcaller, vcf_path in find_vcfs_for_caller(args.sarek_outdir, caller, overrides.get(caller, "")):
            if subset_re and not subset_re.search(contrast):
                continue
            key = (contrast, subcaller)
            if key in seen:
                continue
            seen.add(key)

            samples = vcf_samples(vcf_path)
            t_id, n_id, vcf_t_id, vcf_n_id = resolve_tumor_normal(contrast, samples, sample_meta)
            ctype = "tumor_normal" if (n_id and "_vs_" in contrast) else "tumor_only"

            rows.append({
                "contrast": contrast,
                "caller": subcaller,
                "type": ctype,
                "tumor_id": t_id or "",
                "normal_id": n_id or "",
                "vcf_tumor_id": vcf_t_id or "",
                "vcf_normal_id": vcf_n_id or "",
                "vcf_path": str(vcf_path),
                "patient": sample_meta.get(t_id, {}).get("patient", ""),
                "sex": sample_meta.get(t_id, {}).get("sex", ""),
            })

    if not rows:
        print("WARNING: no contrasts discovered - check --sarek-outdir and --callers", file=sys.stderr)

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(args.out, "wt", newline="") as fh:
        fieldnames = ["contrast", "caller", "type", "tumor_id", "normal_id", "vcf_tumor_id", "vcf_normal_id",
                      "vcf_path", "patient", "sex"]
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for r in rows:
            writer.writerow(r)

    print(f"Discovered {len(rows)} (contrast, caller) entries -> {args.out}", file=sys.stderr)

    if args.cnv_out:
        cnvkit_dir = args.cnvkit_dir or str(Path(args.sarek_outdir) / "cnvkit")
        # One (contrast, tumor_id) pair per contrast, regardless of how many
        # snv/indel/sv callers produced it.
        pairs = sorted({(r["contrast"], r["tumor_id"]) for r in rows if r["type"] != ""})
        cnv_rows = []
        for contrast, tumor_id, cns_path, genemetrics_path in discover_cnvkit(cnvkit_dir, pairs):
            cnv_rows.append({
                "contrast": contrast,
                "tumor_id": tumor_id,
                "cns_path": str(cns_path),
                "genemetrics_path": str(genemetrics_path) if genemetrics_path else "",
            })

        Path(args.cnv_out).parent.mkdir(parents=True, exist_ok=True)
        with gzip.open(args.cnv_out, "wt", newline="") as fh:
            fieldnames = ["contrast", "tumor_id", "cns_path", "genemetrics_path"]
            writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            for r in cnv_rows:
                writer.writerow(r)
        print(f"Discovered {len(cnv_rows)} CNVkit contrast entries -> {args.cnv_out}", file=sys.stderr)


if __name__ == "__main__":
    main()
