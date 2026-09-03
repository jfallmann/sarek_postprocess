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

Usage:
    discover_contrasts.py --sarek-outdir DIR --samplesheet CSV \
        --callers mutect2,strelka --out contrasts.tsv \
        [--sample-subset-regex REGEX]
"""
import argparse
import csv
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


def resolve_tumor_normal(contrast, samples, sample_meta):
    """Resolve which VCF sample column is tumor vs normal.

    Preference order: exact match against contrast tokens, then substring
    match, then samplesheet 'status' (1 = tumor, 0 = normal), then
    positional fallback.
    """
    if "_vs_" in contrast:
        tumor_tok, normal_tok = contrast.split("_vs_", 1)
    else:
        tumor_tok, normal_tok = contrast, None

    if len(samples) == 1:
        return samples[0], None

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
    return t_id, n_id


def find_vcfs_for_caller(outdir, caller):
    """Yield (contrast, subcaller, vcf_path) for a given caller directory,
    checking annotation/ first (VEP+SnpEff annotated, preferred for
    downstream vcf2maf) and falling back to variant_calling/ if absent."""
    patterns = CALLER_FILE_PATTERNS[caller]
    for base in ("annotation", "variant_calling"):
        caller_dir = Path(outdir) / base / caller
        if not caller_dir.is_dir():
            continue
        for contrast_dir in sorted(caller_dir.iterdir()):
            if not contrast_dir.is_dir():
                continue
            contrast = contrast_dir.name
            for vcf in sorted(contrast_dir.glob("*.vcf.gz")):
                for pattern, subcaller in patterns:
                    if pattern.search(vcf.name):
                        yield contrast, subcaller, vcf
                        break
        # If annotation/ existed and yielded any contrast dirs at all, don't
        # also fall back to variant_calling/ for this caller (avoid dupes).
        if any(True for _ in caller_dir.iterdir()):
            break


def main():
    args = parse_args()
    callers = [c.strip() for c in args.callers.split(",") if c.strip()]
    if args.include_manta and "manta" not in callers:
        callers.append("manta")

    sample_meta = read_samplesheet(args.samplesheet)
    subset_re = re.compile(args.sample_subset_regex) if args.sample_subset_regex else None

    rows = []
    seen = set()
    for caller in callers:
        if caller not in CALLER_FILE_PATTERNS:
            print(f"WARNING: unknown caller '{caller}', skipping", file=sys.stderr)
            continue
        for contrast, subcaller, vcf_path in find_vcfs_for_caller(args.sarek_outdir, caller):
            if subset_re and not subset_re.search(contrast):
                continue
            key = (contrast, subcaller)
            if key in seen:
                continue
            seen.add(key)

            samples = vcf_samples(vcf_path)
            t_id, n_id = resolve_tumor_normal(contrast, samples, sample_meta)
            ctype = "tumor_normal" if (n_id and "_vs_" in contrast) else "tumor_only"

            rows.append({
                "contrast": contrast,
                "caller": subcaller,
                "type": ctype,
                "tumor_id": t_id or "",
                "normal_id": n_id or "",
                "vcf_path": str(vcf_path),
                "patient": sample_meta.get(t_id, {}).get("patient", ""),
                "sex": sample_meta.get(t_id, {}).get("sex", ""),
            })

    if not rows:
        print("WARNING: no contrasts discovered - check --sarek-outdir and --callers", file=sys.stderr)

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w", newline="") as fh:
        fieldnames = ["contrast", "caller", "type", "tumor_id", "normal_id", "vcf_path", "patient", "sex"]
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for r in rows:
            writer.writerow(r)

    print(f"Discovered {len(rows)} (contrast, caller) entries -> {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
