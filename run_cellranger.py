#!/usr/bin/env python3
#quick early note: I had already run cellranger prior to writing this,
#didn't want to re-run so this is still thus untested!
"""
run_cellranger.py — stage fastqs and run cellranger-arc, driven entirely by
a manifest.csv (see manifest.csv.template for the exact format).

Usage:
    run_cellranger.py <manifest.csv> [--origin DIR] [--localcores N] [--localmem N]
                       [--dry-run] [--id-prefix STR]

What it does, in order, per unique real_sample_name in the manifest:
    1. Glob ORIGIN/*.fastq.gz, match each file's 10x sample-name field
       against the manifest's fastq_sample_name column.
    2. Symlink matched files into input/<real_sample_name>/<modality>/,
       renamed so the visible name uses real_sample_name instead of the
       raw fastq_sample_name (keeps the eyeball-check you already do).
    3. Write libraries/<real_sample_name>libraries.csv in the format
       cellranger-arc --libraries wants.
    4. Run cellranger-arc count --id=<real_sample_name> ... and on success,
       mv the run dir into cellranger_out/ and symlink
       cra_out/<real_sample_name> -> ../cellranger_out/<real_sample_name>/outs/

Anything that doesn't match cleanly (a fastq with no manifest entry, a
manifest entry with no matching fastq, a sample dir that already exists)
is reported and skipped — never guessed. Fix it by hand and re-run; this
script is idempotent for unfinished/skipped samples.
"""

import argparse
import csv
import os
import re
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

FASTQ_RE = re.compile(r"^(?P<sample>.+)_S\d+_L\d{3}_(?P<read>R[123]|I[12])_001\.fastq\.gz$")


def read_manifest(path: Path):
    with open(path) as f:
        ref_line = f.readline().strip()
        tissue_line = f.readline().strip()
        if not ref_line.lower().startswith("reference:"):
            sys.exit(f"ERROR: expected 'reference: ...' on line 1 of {path}, got: {ref_line!r}")
        if not tissue_line.lower().startswith("tissue:"):
            sys.exit(f"ERROR: expected 'tissue: ...' on line 2 of {path}, got: {tissue_line!r}")
        reference = ref_line.split(":", 1)[1].strip()
        tissue = tissue_line.split(":", 1)[1].strip()
        rows = list(csv.DictReader(f))

    required_cols = {"fastq_sample_name", "real_sample_name", "modality"}
    if not rows or not required_cols.issubset(rows[0].keys()):
        sys.exit(f"ERROR: manifest rows must have columns {required_cols}, got {rows[0].keys() if rows else 'no rows'}")

    return reference, tissue, rows


def group_fastqs(origin_dir: Path, manifest_rows):
    """Match files in origin_dir against manifest fastq_sample_name entries.

    Returns (matches, unmatched_files, unmatched_manifest_rows) where matches
    is {fastq_sample_name: [Path, ...]}.
    """
    wanted = {row["fastq_sample_name"] for row in manifest_rows}
    matches = defaultdict(list)
    unmatched_files = []

    for fq in sorted(origin_dir.glob("*.fastq.gz")):
        m = FASTQ_RE.match(fq.name)
        if not m:
            unmatched_files.append(fq)
            continue
        sample_field = m.group("sample")
        if sample_field in wanted:
            matches[sample_field].append(fq)
        else:
            unmatched_files.append(fq)

    found_sample_fields = set(matches.keys())
    unmatched_manifest_rows = [row for row in manifest_rows if row["fastq_sample_name"] not in found_sample_fields]

    return matches, unmatched_files, unmatched_manifest_rows


def stage_sample(real_sample_name, modality, fastq_sample_name, files, input_root: Path, dry_run: bool):
    target_dir = input_root / real_sample_name / modality
    if target_dir.exists() and any(target_dir.iterdir()):
        print(f"  SKIP (already staged): {target_dir}")
        return target_dir

    print(f"  staging {len(files)} file(s) -> {target_dir}")
    if dry_run:
        for fq in files:
            print(f"    [dry-run] would link {fq.name}")
        return target_dir

    target_dir.mkdir(parents=True, exist_ok=True)
    for fq in files:
        # Replace the raw fastq_sample_name prefix with real_sample_name so the
        # staged filename is human-checkable and never shows the raw ID again.
        new_name = fq.name.replace(fastq_sample_name, f"{real_sample_name}_{modality}", 1)
        link_path = target_dir / new_name
        if link_path.exists():
            print(f"    SKIP (exists): {link_path.name}")
            continue
        rel_target = os.path.relpath(fq.resolve(), target_dir.resolve())
        link_path.symlink_to(rel_target)
        print(f"    linked {link_path.name} -> {rel_target}")

    return target_dir


def write_library_csv(real_sample_name, modality_dirs, libraries_dir: Path, dry_run: bool):
    """modality_dirs: {modality: Path} for this sample, e.g. {'GEX': ..., 'ATAC': ...}"""
    library_type = {"GEX": "Gene Expression", "ATAC": "Chromatin Accessibility"}
    out_path = libraries_dir / f"{real_sample_name}libraries.csv"

    if dry_run:
        print(f"  [dry-run] would write {out_path}")
        return out_path

    libraries_dir.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", newline="\n") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["fastqs", "sample", "library_type"])
        for modality, fq_dir in modality_dirs.items():
            ltype = library_type.get(modality)
            if ltype is None:
                print(f"  WARNING: unknown modality '{modality}' for {real_sample_name}, skipping in libraries.csv")
                continue
            w.writerow([str(fq_dir.resolve()), f"{real_sample_name}_{modality}", ltype])
    print(f"  wrote {out_path}")
    return out_path


def run_cellranger_arc(real_sample_name, library_csv: Path, reference: str,
                        localcores: int, localmem: int, dry_run: bool,
                        cellranger_out_dir: Path, cra_out_dir: Path):
    run_outs = cellranger_out_dir / real_sample_name / "outs"
    cra_link = cra_out_dir / real_sample_name

    if run_outs.exists():
        print(f"  SKIP cellranger-arc (outs already exists): {run_outs}")
    else:
        cmd = [
            "cellranger-arc", "count",
            f"--id={real_sample_name}",
            f"--reference={reference}",
            f"--libraries={library_csv.resolve()}",
            f"--localcores={localcores}",
            f"--localmem={localmem}",
            "--create-bam=false",
        ]
        print(f"  {'[dry-run] would run' if dry_run else 'running'}: {' '.join(cmd)}")
        if not dry_run:
            try:
                result = subprocess.run(cmd)
            except FileNotFoundError:
                print(f"  ERROR: 'cellranger-arc' not found on PATH; skipping {real_sample_name}.")
                return False
            if result.returncode != 0:
                print(f"  ERROR: cellranger-arc failed for {real_sample_name} (exit {result.returncode}); skipping symlink step.")
                return False
            cellranger_out_dir.mkdir(parents=True, exist_ok=True)
            Path(real_sample_name).rename(cellranger_out_dir / real_sample_name)

    if dry_run:
        print(f"  [dry-run] would symlink {cra_link} -> ../cellranger_out/{real_sample_name}/outs/")
        return True

    cra_out_dir.mkdir(parents=True, exist_ok=True)
    if cra_link.exists() or cra_link.is_symlink():
        print(f"  SKIP symlink (already exists): {cra_link}")
    else:
        cra_link.symlink_to(Path("..") / "cellranger_out" / real_sample_name / "outs")
        print(f"  symlinked {cra_link} -> ../cellranger_out/{real_sample_name}/outs/")
    return True


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("manifest", type=Path)
    ap.add_argument("--origin", type=Path, default=Path("originfiles"), help="directory containing raw fastq.gz files")
    ap.add_argument("--input-root", type=Path, default=Path("input"), help="directory to stage symlinked fastqs into")
    ap.add_argument("--libraries-dir", type=Path, default=Path("libraries"))
    ap.add_argument("--cellranger-out-dir", type=Path, default=Path("cellranger_out"))
    ap.add_argument("--cra-out-dir", type=Path, default=Path("cra_out"))
    ap.add_argument("--localcores", type=int, default=16)
    ap.add_argument("--localmem", type=int, default=128)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    reference, tissue, rows = read_manifest(args.manifest)
    print(f"Manifest: reference={reference}  tissue={tissue}  rows={len(rows)}")

    matches, unmatched_files, unmatched_manifest_rows = group_fastqs(args.origin, rows)

    if unmatched_files:
        print(f"\nWARNING: {len(unmatched_files)} fastq file(s) in {args.origin} did not match any manifest entry (or filename pattern); left untouched:")
        for f in unmatched_files:
            print(f"  {f.name}")

    if unmatched_manifest_rows:
        print(f"\nWARNING: {len(unmatched_manifest_rows)} manifest row(s) had no matching fastq files; skipped:")
        for row in unmatched_manifest_rows:
            print(f"  {row['fastq_sample_name']} -> {row['real_sample_name']} ({row['modality']})")

    # Group manifest rows by real_sample_name -> {modality: fastq_sample_name}
    by_sample = defaultdict(dict)
    for row in rows:
        if row["fastq_sample_name"] in matches:
            by_sample[row["real_sample_name"]][row["modality"]] = row["fastq_sample_name"]

    print(f"\n{len(by_sample)} sample(s) ready to process: {sorted(by_sample.keys())}\n")

    for real_sample_name, modality_to_fqsample in sorted(by_sample.items()):
        print(f"=== {real_sample_name} ===")
        try:
            modality_dirs = {}
            for modality, fq_sample_name in sorted(modality_to_fqsample.items()):
                files = matches[fq_sample_name]
                modality_dirs[modality] = stage_sample(
                    real_sample_name, modality, fq_sample_name, files, args.input_root, args.dry_run
                )

            lib_csv = write_library_csv(real_sample_name, modality_dirs, args.libraries_dir, args.dry_run)
            run_cellranger_arc(
                real_sample_name, lib_csv, reference, args.localcores, args.localmem,
                args.dry_run, args.cellranger_out_dir, args.cra_out_dir,
            )
        except Exception as e:
            print(f"  ERROR: unexpected failure on {real_sample_name}: {e}; moving on to next sample.")
        print()

    print("Done. Anything reported as SKIP or WARNING above needs a manual look before trusting results.")


if __name__ == "__main__":
    main()
