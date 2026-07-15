#!/usr/bin/env python3
"""
run_cellbender.py — run CellBender remove-background per sample, driven by
manifest.csv.

For each unique real_sample_name in the manifest:
    1. Read cra_out/<sample>/summary.csv for CellRanger's own "Estimated
       number of cells", and use it as --expected-cells (with
       --total-droplets-included set to 3x that). CellBender's own
       self-inference under-calls on this lab's overloaded runs.

       Per-sample overrides for expected_cells and/or total_droplets_included
       may be set in manifest.csv (add those columns; leave blank for the
       standard behaviour). Overrides are recorded in the run log so the
       deviation is visible without reading the manifest separately.

    2. Extract Gene Expression features from
       cra_out/<sample>/raw_feature_bc_matrix.h5 into
       cb_data/<sample>/raw_gex.h5 (CellRanger v3 h5 format). CellBender
       should never see the combined GEX+ATAC matrix directly.
    3. Run `cellbender remove-background` on the GEX-only matrix.
    4. Move the resulting *_filtered.h5 to
       cb_data/<sample>/cellbender_gex_filtered.h5 (the path
       system_settings.R expects).

All samples run with the same code path. Per-sample parameter differences
live in the manifest, not in script logic.
"""

import csv
import sys
from pathlib import Path

CRA_OUT_DIR = Path("../cra_out")
CB_DATA_DIR = Path("../cb_data")
RAW_H5_NAME = "raw_feature_bc_matrix.h5"
GEX_H5_NAME = "raw_gex.h5"
FINAL_H5_NAME = "cellbender_gex_filtered.h5"  # must match cellbender_rna_h5filename in system_settings.R


def get_unique_samples(manifest_path: Path):
    with open(manifest_path) as f:
        f.readline()  # reference: ...
        f.readline()  # Tissue: ...
        rows = list(csv.DictReader(f))
    return sorted({row["real_sample_name"] for row in rows})


def read_manifest_overrides(manifest_path: Path) -> dict:
    """Return per-sample CellBender parameter overrides from manifest.csv.

    Recognises two optional columns: expected_cells and
    total_droplets_included. Either or both may be present; either or both
    may be blank for any given sample (blank == use standard behaviour).

    Because each sample appears on multiple rows (one per modality), the
    first non-empty value found for each column wins. Conflicting values
    across rows for the same sample are flagged as a warning -- the manifest
    should be the single source of truth, so ambiguity there is worth
    surfacing rather than silently resolving.

    Returns a dict keyed by real_sample_name whose values are dicts with
    only the keys that actually have overrides set, e.g.:
        {
            "LG38": {"expected_cells": 56664, "total_droplets_included": 100000},
        }
    Samples with no overrides do not appear in the dict at all, so a simple
    `overrides.get(sample, {})` in the caller is always safe.
    """
    override_columns = {"expected_cells", "total_droplets_included"}

    with open(manifest_path) as f:
        f.readline()  # reference: ...
        f.readline()  # Tissue: ...
        rows = list(csv.DictReader(f))

    # Which override columns are actually present in this manifest?
    if not rows:
        return {}
    present = override_columns & rows[0].keys()
    if not present:
        return {}

    overrides = {}          # real_sample_name -> {col: int}
    seen_values = {}        # real_sample_name -> {col: value_from_first_row}

    for row in rows:
        sample = row["real_sample_name"]
        for col in present:
            raw = row.get(col, "").strip()
            if not raw:
                continue
            try:
                value = int(raw)
            except ValueError:
                sys.exit(
                    f"ERROR: manifest column '{col}' for sample '{sample}' "
                    f"must be an integer or blank, got: {raw!r}"
                )
            if sample not in seen_values:
                seen_values[sample] = {}
            if col in seen_values[sample]:
                if seen_values[sample][col] != value:
                    print(
                        f"WARNING: manifest has conflicting values for "
                        f"'{col}' on sample '{sample}': "
                        f"{seen_values[sample][col]} vs {value}. "
                        f"Using the first value seen ({seen_values[sample][col]})."
                    )
            else:
                seen_values[sample][col] = value
                overrides.setdefault(sample, {})[col] = value

    return overrides


def read_cellranger_estimate(cra_sample_dir: Path):
    summary_csv = cra_sample_dir / "summary.csv"
    if not summary_csv.exists():
        sys.exit(f"FATAL: {summary_csv} not found. Has cellranger-arc been run for this sample?")
    with open(summary_csv) as f:
        row = next(csv.DictReader(f))
    if "Estimated number of cells" not in row:
        sys.exit(f"FATAL: {summary_csv} has no 'Estimated number of cells' column. Format may have changed.")
    return int(row["Estimated number of cells"])


def extract_gex(raw_h5_path: Path, out_h5_path: Path):
    """Extract Gene Expression features only, write as CellRanger v3 h5."""
    if out_h5_path.exists():
        print(f"  SKIP GEX extraction (already exists): {out_h5_path}")
        return

    import scanpy as sc
    import h5py
    import numpy as np
    from scipy.sparse import csc_matrix

    print(f"  reading {raw_h5_path}...")
    adata = sc.read_10x_h5(str(raw_h5_path), gex_only=False)
    gex = adata[:, adata.var["feature_types"] == "Gene Expression"].copy()
    print(f"  extracted {gex.n_vars} Gene Expression features x {gex.n_obs} barcodes "
          f"(out of {adata.n_vars} total features read)")

    dup_mask = gex.var_names.duplicated(keep=False)
    if dup_mask.any():
        dup_names = sorted(set(gex.var_names[dup_mask]))
        print(f"  WARNING: {len(dup_names)} duplicate gene symbol(s): {dup_names[:20]}"
              f"{' ...' if len(dup_names) > 20 else ''}")

    out_h5_path.parent.mkdir(parents=True, exist_ok=True)

    X = gex.X
    X = (X.tocsc() if hasattr(X, "tocsc") else csc_matrix(X))
    X = X.T.tocsc()

    gene_ids = (gex.var["gene_ids"].astype(str).to_numpy()
                if "gene_ids" in gex.var.columns else gex.var_names.astype(str).to_numpy())
    gene_names = gex.var_names.astype(str).to_numpy()
    n_genes = gex.n_vars
    genome = (gex.var["genome"].astype(str).to_numpy()
              if "genome" in gex.var.columns else np.array([""] * n_genes, dtype=str))
    barcodes = gex.obs_names.astype(str).to_numpy()

    with h5py.File(out_h5_path, "w") as f:
        m = f.create_group("matrix")
        m.create_dataset("barcodes", data=barcodes.astype("S"))
        m.create_dataset("data", data=X.data.astype("uint32"))
        m.create_dataset("indices", data=X.indices.astype("uint32"))
        m.create_dataset("indptr", data=X.indptr.astype("uint32"))
        m.create_dataset("shape", data=np.array(X.shape, dtype="uint64"))
        feat = m.create_group("features")
        feat.create_dataset("id", data=gene_ids.astype("S"))
        feat.create_dataset("name", data=gene_names.astype("S"))
        feat.create_dataset("feature_type", data=np.array(["Gene Expression"] * n_genes, dtype="S"))
        feat.create_dataset("genome", data=genome.astype("S"))
        feat.create_dataset("_all_tag_keys", data=np.array([], dtype="S20"))

    print(f"  wrote {out_h5_path}")


def run_cellbender(input_h5: Path, output_h5: Path, expected_cells: int, total_droplets: int):
    import subprocess
    if output_h5.exists():
        print(f"  SKIP cellbender (output already exists): {output_h5}")
        return
    output_h5.parent.mkdir(parents=True, exist_ok=True)
    cmd = ["cellbender", "remove-background",
           "--input", str(input_h5.resolve()),
           "--output", str(output_h5.resolve()),
           "--expected-cells", str(expected_cells),
           "--total-droplets-included", str(total_droplets),
           "--cuda"]
    print(f"  running: {' '.join(cmd)}")
    try:
        result = subprocess.run(cmd, cwd=str(output_h5.parent))
    except FileNotFoundError:
        sys.exit("FATAL: 'cellbender' not found on PATH.")
    if result.returncode != 0:
        sys.exit(f"FATAL: cellbender failed (exit {result.returncode}) on {input_h5}.")


def main():
    if len(sys.argv) != 2:
        sys.exit(f"usage: {sys.argv[0]} manifest.csv")
    manifest_path = Path(sys.argv[1])

    samples = get_unique_samples(manifest_path)
    overrides = read_manifest_overrides(manifest_path)

    print(f"{len(samples)} unique sample(s) in manifest: {samples}")
    if overrides:
        print(f"Per-sample manifest overrides found for: {sorted(overrides.keys())}")
    print()

    for sample in samples:
        print(f"=== {sample} ===")
        cra_sample_dir = CRA_OUT_DIR / sample
        raw_h5 = cra_sample_dir / RAW_H5_NAME

        if not cra_sample_dir.exists():
            sys.exit(f"FATAL [{sample}]: {cra_sample_dir} does not exist. Run run_cellranger.py first.")
        if not raw_h5.exists():
            sys.exit(f"FATAL [{sample}]: {raw_h5} not found.")

        sample_overrides = overrides.get(sample, {})

        if "expected_cells" in sample_overrides:
            expected_cells = sample_overrides["expected_cells"]
            print(f"  expected_cells: {expected_cells} (manifest override)")
        else:
            expected_cells = read_cellranger_estimate(cra_sample_dir)
            print(f"  expected_cells: {expected_cells} (CellRanger estimate)")

        if "total_droplets_included" in sample_overrides:
            total_droplets = sample_overrides["total_droplets_included"]
            print(f"  total_droplets_included: {total_droplets} (manifest override)")
        else:
            total_droplets = expected_cells * 3
            print(f"  total_droplets_included: {total_droplets} (3x expected_cells)")

        sample_cb_dir = CB_DATA_DIR / sample
        gex_h5 = sample_cb_dir / GEX_H5_NAME
        extract_gex(raw_h5, gex_h5)

        cb_output = sample_cb_dir / "cellbender_gex.h5"
        run_cellbender(gex_h5, cb_output, expected_cells, total_droplets)

        filtered_output = sample_cb_dir / "cellbender_gex_filtered.h5"
        final_path = sample_cb_dir / FINAL_H5_NAME
        if filtered_output.exists() and filtered_output.resolve() != final_path.resolve():
            filtered_output.rename(final_path)
        print(f"  OK -> {final_path}")
        print()

    print("Done.")


if __name__ == "__main__":
    main()
