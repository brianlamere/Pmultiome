#!/usr/bin/env python3
"""
subset_raw_matrix.py — subset LG38's raw_feature_bc_matrix.h5 to the
barcodes in selected_barcodes.csv, renaming the original file first so it
is preserved but out of the way. The new file takes the original name so
that run_cellbender.py (and everything downstream) sees no difference.

RUN THIS ONLY AFTER cellranger-arc reanalyze has validated that the
selected barcodes are a representative subsample. This script is a
one-way transformation of the input data; it does not touch the original
CellRanger run directory, only the cb_data staging area used by
run_cellbender.py. The original raw_feature_bc_matrix.h5 is renamed to
orig_raw_feature_bc_matrix.h5 in the same directory, not deleted.

See select_barcodes.py for the reasoning behind this subsampling.
"""

import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import h5py
import numpy as np
from scipy.sparse import csc_matrix

# ---------------------------------------------------------------------------
# Configuration — hard-coded for LG38; this is a one-off tool
# ---------------------------------------------------------------------------

# The raw GEX h5 that run_cellbender.py extracted from CellRanger output.
# This is the file in cb_data/, not the original CellRanger output directory.
RAW_H5       = Path("../cb_data/LG38/raw_gex.h5")

BARCODES_CSV = Path("../project_export/LG38_selected_barcodes.csv")   # from select_barcodes.py
METADATA_OUT = Path("../project_export/LG38_subset_raw_matrix_metadata.json")

# ---------------------------------------------------------------------------

def load_selected_barcodes(csv_path: Path) -> set:
    lines = csv_path.read_text().splitlines()
    # Strip optional header if it doesn't look like a barcode sequence
    # (barcodes are 16 uppercase ACGT chars followed by -1)
    barcodes = []
    for line in lines:
        line = line.strip()
        if line and not line.lower().startswith("barcode"):
            barcodes.append(line)
    return set(barcodes)


def main():
    if not RAW_H5.exists():
        sys.exit(f"FATAL: input h5 not found: {RAW_H5}\n"
                 f"Has run_cellbender.py been run to extract the GEX h5 from "
                 f"the CellRanger output?")

    if not BARCODES_CSV.exists():
        sys.exit(f"FATAL: barcode list not found: {BARCODES_CSV}\n"
                 f"Run select_barcodes.py first, then validate with "
                 f"cellranger-arc reanalyze before running this script.")

    selected = load_selected_barcodes(BARCODES_CSV)
    print(f"Loaded {len(selected):,} selected barcodes from {BARCODES_CSV}")

    # Rename original before touching anything else — if anything below
    # fails, the original is still intact under its new name.
    orig_path = RAW_H5.parent / f"orig_{RAW_H5.name}"
    if orig_path.exists():
        sys.exit(
            f"FATAL: {orig_path} already exists.\n"
            f"This script has already been run, or a previous run was "
            f"interrupted. If you want to re-run, restore the original:\n"
            f"  mv {orig_path} {RAW_H5}"
        )

    print(f"Renaming {RAW_H5} -> {orig_path} ...")
    RAW_H5.rename(orig_path)

    print(f"Reading original h5 ...")
    with h5py.File(orig_path, "r") as f:
        all_barcodes  = np.array([b.decode("utf-8") for b in f["matrix/barcodes"][:]])
        gene_ids      = f["matrix/features/id"][:]
        gene_names    = f["matrix/features/name"][:]
        feature_types = f["matrix/features/feature_type"][:]
        genome        = f["matrix/features/genome"][:]
        data          = f["matrix/data"][:]
        indices       = f["matrix/indices"][:]
        indptr        = f["matrix/indptr"][:]
        shape         = f["matrix/shape"][:]  # [n_genes, n_barcodes]

    n_genes, n_barcodes_orig = int(shape[0]), int(shape[1])
    print(f"Original matrix: {n_genes:,} genes x {n_barcodes_orig:,} barcodes")

    # Build boolean mask over all barcodes in the original h5
    keep_mask = np.array([bc in selected for bc in all_barcodes])
    n_kept = keep_mask.sum()

    if n_kept == 0:
        # Restore original before dying so the user isn't left with nothing
        orig_path.rename(RAW_H5)
        sys.exit(
            "FATAL: 0 barcodes matched between the h5 and selected_barcodes.csv.\n"
            "Check that the barcode format (e.g. ACGTACGT...-1) matches between "
            "the two files. Original h5 has been restored."
        )

    not_found = len(selected) - n_kept
    if not_found > 0:
        print(f"WARNING: {not_found:,} barcodes in {BARCODES_CSV} were not "
              f"found in the h5 (may be in whitelist but have zero counts -- "
              f"this is normal and expected for a raw matrix).")

    print(f"Keeping {n_kept:,} of {n_barcodes_orig:,} barcodes "
          f"({n_kept/n_barcodes_orig*100:.1f}%)")

    # Reconstruct sparse matrix and subset columns (barcodes are columns)
    # h5 stores genes x barcodes in CSC layout (indptr over barcodes)
    mat = csc_matrix((data, indices, indptr), shape=(n_genes, n_barcodes_orig))
    mat_subset = mat[:, keep_mask].tocsc()

    kept_barcodes = all_barcodes[keep_mask]
    print(f"Writing subsetted h5 to {RAW_H5} ...")

    with h5py.File(RAW_H5, "w") as f:
        m = f.create_group("matrix")
        m.create_dataset("barcodes",
                         data=np.array([b.encode("utf-8") for b in kept_barcodes],
                                       dtype="S"))
        m.create_dataset("data",    data=mat_subset.data.astype("uint32"))
        m.create_dataset("indices", data=mat_subset.indices.astype("uint32"))
        m.create_dataset("indptr",  data=mat_subset.indptr.astype("uint32"))
        m.create_dataset("shape",
                         data=np.array([n_genes, n_kept], dtype="uint64"))
        feat = m.create_group("features")
        feat.create_dataset("id",           data=gene_ids)
        feat.create_dataset("name",         data=gene_names)
        feat.create_dataset("feature_type", data=feature_types)
        feat.create_dataset("genome",       data=genome)
        feat.create_dataset("_all_tag_keys",
                            data=np.array([], dtype="S20"))

    print(f"Done. Subsetted matrix: {n_genes:,} genes x {n_kept:,} barcodes")

    metadata = {
        "tool"              : "subset_raw_matrix.py",
        "run_utc"           : datetime.now(timezone.utc).isoformat(),
        "input_h5_original" : str(orig_path.resolve()),
        "input_barcodes"    : str(BARCODES_CSV.resolve()),
        "output_h5"         : str(RAW_H5.resolve()),
        "n_genes"           : int(n_genes),
        "n_barcodes_orig"   : int(n_barcodes_orig),
        "n_barcodes_kept"   : int(n_kept),
        "n_barcodes_csv"    : int(len(selected)),
        "n_barcodes_not_in_h5": int(not_found),
        "fraction_kept"     : round(float(n_kept) / float(n_barcodes_orig), 6),
        "sample"            : "LG38",
        "note"              : (
            "Original h5 preserved as orig_raw_gex.h5 in the same directory. "
            "Subsetted h5 takes the original filename so run_cellbender.py "
            "requires no changes. This transformation was applied after "
            "cellranger-arc reanalyze confirmed the selected barcodes are a "
            "representative subsample of the full dataset."
        ),
    }
    METADATA_OUT.write_text(json.dumps(metadata, indent=2) + "\n")
    print(f"Wrote metadata to {METADATA_OUT}")
    print()
    print("run_cellbender.py can now be run normally for LG38.")
    print(f"Original h5 preserved at: {orig_path}")


if __name__ == "__main__":
    main()
