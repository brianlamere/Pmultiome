#!/usr/bin/env python3
"""
select_barcodes.py — randomly select a fraction of barcodes from LG38's
raw_feature_bc_matrix.h5, producing a selected_barcodes.csv for use with:

    1. cellranger-arc reanalyze --barcodes=selected_barcodes.csv ...
       (validation: confirm per-cell metrics are preserved at ~FRACTION of
       the original estimated cell count before running Tool 2)

    2. subset_raw_matrix.py
       (transformation: produce a new raw_feature_bc_matrix.h5 containing
       only the selected barcodes, for CellBender)

WHY THIS EXISTS
---------------
LG38 has an unusually high CellRanger-estimated cell count (56,664), which
causes CellBender 0.3.2 to fail with a ZeroDivisionError in
_trim_noiseless_features() regardless of --total-droplets-included or
--low-count-threshold settings above a certain threshold. This is a known
CellBender edge case for overloaded runs. See:
  https://github.com/broadinstitute/CellBender/issues/363

The fix is to reduce the dataset to a representative ~40% random subsample
so that CellRanger's estimated cell count (~22k) allows CellBender to run
with standard 3x total-droplets-included parameters, giving it a proper
empty-droplet training pool. After CellBender, the pipeline continues
identically to all other samples.

REPRODUCTION
------------
This script is fully deterministic given the same input h5 and the same
SEED and FRACTION below. Anyone re-running this script against the same
raw_feature_bc_matrix.h5 will get the identical selected_barcodes.csv.
"""

import csv
import random
import json
from datetime import datetime, timezone
from pathlib import Path

import h5py

# ---------------------------------------------------------------------------
# Configuration — hard-coded for LG38; this is a one-off tool
# ---------------------------------------------------------------------------

RAW_H5       = Path("../cra_out/LG38/raw_feature_bc_matrix.h5")
OUTPUT_CSV   = Path("../project_export/LG38_selected_barcodes.csv")
METADATA_OUT = Path("../project_export/LG38_selected_barcodes_metadata.json")

FRACTION = 0.40   # fraction of all barcodes to select
SEED     = 42     # fixed seed — do not change after first run

# ---------------------------------------------------------------------------

def main():
    if not RAW_H5.exists():
        raise FileNotFoundError(f"Input h5 not found: {RAW_H5}")

    print(f"Reading barcodes from {RAW_H5} ...")
    with h5py.File(RAW_H5, "r") as f:
        raw_barcodes = [b.decode("utf-8") for b in f["matrix/barcodes"][:]]

    total = len(raw_barcodes)
    n_select = round(total * FRACTION)

    print(f"Total barcodes : {total:,}")
    print(f"Fraction       : {FRACTION} ({FRACTION*100:.1f}%)")
    print(f"Seed           : {SEED}")
    print(f"Selecting      : {n_select:,} barcodes")

    rng = random.Random(SEED)
    selected = sorted(rng.sample(raw_barcodes, n_select))

    OUTPUT_CSV.write_text("\n".join(selected) + "\n")
    print(f"Wrote {n_select:,} barcodes to {OUTPUT_CSV}")

    metadata = {
        "tool"           : "select_barcodes.py",
        "run_utc"        : datetime.now(timezone.utc).isoformat(),
        "input_h5"       : str(RAW_H5.resolve()),
        "output_csv"     : str(OUTPUT_CSV.resolve()),
        "seed"           : SEED,
        "fraction"       : FRACTION,
        "total_barcodes" : total,
        "selected_count" : n_select,
        "sample"         : "LG38",
        "reason"         : (
            "CellBender 0.3.2 ZeroDivisionError at high total-droplets-included "
            "values for overloaded run (56,664 estimated cells). Random subsample "
            "to ~40% reduces expected cell count to ~22k, allowing standard 3x "
            "CellBender parameters with a proper empty-droplet training pool."
        ),
    }
    METADATA_OUT.write_text(json.dumps(metadata, indent=2) + "\n")
    print(f"Wrote metadata to {METADATA_OUT}")

    print()
    print("Next step: validate the subsample with cellranger-arc reanalyze.")
    print("Only proceed to subset_raw_matrix.py if per-cell metrics are")
    print("consistent between the original run and the reanalysis.")
    print()
    print("cellranger-arc reanalyze \\")
    print(f"  --id=LG38_subsample_validation \\")
    print(f"  --barcodes={OUTPUT_CSV} \\")
    print(f"  --matrix={RAW_H5.resolve()} \\")
    print(f"  --reference=<path_to_reference> \\")
    print(f"  --atac-fragments=../cra_out/LG38/atac_fragments.tsv.gz")


if __name__ == "__main__":
    main()
