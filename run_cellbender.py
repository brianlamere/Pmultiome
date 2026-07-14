#!/usr/bin/env python3
"""
run_cellbender.py — run CellBender remove-background per sample, driven by
manifest.csv.

For each unique real_sample_name in the manifest:
    1. Read cra_out/<sample>/summary.csv for CellRanger's own "Estimated
       number of cells", and use it as --expected-cells (with
       --total-droplets-included set to 3x that). CellBender's own
       self-inference under-calls on this lab's overloaded runs.
    2. Extract Gene Expression features from
       cra_out/<sample>/raw_feature_bc_matrix.h5 into
       cb_data/<sample>/raw_gex.h5 (CellRanger v3 h5 format). CellBender
       should never see the combined GEX+ATAC matrix directly.
    3. Run `cellbender remove-background` on the GEX-only matrix.
    4. Move the resulting *_filtered.h5 to
       cb_data/<sample>/cellbender_gex_filtered.h5 (the path
       system_settings.R expects).

If a sample needs different handling, run cellbender by hand for that one
sample. This script does the common case only.
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
    # gex_only=False: scanpy's default silently drops Peaks before our
    # own feature_types filter below would ever run.
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
    X = X.T.tocsc()  # CellRanger h5 stores genes x barcodes

    gene_ids = (gex.var["gene_ids"].astype(str).to_numpy()
                if "gene_ids" in gex.var.columns else gex.var_names.astype(str).to_numpy())
    gene_names = gex.var_names.astype(str).to_numpy()
    n_genes = gex.n_vars
    genome = (gex.var["genome"].astype(str).to_numpy()
              if "genome" in gex.var.columns else np.array([""] * n_genes, dtype=str))
    barcodes = gex.obs_names.astype(str).to_numpy()

    # Fixed-width byte strings, not h5py's vlen dtype: PyTables (which
    # CellBender's loader uses) can't read vlen strings.
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
    print(f"{len(samples)} unique sample(s) in manifest: {samples}\n")

    for sample in samples:
        print(f"=== {sample} ===")
        cra_sample_dir = CRA_OUT_DIR / sample
        raw_h5 = cra_sample_dir / RAW_H5_NAME

        if not cra_sample_dir.exists():
            sys.exit(f"FATAL [{sample}]: {cra_sample_dir} does not exist. Run run_cellranger.py first.")
        if not raw_h5.exists():
            sys.exit(f"FATAL [{sample}]: {raw_h5} not found.")

        expected_cells_estimate = overrides.get(sample, {}).get("expected_cells") or read_cellranger_estimate(cra_sample_dir)
        total_droplets = overrides.get(sample, {}).get("total_droplets_included") or expected_cells_estimate * 3
        source = "manifest override" if sample in overrides else "CellRanger estimate"
        print(f"  Estimated cells: {expected_cells_estimate} ({source}, "
              f"--total-droplets-included {total_droplets})")

        sample_cb_dir = CB_DATA_DIR / sample
        gex_h5 = sample_cb_dir / GEX_H5_NAME
        extract_gex(raw_h5, gex_h5)

        cb_output = sample_cb_dir / "cellbender_gex.h5"
        run_cellbender(gex_h5, cb_output, expected_cells_estimate, total_droplets)

        filtered_output = sample_cb_dir / "cellbender_gex_filtered.h5"
        final_path = sample_cb_dir / FINAL_H5_NAME
        if filtered_output.exists() and filtered_output.resolve() != final_path.resolve():
            filtered_output.rename(final_path)
        print(f"  OK -> {final_path}")
        print()

    print("Done.")


if __name__ == "__main__":
    main()
