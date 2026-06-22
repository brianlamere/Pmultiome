import os
import gzip
import shutil

import scanpy as sc
import pandas as pd
import numpy as np
from scipy import io
from scipy.sparse import csr_matrix

h5_in = "rawdata/LG05/raw_feature_bc_matrix.h5"
outdir = "cb_data/LG05/raw_gex_mex"

os.makedirs(outdir, exist_ok=True)
if os.listdir(outdir):
    raise RuntimeError(f"{outdir} is not empty; delete it first.")

adata = sc.read_10x_h5(h5_in)
gex = adata[:, adata.var["feature_types"] == "Gene Expression"].copy()

# Grab IDs/names
gene_ids = gex.var["gene_ids"].astype(str).to_numpy() if "gene_ids" in gex.var.columns else gex.var_names.astype(str).to_numpy()
gene_names = gex.var_names.astype(str).to_numpy()
feature_types = np.array(["Gene Expression"] * gex.n_vars, dtype=str)

# Ensure sparse
X = gex.X
if not hasattr(X, "tocoo"):
    X = csr_matrix(X)
X = X.tocoo()

# AnnData is typically cells x genes. 10x wants genes x cells.
mtx = X.transpose().tocoo()   # (n_genes, n_barcodes)

# Sanity checks
if mtx.shape != (gex.n_vars, gex.n_obs):
    raise RuntimeError(f"Unexpected shape: wrote {mtx.shape}, expected {(gex.n_vars, gex.n_obs)}")

# Write uncompressed
io.mmwrite(os.path.join(outdir, "matrix.mtx"), mtx)
pd.Series(gex.obs_names.astype(str)).to_csv(os.path.join(outdir, "barcodes.tsv"), sep="\t", header=False, index=False)
pd.DataFrame({0: gene_ids, 1: gene_names, 2: feature_types}).to_csv(
    os.path.join(outdir, "features.tsv"), sep="\t", header=False, index=False
)

# gzip
for fn in ("matrix.mtx", "barcodes.tsv", "features.tsv"):
    p = os.path.join(outdir, fn)
    with open(p, "rb") as f_in, gzip.open(p + ".gz", "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)
    os.remove(p)

print(f"Wrote 10x v3 MEX to {outdir}: features={gex.n_vars}, barcodes={gex.n_obs}")
