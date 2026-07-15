import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

pbm = pd.read_csv("../cellranger_out/LG38_HXB2/outs/per_barcode_metrics.csv")
selected = set(Path("../project_export/LG38_selected_barcodes.csv").read_text().splitlines())

orig_cells   = pbm[pbm["is_cell"] == 1]
subset_cells = orig_cells[orig_cells["barcode"].isin(selected)]
rest_cells   = orig_cells[~orig_cells["barcode"].isin(selected)]

fig, axes = plt.subplots(2, 2, figsize=(12, 8))
fig.suptitle("LG38 subsample validation: distribution shape comparison\n"
             f"Original: {len(orig_cells):,} cells  |  "
             f"Subset (40%): {len(subset_cells):,} cells  |  "
             f"Excluded (60%): {len(rest_cells):,} cells")

metrics = [
    ("gex_umis_count",       "GEX UMI counts per cell",       50),
    ("atac_fragments",       "ATAC fragments per cell",        50),
    ("gex_genes_count",      "GEX genes detected per cell",    50),
    ("atac_TSS_fragments",   "ATAC TSS fragments per cell",    50),
]

for ax, (col, title, bins) in zip(axes.flat, metrics):
    if col not in pbm.columns:
        ax.text(0.5, 0.5, f"{col}\nnot found", ha="center", va="center",
                transform=ax.transAxes)
        continue

    # Use log scale for count metrics, linear for enrichment scores
    use_log = col in ("gex_umis_count", "atac_fragments", "gex_genes_detected")
    data_orig   = orig_cells[col].dropna()
    data_subset = subset_cells[col].dropna()

    if use_log:
        data_orig   = np.log10(data_orig.clip(lower=1))
        data_subset = np.log10(data_subset.clip(lower=1))
        ax.set_xlabel(f"log10({col})")
    else:
        ax.set_xlabel(col)

    bin_edges = np.linspace(
        min(data_orig.min(), data_subset.min()),
        max(data_orig.quantile(0.999), data_subset.quantile(0.999)),
        bins + 1
    )

    ax.hist(data_orig,   bins=bin_edges, alpha=0.5, density=True,
            color="steelblue", label=f"Original (n={len(data_orig):,})")
    ax.hist(data_subset, bins=bin_edges, alpha=0.5, density=True,
            color="salmon",    label=f"Subset (n={len(data_subset):,})")
    ax.set_title(title)
    ax.set_ylabel("Density")
    ax.legend(fontsize=8)

plt.tight_layout()
plt.savefig("../project_export/LG38_subsample_validation.png", dpi=150,
            bbox_inches="tight")
print("Saved to ../project_export/LG38_subsample_validation.png")
