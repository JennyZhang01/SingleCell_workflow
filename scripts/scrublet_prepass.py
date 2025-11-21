#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#!/usr/bin/env python3
import argparse, os, numpy as np
import scrublet as scr
import scanpy as sc
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description="Scrublet pre-pass (mask only).")
parser.add_argument("--filtered_dir", required=True, help="Cell Ranger filtered_feature_bc_matrix directory")
parser.add_argument("--out_dir", required=True)
parser.add_argument("--expected_doublet_rate", type=float, default=0.04)
parser.add_argument("--prepass_threshold", type=float, default=0.50)
parser.add_argument("--var_names", default="gene_ids", choices=["gene_symbols","gene_ids"])
args = parser.parse_args()

os.makedirs(args.out_dir, exist_ok=True)

adata = sc.read_10x_mtx(args.filtered_dir, var_names=args.var_names)
adata.var_names_make_unique()

scrub = scr.Scrublet(adata.X, expected_doublet_rate=args.expected_doublet_rate)
scores, calls_auto = scrub.scrub_doublets(n_prin_comps=30)
mask = (scores > args.prepass_threshold).astype(int)

# Save WITH barcodes for safe merging
out_path = os.path.join(args.out_dir, "scrublet_prepass.tsv")
np.savetxt(
    out_path,
    np.c_[np.array(adata.obs_names, dtype=str), scores, calls_auto.astype(int), mask],
    fmt="%s",
    delimiter="\t",
    header="barcode\tdoublet_score\tcall_auto\tcall_prepass_mask",
    comments=""
)

# Histogram
fig, ax = scrub.plot_histogram()
fig.savefig(os.path.join(args.out_dir, "scrublet_prepass_hist.png"), dpi=200, bbox_inches="tight")
print("Wrote:", out_path)


