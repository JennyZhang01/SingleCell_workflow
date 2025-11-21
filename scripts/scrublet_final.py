#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import argparse, os
import scrublet as scr
import scanpy as sc

parser = argparse.ArgumentParser(description="Scrublet final pass on SoupX-corrected counts.")
parser.add_argument("--soupx_dir", required=True, help="soupx_corrected directory with matrix.mtx & tsv(.gz)")
parser.add_argument("--out_dir", required=True)
parser.add_argument("--expected_doublet_rate", type=float, default=0.04)
parser.add_argument("--threshold", type=float, default=0.25)
parser.add_argument("--var_names", default="gene_ids", choices=["gene_symbols","gene_ids"])
args = parser.parse_args()

os.makedirs(args.out_dir, exist_ok=True)

adata = sc.read_10x_mtx(args.soupx_dir, var_names=args.var_names)
adata.var_names_make_unique()

scrub = scr.Scrublet(adata.X, expected_doublet_rate=args.expected_doublet_rate)
scores, calls = scrub.scrub_doublets(n_prin_comps=30)
scrub.call_doublets(threshold=args.threshold)

adata.obs["doublet_score"] = scores
adata.obs["predicted_doublet"] = scrub.predicted_doublets_

adata.write_h5ad(os.path.join(args.out_dir, "soupx_scrublet_all.h5ad"))
adata_filt = adata[~adata.obs["predicted_doublet"], :].copy()
adata_filt.write_h5ad(os.path.join(args.out_dir, "soupx_scrublet_filtered.h5ad"))

fig, ax = scrub.plot_histogram()
fig.savefig(os.path.join(args.out_dir, "scrublet_final_hist.png"), dpi=200, bbox_inches="tight")
print("Saved:", os.path.join(args.out_dir, "soupx_scrublet_filtered.h5ad"))

