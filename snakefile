shell.executable("/bin/bash")

import os
from snakemake.io import glob_wildcards

# Discover samples by the presence of filtered matrix
data_dir = "/rhome/jzhan413/bigdata/proj/alleleRNA/24-12-17_SingleCell/25-06-10_all_samples/"
raw_samples = glob_wildcards(f"{data_dir}{{sample}}/outs/").sample
SAMPLES = sorted({ s.split("/")[-1] for s in raw_samples if 'c' not in s })   # keep only the last path part

rule all:
    input:
        expand("qc/{sample}/prepass/scrublet_prepass.tsv", sample=SAMPLES),
        expand("qc/{sample}/soupx/matrix.mtx.gz", sample=SAMPLES)
    	"seurat_qc/umap_before_after_harmony.pdf",
	"Processed_data/seurat_obj_harm_with_celltype.rds",
	"seurat_qc/cluster_before_tra.pdf"
	"Processed_data/average_pseudotime.pdf"
    threads: 1
    resources:
        mem_mb = 1000
rule scrublet_prepass:
    input:
        matrix = "data/{sample}/outs/filtered_feature_bc_matrix/matrix.mtx.gz"
    output:
        tsv   = "qc/{sample}/prepass/scrublet_prepass.tsv",
        hist  = "qc/{sample}/prepass/scrublet_prepass_hist.png"
    params:
        filtered_dir = "data/{sample}/outs/filtered_feature_bc_matrix/"
    threads: 1
    resources:
        mem_mb = 60000
    shell:
        r"""
        
        source ~/.bashrc || true
        source "$(conda info --base)/etc/profile.d/conda.sh"  
        
        conda activate /bigdata/wmalab/jzhan413/proj/env/scrublet_pre

        python scripts/scrublet_prepass.py \
          --filtered_dir {params.filtered_dir} \
          --out_dir qc/{wildcards.sample}/prepass \
          --prepass_threshold 0.5
        """

rule soupx_with_prepass:
    input:
        data_dir=lambda wildcards: f"data/{wildcards.sample}/",
        prepass_tsv  = "qc/{sample}/prepass/scrublet_prepass.tsv"
    output:
        mtx    = "qc/{sample}/soupx/matrix.mtx.gz",
        feats  = "qc/{sample}/soupx/features.tsv.gz",
        bcs    = "qc/{sample}/soupx/barcodes.tsv.gz",
    params:
        out_dir="qc/{sample}/soupx"
    threads: 1
    resources:
        mem_mb = 120000
    shell:
        r"""
        
        source ~/.bashrc || true
        source "$(conda info --base)/etc/profile.d/conda.sh"
        conda activate /bigdata/wmalab/jzhan413/proj/env/soupx

        Rscript scripts/soupx_with_prepass.R \
          {input.data_dir} \
          {input.prepass_tsv} \
          {params.out_dir}
        """

rule seurat_cluster:
    input:
    	expand("qc/{sample}/soupx/matrix.mtx.gz", sample=SAMPLES),
	"batch_info.txt"
    output:
	"seurat_qc/umap_before_after_harmony.pdf",
        "Processed_data/seurat_obj_harm_with_celltype.rds"
    threads: 4
    resources:
	mem_mb = 160000
    params:
        project_dir=data_dir.rstrip("/")
    shell:
        r"""
	source ~/.bashrc || true
	source "$(conda info --base)/etc/profile.d/conda.sh"
	conda activate /bigdata/wmalab/jzhan413/proj/env/seurat5

	Rscript scripts/seurat_cluster.R \
	{params.project_dir} \
	20
	"""


rule monocle3_trajectory:
    input:
        "Processed_data/seurat_obj_harm_with_celltype.rds"
    output:
        "seurat_qc/cluster_before_tra.pdf",
	"Processed_data/average_pseudotime.pdf"
    threads: 4
    resources:
        mem_mb = 200000
    shell:
        r"""
        source ~/.bashrc || true
        source "$(conda info --base)/etc/profile.d/conda.sh"
        conda activate /bigdata/wmalab/jzhan413/proj/env/monocle3

        Rscript scripts/monocles.R \
        {data_dir.rstrip("/")}
        """
