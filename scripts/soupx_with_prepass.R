#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(SoupX)
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(DropletUtils)
  library(dplyr)
})

## ---------------------------
## Parse arguments
## ---------------------------

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: soupx_scrublet_soupx.R data_dir prepass_tsv out_dir\n",
       "Example:\n",
       "  Rscript soupx_scrublet_soupx.R \\\n",
       "    data/M1/ \\\n",
       "    qc/M1/prepass/scrublet_prepass.tsv \\\n",
       "    qc/M1/soupx_soupx10x\n")
}
data_dir     <- args[1]
prepass_tsv  <- args[2]
out_dir      <- args[3]


filtered_dir <- paste0(data_dir,"/outs/filtered_feature_bc_matrix/")
raw_dir      <- paste0(data_dir,"/outs/raw_feature_bc_matrix/")

message("prepass_tsv : ", prepass_tsv)
message("data_dir    : ", data_dir)
message("out_dir     : ", out_dir)

## -------------------------------------
## Seruat object and prepass of scrublet 
## -------------------------------------

toc <- Read10X(filtered_dir, gene.column = 1)  # genes x cells
seurat_object <- CreateSeuratObject(counts = toc)

scrub <- fread(prepass_tsv, sep = "\t", header = TRUE)
scrub$barcode <- as.character(scrub$barcode)

meta <- data.frame(
  barcode = colnames(seurat_object),
  row.names = colnames(seurat_object),
  stringsAsFactors = FALSE
)

meta <- merge(
  meta,
  scrub,
  by    = "barcode",
  all.x = TRUE,
  sort  = FALSE
)

rownames(meta) <- meta$barcode

seurat_object <- AddMetaData(
  object   = seurat_object,
  metadata = meta[, c("doublet_score", "call_auto", "call_prepass_mask"), drop = FALSE]
)

## keep singlets only
seurat_singlet <- subset(
  seurat_object,
  subset = is.na(call_prepass_mask) | call_prepass_mask == 0
)

cell_list <- Cells(seurat_singlet)
message("Number of singlets selected: ", length(cell_list))


## ---------------------------
## 2. Build SoupChannel
## ---------------------------

toc <- Read10X(filtered_dir, gene.column = 1)  # genes x cells
tod <- Read10X(raw_dir,      gene.column = 1)

toc_scrublet <- toc[, cell_list]
sc <- SoupChannel(tod, toc_scrublet)

## ---------------------------
## 3. 10X clustering + UMAP
## ---------------------------
TenXCRclusters <- read.csv(
  file = file.path(data_dir, "outs/analysis/clustering/gene_expression_graphclust/clusters.csv"),
  stringsAsFactors = FALSE
)

UMAP <- read.csv(
  file = file.path(data_dir, "outs/analysis/umap/gene_expression_2_components/projection.csv"),
  stringsAsFactors = FALSE
)

sc_meta <- dplyr::inner_join(TenXCRclusters, UMAP, by = "Barcode") %>%
  dplyr::filter(Barcode %in% cell_list)

rownames(sc_meta) <- sc_meta$Barcode

## ---------------------------
## 4. Set clusters
## ---------------------------

common_barcodes <- intersect(colnames(sc$toc), rownames(sc_meta))

if (length(common_barcodes) != ncol(sc$toc)) {
  warning("Not all barcodes in sc$toc have cluster/UMAP info; subsetting to common barcodes.")
  sc$toc      <- sc$toc[, common_barcodes, drop = FALSE]
  sc$metaData <- sc$metaData[common_barcodes, , drop = FALSE]
}

cluster_vec <- setNames(
  sc_meta[colnames(sc$toc), "Cluster"],
  colnames(sc$toc)
)

sc <- setClusters(sc, cluster_vec)

## ---------------------------
## 5. Set UMAP as DR
## ---------------------------


dr_mat <- as.matrix(sc_meta[colnames(sc$toc), c("UMAP.1", "UMAP.2")])
rownames(dr_mat) <- colnames(sc$toc)

sc  <- setDR(sc, dr_mat)


## ---------------------------
## 6. Estimate and adjust counts
## ---------------------------

sc  <- autoEstCont(sc, verbose = FALSE, doPlot = FALSE)
out <- adjustCounts(sc)          # dgCMatrix, genes x cells

## ---------------------------
## 7. Save as 10x-like folder 
## ---------------------------

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

feat_path_v3 <- file.path(filtered_dir, "features.tsv.gz")
feat_path_v2 <- file.path(filtered_dir, "genes.tsv.gz")

if (file.exists(feat_path_v3)) {
  features <- data.table::fread(feat_path_v3, header = FALSE)
} else if (file.exists(feat_path_v2)) {
  features <- data.table::fread(feat_path_v2, header = FALSE)
} else {
  stop("Could not find features.tsv.gz or genes.tsv.gz in ", filtered_dir)
}

# 10x v3: columns = gene_id, gene_name, feature_type
# We assume rownames(out) match features[[1]] (gene_id or gene_name depending on how CellRanger was run)
idx <- match(rownames(out), features[[1]])

if (any(is.na(idx))) {
  stop("Some genes in 'out' not found in features.tsv/genes.tsv; check gene identifiers.")
}

features_sub <- features[idx, ]

## 7b. Write matrix.mtx.gz
mtx_path <- file.path(out_dir, "matrix.mtx.gz")
Matrix::writeMM(out, mtx_path)  # writes plain text to that path

## 7c. Write barcodes.tsv.gz
bc_con <- gzfile(file.path(out_dir, "barcodes.tsv.gz"), "w")
writeLines(colnames(out), bc_con)
close(bc_con)

## 7d. Write features.tsv.gz (or genes.tsv.gz)
feat_out_name <- if (file.exists(feat_path_v3)) "features.tsv.gz" else "genes.tsv.gz"
feat_con <- gzfile(file.path(out_dir, feat_out_name), "w")
write.table(
  features_sub,
  file      = feat_con,
  sep       = "\t",
  row.names = FALSE,
  col.names = FALSE,
  quote     = FALSE
)
close(feat_con)

message("Saved 10x-style SoupX output to: ", out_dir)
message(" - matrix.mtx.gz")
message(" - barcodes.tsv.gz")
message(" - ", feat_out_name)
