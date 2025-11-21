suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(dplyr)
  library(harmony)
  library(tibble)
  library(sctransform)
  library(ggplot2)
  library(gridExtra)
  library(RColorBrewer)
  library(patchwork) 
  library(future)
})
options(future.globals.maxSize = 10 * 1024^3)

setwd("/rhome/jzhan413/bigdata/proj/alleleRNA/24-12-17_SingleCell/25-06-10_all_samples/")

############### reading batch information ############## 
batch_table <- fread("batch_info.txt")
batch_map <- setNames(batch_table$batch, batch_table$sample)



mt_threshold <- 20   # change to 25 or whatever you decide
out_dir <- "qc/"
out_seurat_dir <- "Processed_data"
if (!dir.exists(out_seurat_dir)){
  dir.create(out_seurat_dir,recursive = TRUE)
} 


soupx_dirs <- list.dirs("qc", full.names = TRUE, recursive = TRUE)
soupx_dirs <- soupx_dirs[grepl("soupx$", soupx_dirs)]

sample_names <- basename(dirname(soupx_dirs))

ord <- order(sample_names)
soupx_dirs   <- soupx_dirs[ord]
sample_names <- sample_names[ord]

message("Detected SoupX dirs:")
for (i in seq_along(soupx_dirs)) {
  message("  sample = ", sample_names[i], "  dir = ", soupx_dirs[i])
}
batch_map <- setNames(batch_table$batch, batch_table$sample)






## ---------------------------
## 2. Read each SoupX-corrected matrix as a Seurat object
## ---------------------------

seurat_list <- vector("list", length(soupx_dirs))
names(seurat_list) <- sample_names

for (i in seq_along(soupx_dirs)) {
  dir_i <- soupx_dirs[i]
  samp  <- sample_names[i]
  batch_label <- unname(batch_map[samp])
  
  message("[", i, "/", length(soupx_dirs), "] Reading sample '", samp,
          "' from: ", dir_i)
  
  counts <- Read10X(data.dir = dir_i)  # genes x cells
  obj <- CreateSeuratObject(counts = counts, project = samp)
  obj$sample <- samp                  
  obj$batch <- batch_label
  seurat_list[[samp]] <- obj
}

## ---------------------------
## 3. Merge all samples into one Seurat object
## ---------------------------

seurat_obj <- seurat_list[[1]]
if (length(seurat_list) > 1) {
  for (k in 2:length(seurat_list)) {
    seurat_obj <- merge(seurat_obj, y = seurat_list[[k]])
  }
}

message("Total merged cells (before any QC filtering): ", ncol(seurat_obj))


####### adding log1p_total_count in seurat object###### 
seurat_obj$log1p_total_count <- log1p(seurat_obj$nCount_RNA)
 
## ---------------------------
## 4. Compute QC metrics
## ---------------------------

# Add mitochondrial percentage (human MT genes)
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

qc_raw <- seurat_obj@meta.data %>%
  rownames_to_column("barcode")

write.csv(
  qc_raw,
  file = file.path(out_dir, "qc_metrics_raw.csv"),
  row.names = FALSE
)

message("Wrote QC metrics table: ", file.path(out_dir, "qc_metrics_raw.csv"))

## ---------------------------
## 5. QC violin plots grouped by sample
## ---------------------------
Idents(seurat_obj) <- seurat_obj$sample

qc_features <- c("nFeature_RNA", "nCount_RNA", "percent.mt")

qc_vln <- VlnPlot(
  seurat_obj,
  features = qc_features,
  pt.size = 0.1,
  ncol = 3
) +
  ggtitle("QC metrics per sample") +
  theme_bw()

pdf(file.path(out_dir, "qc_violin_by_sample.pdf"), width = 12, height = 4)
print(qc_vln)
dev.off()

## ---------------------------
## 6. QC distributions of percent.mt (mitochondrial RNA)
## ---------------------------

qc_df <- seurat_obj@meta.data %>% tibble::rownames_to_column("barcode")

mt_hist <- ggplot(qc_df, aes(x = percent.mt, fill = sample)) +
  geom_histogram(bins = 50, alpha = 0.6) +
  facet_wrap(~sample, scales = "free_y") +
  theme_bw() +
  ggtitle("Distribution of percent.mt per sample")

pdf(file.path(out_dir, "qc_percent_mt_histogram.pdf"), width = 10, height = 6)
print(mt_hist)
dev.off()

###################### QC filtering ##########################
## ---------------------------------------------------
## filter samples with median mitochrodra percentage 
## ---------------------------------------------------

# 1) Compute median percent.mt per sample
mt_stats <- seurat_obj@meta.data %>%
  rownames_to_column("barcode") %>%
  group_by(sample) %>%
  summarize(
    n_cells     = n(),
    mt_median   = median(percent.mt, na.rm = TRUE),
    mt_mean     = mean(percent.mt, na.rm = TRUE)
  ) %>%
  arrange(desc(mt_median))

bad_samples <- mt_stats %>%
  filter(mt_median > mt_threshold) %>%
  pull(sample)

meta <- seurat_obj@meta.data %>%
  rownames_to_column("barcode")

keep_cells <- meta %>%
  filter(!(sample %in% bad_samples)) %>%
  pull(barcode)

seurat_obj <- subset(seurat_obj, cells = keep_cells)


################### QC filtering for cells #################### 

qc_filter_mad_samplewise <- function(
    seurat_obj,
    metric,
    nmads = 3,
    upper_only = FALSE,
    sample_col = "sample"
) {
  # seurat_obj: merged Seurat object
  # metric:     column name in meta.data, e.g. "percent.mt", "log1p_total_count"
  # nmads:      number of MADs for cutoff
  # upper_only: FALSE -> two-sided outlier (like Python function)
  #             TRUE  -> only high outliers
  # sample_col: column in meta.data that stores sample IDs
  
  if (!inherits(seurat_obj, "Seurat")) {
    stop("seurat_obj must be a Seurat object.")
  }
  
  md <- seurat_obj@meta.data %>%
    rownames_to_column("barcode")
  
  if (!metric %in% colnames(md)) {
    stop("Metric '", metric, "' not found in meta.data. Available: ",
         paste(colnames(md), collapse = ", "))
  }
  if (!sample_col %in% colnames(md)) {
    stop("Sample column '", sample_col, "' not found in meta.data.")
  }
  
  # Compute sample-wise median and MAD for the metric
  # Use constant=1 to be closer to typical Python mad() behavior
  stats <- md %>%
    group_by(.data[[sample_col]]) %>%
    summarize(
      median_val = median(.data[[metric]], na.rm = TRUE),
      mad_val    = mad(
        x        = .data[[metric]],
        center   = median(.data[[metric]], na.rm = TRUE),
        constant = 1,
        na.rm    = TRUE
      ),
      .groups = "drop"
    ) %>%
    mutate(
      lower_cutoff = median_val - nmads * mad_val,
      upper_cutoff = median_val + nmads * mad_val
    )
  
  message("Sample-wise ", metric, " median & MAD stats:")
  print(stats)
  
  # Join back to per-cell metadata
  md2 <- md %>%
    left_join(stats, by = sample_col)
  
  x <- md2[[metric]]
  lower <- md2$lower_cutoff
  upper <- md2$upper_cutoff
  
  if (!upper_only) {
    outlier <- (x < lower) | (x > upper)
    direction <- "two-sided"
  } else {
    outlier <- (x > upper)
    direction <- "upper-only"
  }
  
  outlier[is.na(outlier)] <- FALSE
  
  n_before <- nrow(md2)
  n_out    <- sum(outlier)
  n_after  <- n_before - n_out
  
  message("MAD-based QC on metric '", metric, "' (", direction, ", sample-wise)")
  message("  nmads      : ", nmads)
  message("  Cells total: ", n_before)
  message("  Outliers   : ", n_out, " (", round(100 * n_out / n_before, 1), "%)")
  message("  Kept cells : ", n_after, " (", round(100 * n_after / n_before, 1), "%)")
  
  keep_cells <- md2$barcode[!outlier]
  seurat_filtered <- subset(seurat_obj, cells = keep_cells)
  
  return(seurat_filtered)
}

seurat_obj <- qc_filter_mad_samplewise(
  seurat_obj,
  metric     = "log1p_total_count",
  nmads      = 5,
  upper_only = FALSE,  # two-sided; but practically removes low + very high weirdos
  sample_col = "sample"
)

seurat_obj <- qc_filter_mad_samplewise(
  seurat_obj,
  metric     = "nFeature_RNA",
  nmads      = 5,
  upper_only = FALSE,
  sample_col = "sample"
)

low.nFeature_RNA  <- quantile(seurat_obj$nFeature_RNA, 0.01)
high.nFeature_RNA <- quantile(seurat_obj$nFeature_RNA, 0.99)
low.nCount_RNA    <- quantile(seurat_obj$nCount_RNA, 0.01)
high.nCount_RNA   <- quantile(seurat_obj$nCount_RNA, 0.99)


n_before <- ncol(seurat_obj)

# Apply global QC filter
seurat_obj <- subset( seurat_obj,subset = nFeature_RNA > low.nFeature_RNA &
    nFeature_RNA < high.nFeature_RNA &
    nCount_RNA   > low.nCount_RNA &
   nCount_RNA   < high.nCount_RNA &
    percent.mt   < 5)



# cells after
n_after <- ncol(seurat_obj)
message("Global QC filter applied:")
message("Cells before: ", n_before)
message("Cells after : ", n_after,
        " (", round(100 * n_after / n_before, 1), "% kept)")



################### Normalization ####################### 

#seurat_obj_filtered <- NormalizeData(seurat_obj)
#seurat_obj_filtered <- FindVariableFeatures(seurat_obj_filtered)
#seurat_obj_filtered <- ScaleData(seurat_obj_filtered)
#seurat_obj_filtered <- RunPCA(seurat_obj_filtered)
#seurat_obj_filtered <- FindNeighbors(seurat_obj_filtered,dims=1:20)
#seurat_obj_filtered <- FindClusters((seurat_obj_filtered))
#seurat_obj_filtered <- RunUMAP(seurat_obj_filtered,dims =1:20)

#p1 <- DimPlot(seurat_obj_filtered,reduction = 'umap',group.by="sample")
#p2 <- DimPlot(seurat_obj_filtered,reduction = 'umap',group.by="batch")

#grid.arrange(p1,p2, ncol = 2,nrow = 2)


################### sctransform ####################### 

seurat_obj_norm <-SCTransform(seurat_obj,vars.to.regress = "percent.mt",verbose = FALSE)
seurat_obj_norm <- RunPCA(seurat_obj_norm,assay   = "SCT", verbose = TRUE)

seurat_obj_before <- seurat_obj_norm %>%
  RunUMAP(reduction = "pca", dims = 1:20) %>%
  FindNeighbors(reduction = "pca", dims = 1:20) %>%
  FindClusters(resolution = 0.2)

before <- DimPlot(seurat_obj_before, reduction = "umap",group.by  = "batch") + ggtitle("Before Harmony")
############### Run Harmony ####################
seurat_obj_harm <- seurat_obj_norm %>%
  RunHarmony(group.by.vars = "batch", plot_convergence = FALSE)


seurat_obj_harm <- seurat_obj_harm %>% 
  RunUMAP(reduction = 'harmony', dims = 1:20) %>% 
  FindNeighbors(reduction = 'harmony', dims = 1:20) %>% 
  FindClusters(resolution = 0.2)

after <- DimPlot(
  seurat_obj_harm,
  reduction = 'umap',
  group.by = 'batch'
) + ggtitle("After Harmony")


combined <- before | after

ggsave(filename = paste0(out_dir,"/umap_before_after_harmony.pdf"),plot= combined, width    = 12,height   = 6)

seurat_obj_harm$gender <- ifelse(startsWith(seurat_obj_harm$sample, "F"),"Female","Male")

#########################clustering analysis ################
n_clust <- length(unique(seurat_obj_harm$seurat_clusters))
my_colors <- brewer.pal(n = max(3, min(n_clust, 12)), "Set3")   # Set3 works up to 12 colors


clusters <- DimPlot(seurat_obj_harm,reduction = 'umap', cols = my_colors, group.by = 'seurat_clusters') + ggtitle("clustering")
conditions <- DimPlot(seurat_obj_harm,reduction = 'umap', group.by = 'gender') + ggtitle("Gender")

combined_condition <- clusters|conditions

ggsave(filename = paste0(out_dir,"/umap_clusters.pdf"),plot = combined_condition, width    = 12,height   = 6)


########################### Finder All markers ################### 
Idents(seurat_obj_harm) <- "seurat_clusters"
seurat_obj_harm<- PrepSCTFindMarkers(seurat_obj_harm) 

all_cluster_markers <- FindAllMarkers(
  object = seurat_obj_harm,
  assay          = "SCT",
  only.pos = TRUE,         # Only return genes that are more highly expressed in the cluster
  min.pct = 0.25,          # Minimum percentage of cells expressing the gene
  logfc.threshold = 0.25   # Minimum log2 fold change threshold
)

write.csv(all_cluster_markers, file = paste0(out_dir,"DEG_all_clusters.csv"), row.names = FALSE)


top5 <- all_cluster_markers %>% 
  group_by(cluster) %>% 
  top_n(n = 5, wt = avg_log2FC) 

# Extract just the list of gene names for plotting
top5_genes_to_plot <- unique(top5$gene)

seurat_obj_harm <- ScaleData(seurat_obj_harm, assay = "SCT", features = top5_genes_to_plot)
p <- DoHeatmap(
  object   = seurat_obj_harm,
  features = top5_genes_to_plot,
  group.by = "seurat_clusters",
  assay    = "SCT"
) + 
  scale_fill_gradientn(colors = viridis::inferno(256))
ggsave(
  filename = file.path(out_dir, "test_heatmap.pdf"),
  plot     = p,
  width    = 10,
  height   = 10,
  units    = "in")


##### manually check cell marker and adding celltype information ############# 
new_cell_types <- c(
  "0" = "STB",
  "1" = "CTB",
  "2" = "Stromal",
  "3" = "Hofbauer",
  "4" = "EVT",
  "5" = "CTB",
  "6" = "Stromal",
  "7" = "endothelial",
  "8" = "Hofbauer",
  "9" = "endothelial"
)

Idents(seurat_obj_harm) <- "seurat_clusters"
cell_annotations <- data.frame(cell_type = new_cell_types[as.character(Idents(seurat_obj_harm))])
rownames(cell_annotations) <- colnames(seurat_obj_harm)

seurat_obj_harm <- AddMetaData(
  object   = seurat_obj_harm,
  metadata = cell_annotations)

saveRDS(seurat_obj_harm,
        file = paste0(out_seurat_dir,"/seurat_obj_harm_with_celltype.rds"))

############### Cell annotation ############################## 
lapply(c("dplyr","Seurat","HGNChelper","openxlsx"), library, character.only = T)

source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_wrapper.R"); 
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
db_ <- "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
tissue <- "Placenta" 
gs_list <- gene_sets_prepare(db_, tissue)
seurat_package_v5 <- isFALSE('counts' %in% names(attributes(seurat_obj_harm[["SCT"]])));
scRNAseqData_scaled <- if (seurat_package_v5) as.matrix(seurat_obj_harm[["SCT"]]$scale.data) else as.matrix(seurat_obj_harm[["SCT"]]@scale.data)
es.max <- sctype_score(scRNAseqData = scRNAseqData_scaled, scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)
cL_resutls <- do.call("rbind", lapply(unique(seurat_obj_harm@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(seurat_obj_harm@meta.data[seurat_obj_harm@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(seurat_obj_harm@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores <- cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] <- "Unknown"
print(sctype_scores[,1:3])



#############moncle3 ###########
suppressPackageStartupMessages({
  library(Seurat)
  library(monocle3)
  library(SeuratWrappers)})

DefaultAssay(seurat_obj_harm) <- "RNA"

#### cpmvert seurat object into cell_data_set object in monocles 
cds <- as.cell_data_set(seurat_obj_harm)
#cell meta data 
#colData(cds)
#gene meta data 
#fData(cds)
fData(cds)$gene_short_name <-rownames(fData(cds))
# to get counts of cell_data_set 
#counts(cds)

recreate.partition <- rep(1, length(colnames(cds)))
names(recreate.partition) <- colnames(cds)
recreate.partition <- as.factor(recreate.partition)

cds@clusters$UMAP$partitions <- recreate.partition

list_cluster <- seurat_obj_harm@active.ident
cds@clusters$UMAP$clusters <- list_cluster

reducedDims(cds)$UMAP <- Embeddings(seurat_obj_harm, "umap")
cluster.before.trajectory <- plot_cells(
  cds,
  color_cells_by = "cluster",
  label_groups_by_cluster = FALSE,
  group_label_size = 5
) +
  theme(legend.position = "right")

ggsave(filename = paste0(out_dir,"/cluster_before_tra.pdf"),plot = cluster.before.trajectory, width    = 12,height   = 6)


cds <- preprocess_cds(cds, num_dim = 50)
cds <- cluster_cells(cds)

colData(cds)$redefined_cluster <- colData(cds)$cell_type
colData(cds)$redefined_cluster <- as.factor(colData(cds)$redefined_cluster)
cluster.names <- plot_cells(
  cds,
  color_cells_by = "redefined_cluster",
  label_groups_by_cluster = FALSE,
  group_label_size = 3
) +theme(legend.position = "right")

ggsave(filename = paste0(out_dir,"/cluster_before_tra_celltype.pdf"),plot = cluster.names , width    = 12,height   = 6)



cds <- learn_graph(cds,use_partition = FALSE)


root_cells <- colnames(cds)[colData(cds)$redefined_cluster == "CTB"]

# 
# plot_cells(cds, color_cells_by = "seurat_clusters",
#            label_groups_by_cluster = FALSE,
#            label_leaves = FALSE,
#            label_branch_points = FALSE,
#            label_roots = FALSE,
#            group_label_size = 3)
# 


cds <- order_cells(cds,reduction_method = 'UMAP', root_cells = root_cells)
pseudotime_plt<-plot_cells(cds, color_cells_by = "pseudotime",
           label_groups_by_cluster = FALSE,
           label_leaves = FALSE,
           label_branch_points = FALSE,
           label_roots = FALSE)
ggsave(filename = paste0(out_dir,"/pseudotime_plt.pdf"),plot = pseudotime_plt , width    = 12,height   = 6)

cds$monocle3_pseudotime <- pseudotime(cds)
pseudo_data <- as.data.frame(colData(cds)) 

pseudo_time_cluster <- ggplot(pseudo_data,aes(monocle3_pseudotime, redefined_cluster, fill=redefined_cluster)) +
  geom_boxplot() 

ggsave(filename = paste0(out_seurat_dir,"/average_pseudotime.pdf"),plot = pseudo_time_cluster, width    = 12,height   = 6)


############################ Marker exploration for gender specific ###################


Idents(seurat_obj_harm) <- "seurat_clusters"
seurat_obj_harm<- PrepSCTFindMarkers(seurat_obj_harm) 


marker_cluster1 <- FindConservedMarkers(seurat_obj_harm, 
                  ident.1 = 1,
                  grouping.var = "gender",
                  assay = "SCT",
                  slot = "data",
                  only.pos = TRUE,
                  min.pct = 0.25,
                  logfc.threshold = 0.25)
head(marker_cluster1)


