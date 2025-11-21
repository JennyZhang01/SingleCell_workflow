#############monocle3 ###########
suppressPackageStartupMessages({
  library(Seurat)
  library(monocle3)
  library(SeuratWrappers)})

## ---------------------------
## Parse arguments
## ---------------------------

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: seurat_cluster.R project_dir mt_threhold\n",
       "Example:\n",
       "  Rscript monocle3.R \\\n",
       "  /rhome/jzhan413/bigdata/proj/alleleRNA/24-12-17_SingleCell/25-06-10_all_samples/ \\\n")
}

project_dir     <- args[1]
out_dir <- "monocle/"
setwd(project_dir)
if (!dir.exists(out_dir)){
  dir.create(out_dir,recursive = TRUE)
} 


seurat_dir <- "Processed_data/"
seurat_obj_harm <- readRDS(paste0(seurat_dir,"/seurat_obj_harm_with_celltype.rds"))



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

