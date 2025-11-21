#############moncle3 ###########
suppressPackageStartupMessages({
  library(Seurat)
  library(monocle3)
  library(SeuratWrappers)})

#todo adding reading seurat_object ############ 
seurat_dir <- "Processed_data/"
seurat_obj_harm <- readRDS(paste0(seurat_dir,"/seurat_obj_harm_with_celltype.rds"))


DefaultAssay(seurat_obj_harm) <- "RNA"

cds <- as.cell_data_set(seurat_obj_harm)
cds <- preprocess_cds(cds, num_dim = 50)
cds <- cluster_cells(cds)
cds <- learn_graph(cds)

root_cells <- colnames(cds)[colData(cds)$seurat_clusters == "0"]


cds <- order_cells(cds, root_cells = root_cells)
plot_cells(cds, color_cells_by = "pseudotime", show_trajectory_graph = TRUE)
plot_cells(cds, color_cells_by = "seurat_clusters",
           label_groups_by_cluster = FALSE,
           label_leaves = TRUE,
           label_branch_points = TRUE)




