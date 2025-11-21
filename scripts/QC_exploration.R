

############ QC and marker gene exploration ################# 
plot_and_save_feature <- function(obj, gene, out_dir, width = 12, height = 6) {
  p <- FeaturePlot(
    obj,
    features = gene,
    cols = c("lightgray", "blue")
  ) + ggtitle(gene)
  
  ggsave(
    filename = file.path(out_dir, paste0(gene, ".pdf")),
    plot = p,
    width = width,
    height = height
  )
}



########vlnplots for some genes ######### 
seurat_obj_harm <- subset(seurat_obj_harm, idents = setdiff(levels(seurat_obj_harm), NA))

genes <- c("PAPPA2", "CYP19A1")

p_list <- VlnPlot(
  seurat_obj_harm,
  features = genes,
  group.by = "seurat_clusters",
  pt.size = 0.3,
  combine = FALSE
)

p_list[[1]] <- p_list[[1]] + ggtitle("PAPPA2")
p_list[[2]] <- p_list[[2]] + ggtitle("CYP19A1")

combined <- p_list[[1]] + p_list[[2]] + plot_layout(ncol = 2)

ggsave(filename = paste0(out_dir,"/Differential_expressed_genes.pdf"),plot = combined, width    = 12,height   = 6)

XIST_vln <- VlnPlot(
  seurat_obj_harm,
  features = "XIST",
  group.by = "gender",
  pt.size = 0.1
) + NoLegend()

ggsave(filename = paste0(out_dir,"/XIST_Vln.pdf"),plot = XIST_vln, width    = 12,height   = 6)

DefaultAssay(seurat_obj_harm) <- "SCT"

################## Xist exploration #################### 
male_only <- subset(seurat_obj_harm, subset = gender == "Male")
female_only <- subset(seurat_obj_harm, subset = gender == "Female")

male_xist <- FeaturePlot(male_only,features = "XIST",cols = c("lightgray", "blue"))
female_xist <- FeaturePlot(female_only,features = "XIST",cols = c("lightgray", "blue"))
Xist_plotting <- male_xist |female_xist
ggsave(paste0(filename = out_dir,"/XIST_male_female.pdf"),plot = Xist_plotting, width    = 12,height   = 6)

xist_male_cells <- WhichCells(
  seurat_obj_harm,
  expression = gender == "Male" & XIST > 0)
meta_sub <- seurat_obj_harm@meta.data[xist_male_cells, ]

# count number per sample
table(meta_sub$sample)

genes_to_plot <- c(
  "HLA-G",
  "CYP19A1",
  "EPS8L1",
  "NOTUM",
  "SLC6A4",
  "SLC27A2",
  "KRT7",
  "CSH1"
)
for (g in genes_to_plot) {
  plot_and_save_feature(seurat_obj_harm, g, out_dir)
}