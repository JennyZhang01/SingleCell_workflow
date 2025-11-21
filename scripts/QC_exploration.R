

############ QC and marker gene exploration ################# 

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

HLA_G_umap <- FeaturePlot(seurat_obj_harm,features = "HLA-G",cols = c("lightgray", "blue"))
ggsave(filename = paste0(out_dir,"/HLA-g.pdf"),plot = HLA_G_umap, width    = 12,height   = 6)



CYP19A1_umap <- FeaturePlot(seurat_obj_harm,features = "CYP19A1",cols = c("lightgray", "blue"))
ggsave(filename = paste0(out_dir,"/CYP19A.pdf"),plot =CYP19A1_umap, width    = 12,height   = 6)

EPS8L1_umap <- FeaturePlot(seurat_obj_harm,features = "EPS8L1",cols = c("lightgray", "blue"))
ggsave(filename = paste0(out_dir,"/EPS8L1.pdf"),plot =EPS8L1_umap, width    = 12,height   = 6)


NOTUM_umap <- FeaturePlot(seurat_obj_harm,features = "NOTUM",cols = c("lightgray", "blue"))
ggsave(filename = paste0(out_dir,"/NOTUM.pdf"),plot =NOTUM_umap, width    = 12,height   = 6)


SLC6A4_umap <- FeaturePlot(seurat_obj_harm,features = "SLC6A4",cols = c("lightgray", "blue"))
ggsave(filename = paste0(out_dir,"/SLC6A4.pdf"),plot =SLC6A4_umap, width    = 12,height   = 6)


SLC27A2_umap <- FeaturePlot(seurat_obj_harm,features = "SLC27A2",cols = c("lightgray", "blue"))
ggsave(filename = paste0(out_dir,"/SLC27A2.pdf"),plot = SLC27A2_umap, width    = 12,height   = 6)


KRT7_umap <- FeaturePlot(seurat_obj_harm,features = "KRT7",cols = c("lightgray", "blue"))
ggsave(filename = paste0(out_dir,"/KRT7.pdf"),plot = KRT7_umap, width    = 12,height   = 6)


CSH1_umap <- FeaturePlot(seurat_obj_harm,features = "CSH1",cols = c("lightgray", "blue"))
ggsave(filename = paste0(out_dir,"/SLC27A2.pdf"),plot =CSH1_umap, width    = 12,height   = 6)
