#Normalize and identify highly variable genes
pd_filtered <- NormalizeData(pd_filtered, normalization.method = "LogNormalize", scale.factor = 10000)
pd_filtered <- FindVariableFeatures(pd_filtered, selection.method = "vst", nfeatures = 2000)
#Scale the data (center and standardize)
pd_filtered <- ScaleData(pd_filtered)
#Linear Dimesionality Reduction (PCA)
pd_filtered <- RunPCA(pd_filtered, features = VariableFeatures(object = pd_filtered))
#Construct the K-Nearest Neighbor (KNN) graph and Cluster
#Using the top 30 Principal Components is standard for complex human brain tissue
pd_filtered <- FindNeighbors(pd_filtered, dims = 1:30)
pd_filtered <- FindClusters(pd_filtered, resolution = 0.5)
#Non-linear Dimensionality Reduction (UMAP)
pd_filtered <- RunUMAP(pd_filtered, dims = 1:30)
#Visualization 1: The Global Atlas
library(ggplot2)
global_umap <- DimPlot(pd_filtered, reduction = "umap", label = TRUE, pt.size = 0.5) + 
  ggtitle("Global Substantia Nigra Atlas") +
  NoLegend()
print(global_umap)
print(global_umap)

#Visualization 2: Locating the Microglia
#Using canonical human microglial markers to light up their specific cluster
microglia_markers <- c("CSF1R", "CX3CR1", "P2RY12", "C3")
microglia_umap <- FeaturePlot(pd_filtered, features = microglia_markers, pt.size = 0.5, ncol = 2)
print(microglia_umap)
fig_dir <- "E:/Viral_TLR4_Echo_Project/results/figures"
# Defensively create the folder if it does not already exist
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
# Save Plot 1: The Global Atlas
ggsave(filename = file.path(fig_dir, "01_Global_UMAP_Atlas.png"), 
       plot = global_umap, 
       width = 8, height = 6, dpi = 300, bg = "white")

# Save Plot 2: The Microglial Markers
ggsave(filename = file.path(fig_dir, "02_Microglia_Markers_FeaturePlot.png"), 
       plot = microglia_umap, 
       width = 10, height = 8, dpi = 300, bg = "white")

message("Plots successfully saved at 300 DPI in: ", fig_dir)