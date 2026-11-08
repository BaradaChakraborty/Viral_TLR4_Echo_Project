# 1. Isolate the microglial clusters from the global atlas
microglia_obj <- subset(pd_filtered, idents = c(4, 18, 19))

# 2. Re-process the isolated microglia
# This is a critical mathematical step. We must recalculate the variance now that the neurons are gone.
microglia_obj <- FindVariableFeatures(microglia_obj, selection.method = "vst", nfeatures = 2000)
microglia_obj <- ScaleData(microglia_obj)
microglia_obj <- RunPCA(microglia_obj, features = VariableFeatures(object = microglia_obj))
microglia_obj <- RunUMAP(microglia_obj, dims = 1:15) # We use fewer dimensions (15) because the biological complexity is now much lower

# 3. Visualize the purified microglial dataset
library(ggplot2)
microglia_only_umap <- DimPlot(microglia_obj, reduction = "umap", label = TRUE, pt.size = 1) + 
  ggtitle("Purified Microglial Atlas (Clusters 4, 18, 19)")
print(microglia_only_umap)

# 4. Save the purified object to disk so you never have to run the pipeline from scratch again
saveRDS(microglia_obj, file = "E:/Viral_TLR4_Echo_Project/data/microglia_subset.rds")
message("Microglia successfully isolated and saved as an .rds file.")

# 5. Save the Purified UMAP for your portfolio
fig_dir <- "E:/Viral_TLR4_Echo_Project/results/figures"
ggsave(filename = file.path(fig_dir, "03_Purified_Microglia_UMAP.png"), 
       plot = microglia_only_umap, 
       width = 8, height = 6, dpi = 300, bg = "white")

message("Purified Microglial UMAP successfully saved to figures directory.")