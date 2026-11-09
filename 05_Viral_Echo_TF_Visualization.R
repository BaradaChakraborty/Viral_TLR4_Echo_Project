# 1. Ensure the default assay is set to the decoupleR TF activity scores
DefaultAssay(microglia_obj) <- "tf_ulm"
# 2. Select the Master Regulators of the Viral-Echo (TLR4/MyD88) Pathway
# RELA (NF-kB p65) = The master inflammatory switch downstream of MyD88
# STAT1 = The master anti-viral / interferon response regulator
# IRF1 = Interferon Regulatory Factor 1 (classic viral response)
# SPI1 (PU.1) = The core microglial survival factor (used as a biological control)
viral_echo_tfs <- c("RELA", "STAT1", "IRF1", "SPI1")
# Load patchwork to safely handle the plot_annotation() additions
library(patchwork)
# 3. Spatial Visualization (FeaturePlot)
# Let's remove the q5/q95 cutoffs because decoupleR outputs continuous, normally-distributed statistical scores, not sparse integer counts.
tf_spatial_plot <- FeaturePlot(microglia_obj, features = viral_echo_tfs, 
                               pt.size = 1, ncol = 2) + 
  plot_annotation(title = "Inferred Transcription Factor Activity (Viral-Echo Pathway)")
print(tf_spatial_plot)
# 4. Quantitative Statistical Visualization (Violin Plot)
tf_vln_plot <- VlnPlot(microglia_obj, features = c("RELA", "STAT1"), 
                       pt.size = 0.1, ncol = 2) + 
  plot_annotation(title = "Quantitative TF Activation per Cluster")

print(tf_vln_plot)

install.packages("withr", force = TRUE)
library(withr)
# 4. Quantitative Statistical Visualization (Violin Plot)
tf_vln_plot <- VlnPlot(microglia_obj, features = c("RELA", "STAT1"), 
                       pt.size = 0.1, ncol = 2) + 
  plot_annotation(title = "Quantitative TF Activation per Cluster")

print(tf_vln_plot)

# 5. Save the conclusive figures for your portfolio
fig_dir <- "E:/Viral_TLR4_Echo_Project/results/figures"
ggsave(filename = file.path(fig_dir, "04_Viral_Echo_TF_Spatial.png"), plot = tf_spatial_plot, width = 10, height = 8, dpi = 300, bg = "white")
ggsave(filename = file.path(fig_dir, "05_Viral_Echo_TF_Violin.png"), plot = tf_vln_plot, width = 10, height = 6, dpi = 300, bg = "white")

message("Conclusive TF Activity plots generated and saved!")
