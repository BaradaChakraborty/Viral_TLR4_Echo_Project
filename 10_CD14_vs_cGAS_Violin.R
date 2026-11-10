# 10_CD14_vs_cGAS_Violin.R
library(Seurat)
library(ggplot2)
library(patchwork) # For combining plots side-by-side

# 1. Load your updated Seurat object
microglia <- readRDS("data/microglia_subset.rds")

# 2. Create a clean binary grouping for the plot ("Cluster 18" vs "Rest")
microglia$Comparison_Group <- ifelse(Idents(microglia) == "18", "Cluster 18", "Homeostatic (4 & 19)")
Idents(microglia) <- "Comparison_Group"
# Ensure the order on the x-axis is logical
microglia$Comparison_Group <- factor(microglia$Comparison_Group, levels = c("Homeostatic (4 & 19)", "Cluster 18"))

# 3. Create the CD14 Violin Plot (The TLR4 Smoking Gun)
# I manually add the exact p-value you derived from FindMarkers
p1 <- VlnPlot(microglia, features = "CD14", pt.size = 0.5) +
  scale_fill_manual(values = c("Homeostatic (4 & 19)" = "lightgrey", "Cluster 18" = "#00008B")) +
  labs(title = "CD14 (TLR4 Co-Receptor)",
       subtitle = "Wilcoxon padj = 2.46e-30 ***",
       y = "Expression Level") +
  theme(legend.position = "none", 
        axis.title.x = element_blank(),
        plot.subtitle = element_text(color = "darkred", face = "bold"))

# 4. Create the MB21D1 Violin Plot (The cGAS Cytosolic Sensor)
p2 <- VlnPlot(microglia, features = "MB21D1", pt.size = 0.5) +
  scale_fill_manual(values = c("Homeostatic (4 & 19)" = "lightgrey", "Cluster 18" = "#1E90FF")) +
  labs(title = "MB21D1 (cGAS Sensor)",
       subtitle = "Wilcoxon padj = 1.0 (ns)",
       y = "Expression Level") +
  theme(legend.position = "none", 
        axis.title.x = element_blank(),
        plot.subtitle = element_text(color = "grey40", face = "bold"))

# 5. Combine the plots side-by-side using patchwork
final_plot <- p1 | p2
final_plot <- final_plot + plot_annotation(
  title = "Figure 7: Statistical Divergence of Surface vs. Cytosolic Sensors",
  theme = theme(plot.title = element_text(face = "bold", size = 16))
)

# 6. Save the definitive statistical figure
ggsave("results/figures/08_CD14_vs_cGAS_Violin.png", plot = final_plot, width = 9, height = 5, dpi = 300)