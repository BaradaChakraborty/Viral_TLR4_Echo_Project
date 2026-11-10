# 08_Receptor_Complex_Validation.R
library(Seurat)
library(ggplot2)
microglia <- readRDS("data/microglia_subset.rds")

# 2. Define the diagnostic genes using alternative aliases for the cytosolic sensors
diagnostic_genes <- c(
  "TLR4", "LY96", "CD14",    # The Surface Sensing Complex 
  "MYD88", "TICAM1",         # The TLR4 Adaptors 
  "NLRP3", "IL1B",           # TLR4 Signal 1 
  "MB21D1", "TMEM173"        # The Cytosolic Sensors (CGAS and STING1)
)

# 3. Generate a DotPlot to visualize both expression level and percentage of active cells
validation_plot <- DotPlot(microglia, features = diagnostic_genes, cluster.idents = TRUE) +
  theme_minimal() +
  scale_color_gradientn(colors = c("lightgrey", "#1E90FF", "#00008B")) +
  labs(
    title = "Figure 6: Upstream Structural Validation of the Receptor Complex",
    subtitle = "Contrasting surface-level TLR4 machinery with cytosolic cGAS-STING sensors",
    x = "Diagnostic Target Genes",
    y = "Microglial Cluster",
    color = "Average\nExpression",
    size = "Percent\nExpressed"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold.italic", size = 12, color = "black"),
    axis.text.y = element_text(face = "bold", size = 12, color = "black"),
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 12, color = "grey30"),
    panel.grid.major.x = element_blank()
  )

# 4. Save the figure to your results folder
ggsave("results/figures/07_Upstream_Receptor_DotPlot.png", plot = validation_plot, width = 10, height = 5, dpi = 300)