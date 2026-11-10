# 10_CD14_vs_cGAS_BarPlot.R
library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)

# 1. Load your updated Seurat object
microglia <- readRDS("data/microglia_subset.rds")

# 2. Extract the raw expression data for our two diagnostic genes
expr_data <- FetchData(microglia, vars = c("CD14", "MB21D1"))
expr_data$Cluster <- ifelse(Idents(microglia) == "18", "Cluster 18", "Homeostatic (4 & 19)")
expr_data$Cluster <- factor(expr_data$Cluster, levels = c("Homeostatic (4 & 19)", "Cluster 18"))

# 3. Calculate Mean and Standard Error of the Mean (SEM) for clean plotting
summary_data <- expr_data %>%
  group_by(Cluster) %>%
  summarise(
    CD14_mean = mean(CD14), 
    CD14_se = sd(CD14)/sqrt(n()),
    MB21D1_mean = mean(MB21D1), 
    MB21D1_se = sd(MB21D1)/sqrt(n())
  )

# 4. Create the CD14 Bar Plot (The TLR4 Smoking Gun)
p1 <- ggplot(summary_data, aes(x = Cluster, y = CD14_mean, fill = Cluster)) +
  geom_bar(stat = "identity", color = "black", width = 0.6) +
  geom_errorbar(aes(ymin = CD14_mean - CD14_se, ymax = CD14_mean + CD14_se), width = 0.2, linewidth = 0.8) +
  scale_fill_manual(values = c("Homeostatic (4 & 19)" = "lightgrey", "Cluster 18" = "#00008B")) +
  labs(title = "CD14 (TLR4 Co-Receptor)",
       subtitle = "Wilcoxon padj = 2.46e-30 ***",
       y = "Mean Expression Level") +
  theme_classic() +
  theme(legend.position = "none", 
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(color = "darkred", face = "bold", size = 12))

# 5. Create the MB21D1 Bar Plot (The cGAS Cytosolic Sensor)
p2 <- ggplot(summary_data, aes(x = Cluster, y = MB21D1_mean, fill = Cluster)) +
  geom_bar(stat = "identity", color = "black", width = 0.6) +
  geom_errorbar(aes(ymin = MB21D1_mean - MB21D1_se, ymax = MB21D1_mean + MB21D1_se), width = 0.2, linewidth = 0.8) +
  scale_fill_manual(values = c("Homeostatic (4 & 19)" = "lightgrey", "Cluster 18" = "#1E90FF")) +
  labs(title = "MB21D1 (cGAS Sensor)",
       subtitle = "Wilcoxon padj = 1.0 (ns)",
       y = "Mean Expression Level") +
  theme_classic() +
  theme(legend.position = "none", 
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(color = "grey40", face = "bold", size = 12))

# 6. Combine the plots side-by-side
final_plot <- p1 | p2
final_plot <- final_plot + plot_annotation(
  title = "Figure 8: Statistical Divergence of Surface vs. Cytosolic Sensors",
  theme = theme(plot.title = element_text(face = "bold", size = 16))
)

# 7. Save the definitive, clean statistical figure
ggsave("results/figures/08_CD14_vs_cGAS_BarPlot.png", plot = final_plot, width = 9, height = 5, dpi = 300)
