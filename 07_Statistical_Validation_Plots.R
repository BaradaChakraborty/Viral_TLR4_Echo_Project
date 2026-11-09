# 1. Load required libraries
library(Seurat)
library(ggplot2)
library(dplyr)

# 3. Extract the activity scores for our target TFs
target_tfs <- c("IRF1", "STAT1", "SPI1", "RELA")
tf_data <- as.data.frame(t(as.matrix(GetAssayData(microglia_obj, assay = "tf_ulm")[target_tfs, ])))
tf_data$Cluster <- ifelse(Idents(microglia_obj) == 18, "Cluster 18", "Rest (4 & 19)")

# 4. Calculate mean and standard error for error bars
bar_data <- tf_data %>%
  tidyr::pivot_longer(cols = all_of(target_tfs), names_to = "TF", values_to = "Score") %>%
  group_by(TF, Cluster) %>%
  summarise(Mean = mean(Score), 
            SE = sd(Score) / sqrt(n()), 
            .groups = "drop")

# Ensure the TFs are ordered logically in the plot
bar_data$TF <- factor(bar_data$TF, levels = c("IRF1", "STAT1", "SPI1", "RELA"))

# 5. Build the Standalone, Publication-Ready Bar Plot
final_plot <- ggplot(bar_data, aes(x = TF, y = Mean, fill = Cluster)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7, color = "black") +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), position = position_dodge(width = 0.8), width = 0.25) +
  scale_fill_manual(values = c("Cluster 18" = "#00B050", "Rest (4 & 19)" = "#A6A6A6")) +
  theme_classic(base_size = 16) +
  labs(title = "Figure 5: Quantitative Validation of the TLR4 'Viral-Echo'", 
       subtitle = "Inferred Protein Activity in Parkinson's Disease Microglia",
       y = "Mean TF Activity Score", 
       x = "Transcription Factors") +
  theme(legend.position = "right", 
        plot.title = element_text(face = "bold", size = 18),
        plot.subtitle = element_text(size = 14, color = "grey30", margin = margin(b = 15))) +
  
  # Expand Y-axis to make room for brackets
  coord_cartesian(ylim = c(0, 5.0)) +
  
  # IRF1 Bracket & P-value
  annotate("segment", x = 0.8, xend = 1.2, y = 2.1, yend = 2.1, color = "black") +
  annotate("segment", x = 0.8, xend = 0.8, y = 2.0, yend = 2.1, color = "black") +
  annotate("segment", x = 1.2, xend = 1.2, y = 2.0, yend = 2.1, color = "black") +
  annotate("text", x = 1, y = 2.25, label = "p = 2.85e-64", size = 4.5, fontface = "bold") +
  
  # STAT1 Bracket & P-value
  annotate("segment", x = 1.8, xend = 2.2, y = 4.2, yend = 4.2, color = "black") +
  annotate("segment", x = 1.8, xend = 1.8, y = 4.1, yend = 4.2, color = "black") +
  annotate("segment", x = 2.2, xend = 2.2, y = 4.1, yend = 4.2, color = "black") +
  annotate("text", x = 2, y = 4.35, label = "p = 6.58e-36", size = 4.5, fontface = "bold") +
  
  # SPI1 Bracket & P-value
  annotate("segment", x = 2.8, xend = 3.2, y = 4.6, yend = 4.6, color = "black") +
  annotate("segment", x = 2.8, xend = 2.8, y = 4.5, yend = 4.6, color = "black") +
  annotate("segment", x = 3.2, xend = 3.2, y = 4.5, yend = 4.6, color = "black") +
  annotate("text", x = 3, y = 4.75, label = "p = 1.78e-68", size = 4.5, fontface = "bold") +
  
  # RELA Bracket & P-value
  annotate("segment", x = 3.8, xend = 4.2, y = 4.0, yend = 4.0, color = "black") +
  annotate("segment", x = 3.8, xend = 3.8, y = 3.9, yend = 4.0, color = "black") +
  annotate("segment", x = 4.2, xend = 4.2, y = 3.9, yend = 4.0, color = "black") +
  annotate("text", x = 4, y = 4.15, label = "p = 4.14e-09", size = 4.5, fontface = "bold")

# 6. Print and Save
print(final_plot)

fig_dir <- "E:/Viral_TLR4_Echo_Project/results/figures"
ggsave(filename = file.path(fig_dir, "06_Cluster_18_Final_Statistical_BarPlot.png"), 
       plot = final_plot, width = 9, height = 7, dpi = 300, bg = "white")

message("Final validated statistical figure generated and saved!")