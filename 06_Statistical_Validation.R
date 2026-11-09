# 1. Ensure the active assay is set to the decoupleR TF activity scores
DefaultAssay(microglia_obj) <- "tf_ulm"

# 2. Run the Statistical Test (Wilcoxon Rank Sum)
# Question: "Is Cluster 18 mathematically distinct from Clusters 4 and 19?"
message("Running Wilcoxon Rank Sum test on TF activity scores...")
cluster_18_stats <- FindMarkers(microglia_obj, 
                                ident.1 = 18, 
                                ident.2 = c(4, 19), 
                                test.use = "wilcox", 
                                min.pct = 0,          # TF scores exist for all cells, so we don't filter by percentage
                                logfc.threshold = 0)  # We want to calculate the p-value for all TFs

# 3. Isolate the statistics for our specific Viral-Echo hypothesis
viral_echo_tfs <- c("RELA", "STAT1", "IRF1", "SPI1")
viral_echo_results <- cluster_18_stats[rownames(cluster_18_stats) %in% viral_echo_tfs, ]

# Print the definitive mathematical proof to the console
print("--- STATISTICAL PROOF: VIRAL-ECHO TFs IN CLUSTER 18 ---")
print(viral_echo_results)

# 4. Save the full statistical table as a CSV for your manuscript / portfolio
dir.create("E:/Viral_TLR4_Echo_Project/results/tables", recursive = TRUE, showWarnings = FALSE)
write.csv(cluster_18_stats, "E:/Viral_TLR4_Echo_Project/results/tables/Cluster_18_vs_Rest_TF_Stats.csv")

message("Statistical analysis complete! CSV table saved.")
###
# 1. Load the corrupted stats table
stats_df <- read.csv("E:/Viral_TLR4_Echo_Project/results/tables/Cluster_18_vs_Rest_TF_Stats.csv", row.names = 1)

# 2. Extract the raw numerical activity matrix
DefaultAssay(microglia_obj) <- "tf_ulm"
tf_matrix <- as.matrix(GetAssayData(microglia_obj, assay = "tf_ulm"))

# 3. Calculate the true mathematical means for both groups
mean_18 <- rowMeans(tf_matrix[, Idents(microglia_obj) == 18])
mean_rest <- rowMeans(tf_matrix[, Idents(microglia_obj) != 18])

# 4. Calculate the raw Mean Difference and match it to the stats table
mean_diff_all <- mean_18 - mean_rest
stats_df$True_Mean_Diff <- mean_diff_all[rownames(stats_df)]

# 5. Clean up the columns and replace the corrupted Log2FC
stats_df <- stats_df[, c("True_Mean_Diff", "p_val", "p_val_adj", "pct.1", "pct.2")]

# 6. Save the publication-ready table
write.csv(stats_df, "E:/Viral_TLR4_Echo_Project/results/tables/Final_Corrected_TF_Stats.csv")
message("Corrected Statistical Table saved successfully!")