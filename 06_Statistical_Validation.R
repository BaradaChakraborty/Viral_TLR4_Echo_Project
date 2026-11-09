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