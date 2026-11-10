# 09_Receptor_Statistical_Validation.R
library(Seurat)

# 1. Load my isolated microglia Seurat object
microglia <- readRDS("data/microglia_subset.rds") 

# 2. Define my exact diagnostic genes
diagnostic_genes <- c("TLR4", "LY96", "CD14", "MYD88", "TICAM1", 
                      "NLRP3", "IL1B", "CGAS", "MB21D1", "STING1", "TMEM173")

# 2.5 Join the fragmented data layers
microglia <- JoinLayers(microglia)

# 3. Filter the list to ONLY include genes that actually exist in your dataset
valid_genes <- intersect(diagnostic_genes, rownames(microglia))
print(paste("Valid genes found for testing:", paste(valid_genes, collapse=", ")))

# 4. Now force Seurat to calculate statistics on the valid genes
receptor_stats <- FindMarkers(microglia, 
                              ident.1 = "18", 
                              ident.2 = c("4", "19"), 
                              features = valid_genes,
                              logfc.threshold = 0,
                              min.pct = 0,
                              test.use = "wilcox")

# 5. Format the table for publication
receptor_stats$Gene <- rownames(receptor_stats)
receptor_stats <- receptor_stats[, c("Gene", "avg_log2FC", "p_val", "p_val_adj", "pct.1", "pct.2")]

# 6. Save the definitive mathematical proof to your tables folder
write.csv(receptor_stats, "results/tables/Upstream_Receptor_Stats.csv", row.names = FALSE)

# 7. Print to console to view immediately
print(receptor_stats)