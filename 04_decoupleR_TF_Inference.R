library(Seurat)
library(decoupleR)
library(OmnipathR)
library(dplyr)
library(tidyr)
library(ggplot2)
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c("decoupleR", "OmnipathR", "tidyr"))
library(Seurat)
library(decoupleR)
library(OmnipathR)
library(dplyr)
library(tidyr)
library(ggplot2)

microglia_obj <- readRDS("E:/Viral_TLR4_Echo_Project/data/microglia_subset.rds")

# Download the CollecTRI Gene Regulatory Network
# CollecTRI is the most comprehensively curated database of human Transcription Factor-target interactions.
message("Downloading CollecTRI network...")
net <- get_collectri(organism = 'human', split_complexes = FALSE)
# Extract the normalized RNA matrix from the Seurat object
# decoupleR requires a raw matrix of cells (columns) and genes (rows)
mat <- as.matrix(GetAssayData(microglia_obj, layer = "data"))



# Join the v5 layers and extract the normalized RNA matrix
# We must collapse the 12 patient-specific layers into one unified matrix for decoupleR
microglia_obj <- JoinLayers(microglia_obj)
mat <- as.matrix(GetAssayData(microglia_obj, assay = "RNA", layer = "data"))


# Run decoupleR to infer Transcription Factor (TF) activity
message("Inferring TF activities. This may take a moment...")
tf_acts <- run_ulm(mat = mat, network = net, .source = "source", .target = "target", .mor = "mor", minsize = 5)

# Reshape the mathematical output into a matrix
tf_mat <- tf_acts %>%
  select(statistic, source, condition) %>%
  pivot_wider(names_from = condition, values_from = statistic) %>%
  tibble::column_to_rownames("source") %>%
  as.matrix()

# Inject the TF activity back into the Seurat object
microglia_obj[["tf_ulm"]] <- CreateAssayObject(data = tf_mat)
DefaultAssay(microglia_obj) <- "tf_ulm"

message("decoupleR Inference Complete! The default assay is now TF Activity.")


# 6. Let's Reshape the mathematical output into a matrix (Edit)
# CRITICAL FIX: Let's extract 'score' (the numeric inference) instead of 'statistic' (the text label)
tf_mat <- tf_acts %>%
  select(score, source, condition) %>%
  pivot_wider(names_from = condition, values_from = score) %>%
  tibble::column_to_rownames("source") %>%
  as.matrix()

# 7. Let's inject the numeric TF activity back into the Seurat object
microglia_obj[["tf_ulm"]] <- CreateAssayObject(data = tf_mat)
DefaultAssay(microglia_obj) <- "tf_ulm"

message("Numeric TF matrix successfully injected. Ready for visualization.")



# 1. Reload required libraries
library(Seurat)
library(decoupleR)
library(OmnipathR)
library(dplyr)
library(tidyr)
library(ggplot2)

# 2. Load the purified microglia dataset
microglia_obj <- readRDS("E:/Viral_TLR4_Echo_Project/data/microglia_subset.rds")

# 3. Download the CollecTRI Gene Regulatory Network
message("Downloading CollecTRI network...")
net <- get_collectri(organism = 'human', split_complexes = FALSE)

# 4. Join the v5 layers and extract the normalized RNA matrix
# We must collapse the 12 patient-specific layers into one unified matrix for decoupleR
microglia_obj <- JoinLayers(microglia_obj)
mat <- as.matrix(GetAssayData(microglia_obj, assay = "RNA", layer = "data"))

# 5. Run decoupleR to infer Transcription Factor (TF) activity
message("Inferring TF activities. This may take a moment...")
tf_acts <- run_ulm(mat = mat, network = net, .source = "source", .target = "target", .mor = "mor", minsize = 5)

# 6. Reshape the mathematical output into a matrix
# CRITICAL FIX: Extract 'score' (numeric inference) instead of 'statistic' (text label)
tf_mat <- tf_acts %>%
  select(score, source, condition) %>%
  pivot_wider(names_from = condition, values_from = score) %>%
  tibble::column_to_rownames("source") %>%
  as.matrix()

# 7. Inject the numeric TF activity back into the Seurat object
microglia_obj[["tf_ulm"]] <- CreateAssayObject(data = tf_mat)
DefaultAssay(microglia_obj) <- "tf_ulm"

message("Numeric TF matrix successfully injected. Ready for visualization.")