# 1. Set the directory
data_dir <- "E:/Viral_TLR4_Echo_Project/data/GSE140231_RAW"

# 2. Find all the sample files by looking for "GSM"
archive_files <- list.files(data_dir, pattern = "GSM", full.names = TRUE)

# 3. Extract each one into its own dedicated folder
for (file in archive_files) {
  # Strip away all extensions (.tar.gz, .tar, etc.) for a clean folder name
  folder_name <- sub("\\..*", "", basename(file))
  out_path <- file.path(data_dir, folder_name)
  
  # Create the folder and extract
  dir.create(out_path, showWarnings = FALSE)
  untar(tarfile = file, exdir = out_path)
  
  message(paste("Extracted:", folder_name))
}
# 1. Set the directory
data_dir <- "E:/Viral_TLR4_Echo_Project/data/GSE140231_RAW"

# 2. Find all the sample files by looking for "GSM"
archive_files <- list.files(data_dir, pattern = "GSM", full.names = TRUE)

# 3. Extract each one into its own dedicated folder
for (file in archive_files) {
  # Strip away all extensions (.tar.gz, .tar, etc.) for a clean folder name
  folder_name <- sub("\\..*", "", basename(file))
  out_path <- file.path(data_dir, folder_name)
  
  # Create the folder and extract
  dir.create(out_path, showWarnings = FALSE)
  untar(tarfile = file, exdir = out_path)
  
  message(paste("Extracted:", folder_name))
}

library(Seurat)
library(dplyr)
data_dir <- "E:/Viral_TLR4_Echo_Project/data/GSE140231_RAW"
sample_folders <- list.dirs(data_dir, recursive = FALSE)
seurat_list <- list()

for (folder in sample_folders) {
  sample_name <- basename(folder)
  files <- list.files(folder, full.names = TRUE)
  for (f in files) {
    if (grepl("barcodes|features|genes|matrix", basename(f))) {
      clean_name <- sub(".*_(barcodes|features|genes|matrix)", "\\1", basename(f))
      file.rename(f, file.path(folder, clean_name))
    }
  }
 
  for (folder in sample_folders) {
    sample_name <- basename(folder)
    files <- list.files(folder, full.names = TRUE)
    for (f in files) {
      if (grepl("barcodes|features|genes|matrix", basename(f))) {
        clean_name <- sub(".*_(barcodes|features|genes|matrix)", "\\1", basename(f))
        file.rename(f, file.path(folder, clean_name))
      }
    }
    data_dir <- "E:/Viral_TLR4_Echo_Project/data/GSE140231_RAW"
    sample_folders <- list.dirs(data_dir, recursive = FALSE)
    seurat_list <- list()    
    for (folder in sample_folders) {
      sample_name <- basename(folder)
      files <- list.files(folder, full.names = TRUE)
      for (f in files) {
        if (grepl("barcodes|features|genes|matrix", basename(f))) {
          clean_name <- sub(".*_(barcodes|features|genes|matrix)", "\\1", basename(f))
          file.rename(f, file.path(folder, clean_name))
        }
      }
      counts <- Read10X(data.dir = folder)
      seurat_obj <- CreateSeuratObject(counts = counts, project = sample_name, min.cells = 3, min.features = 200)
      seurat_list[[sample_name]] <- seurat_obj
      message(paste("Loaded:", sample_name))
    }
    
    library(Seurat)
    library(dplyr)
    
    data_dir <- "E:/Viral_TLR4_Echo_Project/data/GSE140231_RAW"
    sample_folders <- list.dirs(data_dir, recursive = FALSE)
    seurat_list <- list()
    
    for (folder in sample_folders) {
      sample_name <- basename(folder)
      
      # Step A: Brute-force rename the files inside the folder
      files <- list.files(folder, full.names = TRUE)
      for (f in files) {
        base_f <- basename(f)
        # Skip if already perfectly named
        if (base_f %in% c("barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz", "barcodes.tsv", "features.tsv", "matrix.mtx")) next
        
        # Rename based on keyword detection
        if (grepl("barcode", base_f, ignore.case = TRUE)) {
          new_name <- ifelse(grepl("\\.gz$", base_f), "barcodes.tsv.gz", "barcodes.tsv")
          file.rename(f, file.path(folder, new_name))
        } else if (grepl("feature|gene", base_f, ignore.case = TRUE)) {
          new_name <- ifelse(grepl("\\.gz$", base_f), "features.tsv.gz", "features.tsv")
          file.rename(f, file.path(folder, new_name))
        } else if (grepl("matrix", base_f, ignore.case = TRUE)) {
          new_name <- ifelse(grepl("\\.gz$", base_f), "matrix.mtx.gz", "matrix.mtx")
          file.rename(f, file.path(folder, new_name))
        }
      }
      
      # Step B: Load into Seurat
      counts <- Read10X(data.dir = folder)
      seurat_obj <- CreateSeuratObject(counts = counts, project = sample_name, min.cells = 3, min.features = 200)
      seurat_list[[sample_name]] <- seurat_obj
      message(paste("Loaded:", sample_name))
    }
    
    # Step C: Merge all 12 samples into one master dataset
    pd_merged <- merge(seurat_list[[1]], y = seurat_list[-1], add.cell.ids = names(seurat_list))
    rm(seurat_list, counts, seurat_obj)
    gc() # Clear RAM
    
    # Step D: Quality Control
    pd_merged[["percent.mt"]] <- PercentageFeatureSet(pd_merged, pattern = "^MT-")
    pd_filtered <- subset(pd_merged, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 5)
    
    message("Integration Complete. Final dataset:")
    print(pd_filtered)
    
library(Seurat)
    data_dir <- "E:/Viral_TLR4_Echo_Project/data/GSE140231_RAW"
    sample_folders <- list.dirs(data_dir, recursive = FALSE)
    seurat_list <- list()
    
    for (folder in sample_folders) {
      sample_name <- basename(folder)
      
      # Load the cleanly named, uncompressed files into Seurat
      counts <- Read10X(data.dir = folder)
      seurat_obj <- CreateSeuratObject(counts = counts, project = sample_name, min.cells = 3, min.features = 200)
      seurat_list[[sample_name]] <- seurat_obj
      message(paste("Loaded:", sample_name))
    }    
    
    library(Seurat)
    data_dir <- "E:/Viral_TLR4_Echo_Project/data/GSE140231_RAW"
    sample_folders <- list.dirs(data_dir, recursive = FALSE)
    seurat_list <- list()
    
    for (folder in sample_folders) {
      sample_name <- basename(folder) 
      # Step 1: Bulletproof Renaming (Bypasses Windows hidden extensions)
      actual_files <- list.files(folder, full.names = TRUE)
      for (f in actual_files) {
        if (grepl("barcode", basename(f), ignore.case = TRUE)) {
          file.rename(f, file.path(folder, "barcodes.tsv"))
        } else if (grepl("feature|gene", basename(f), ignore.case = TRUE)) {
          file.rename(f, file.path(folder, "features.tsv"))
        } else if (grepl("matrix", basename(f), ignore.case = TRUE)) {
          file.rename(f, file.path(folder, "matrix.mtx"))
        }
      }
      
      # Step 2: Load into Seurat
      counts <- Read10X(data.dir = folder)
      seurat_obj <- CreateSeuratObject(counts = counts, project = sample_name, min.cells = 3, min.features = 200)
      seurat_list[[sample_name]] <- seurat_obj
      message(paste("Loaded:", sample_name))
    }
    
    library(Seurat)
    
    data_dir <- "E:/Viral_TLR4_Echo_Project/data/GSE140231_RAW"
    sample_folders <- list.dirs(data_dir, recursive = FALSE)
    seurat_list <- list()
    
    for (folder in sample_folders) {
      sample_name <- basename(folder)
      
      # Step 1: Bulletproof Renaming (Bypasses Windows hidden extensions)
      actual_files <- list.files(folder, full.names = TRUE)
      for (f in actual_files) {
        if (grepl("barcode", basename(f), ignore.case = TRUE)) {
          file.rename(f, file.path(folder, "barcodes.tsv"))
        } else if (grepl("feature|gene", basename(f), ignore.case = TRUE)) {
          file.rename(f, file.path(folder, "features.tsv"))
        } else if (grepl("matrix", basename(f), ignore.case = TRUE)) {
          file.rename(f, file.path(folder, "matrix.mtx"))
        }
      }
      
      # Step 2: Load into Seurat
      counts <- Read10X(data.dir = folder)
      seurat_obj <- CreateSeuratObject(counts = counts, project = sample_name, min.cells = 3, min.features = 200)
      seurat_list[[sample_name]] <- seurat_obj
      message(paste("Loaded:", sample_name))
    }
    
    # Step 3: Merge all 12 samples into one master dataset
    pd_merged <- merge(seurat_list[[1]], y = seurat_list[-1], add.cell.ids = names(seurat_list))
    rm(seurat_list, counts, seurat_obj)
    gc()
    
    # Step 4: Quality Control
    pd_merged[["percent.mt"]] <- PercentageFeatureSet(pd_merged, pattern = "^MT-")
    pd_filtered <- subset(pd_merged, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 5)
    
    message("Integration Complete. Final dataset:")
    print(pd_filtered)
    
    library(Seurat)
    library(Matrix)
    
    data_dir <- "E:/Viral_TLR4_Echo_Project/data/GSE140231_RAW"
    sample_folders <- list.dirs(data_dir, recursive = FALSE)
    seurat_list <- list()
    
    for (folder in sample_folders) {
      sample_name <- basename(folder)
      
      # 1. Dynamically grab the exact paths of the 3 files, regardless of what Windows named them
      all_files <- list.files(folder, full.names = TRUE)
      mtx_file <- all_files[grepl("matrix", all_files, ignore.case = TRUE)][1]
      barcode_file <- all_files[grepl("barcode", all_files, ignore.case = TRUE)][1]
      feature_file <- all_files[grepl("feature|gene", all_files, ignore.case = TRUE)][1]
      
      # 2. Build the matrix manually (Bypassing Read10X completely)
      counts <- readMM(mtx_file)
      barcodes <- read.delim(barcode_file, header = FALSE)[,1]
      features <- read.delim(feature_file, header = FALSE)
      
      # Safely handle the genes whether the original authors used 1 column or 2 columns
      rownames(counts) <- make.unique(as.character(if(ncol(features) >= 2) features[,2] else features[,1]))
      colnames(counts) <- barcodes
      
      # 3. Create the Seurat Object
      seurat_obj <- CreateSeuratObject(counts = counts, project = sample_name, min.cells = 3, min.features = 200)
      seurat_list[[sample_name]] <- seurat_obj
      message(paste("Loaded:", sample_name))
    }
    
    # Step C: Merge all 12 samples into one master dataset
    pd_merged <- merge(seurat_list[[1]], y = seurat_list[-1], add.cell.ids = names(seurat_list))
    rm(seurat_list, counts, seurat_obj)
    gc() # Clear RAM
    
    # Step D: Quality Control
    pd_merged[["percent.mt"]] <- PercentageFeatureSet(pd_merged, pattern = "^MT-")
    pd_filtered <- subset(pd_merged, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 5)
    
    message("Integration Complete. Final dataset:")
    print(pd_filtered)
    
    library(Seurat)
    
    data_dir <- "E:/Viral_TLR4_Echo_Project/data/GSE140231_RAW"
    
    # THE FIX: list.dirs includes the parent folder. 
    # We use grepl to filter the list so it ONLY keeps the 12 actual sample folders.
    all_dirs <- list.dirs(data_dir, recursive = FALSE)
    sample_folders <- all_dirs[grepl("GSM", all_dirs)] 
    
    seurat_list <- list()
    
    for (folder in sample_folders) {
      sample_name <- basename(folder)
      
      # Safely load the data into Seurat
      counts <- Read10X(data.dir = folder)
      seurat_obj <- CreateSeuratObject(counts = counts, project = sample_name, min.cells = 3, min.features = 200)
      seurat_list[[sample_name]] <- seurat_obj
      message(paste("Successfully Loaded:", sample_name))
    }
    
    # Merge all 12 samples into one master human brain atlas
    pd_merged <- merge(seurat_list[[1]], y = seurat_list[-1], add.cell.ids = names(seurat_list))
    
    # Clear out the temporary individual files to free up your computer's RAM
    rm(seurat_list, counts, seurat_obj)
    gc() 
    
    # Standard Quality Control: Remove dying cells (high mitochondrial content) and doublets
    pd_merged[["percent.mt"]] <- PercentageFeatureSet(pd_merged, pattern = "^MT-")
    pd_filtered <- subset(pd_merged, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 5)
    
    message("Integration Complete. Final human PD dataset:")
    print(pd_filtered)
    
    
    library(Seurat)
    library(Matrix) # Required to read the raw matrix directly
    
    data_dir <- "E:/Viral_TLR4_Echo_Project/data/GSE140231_RAW"
    
    # 1. Filter out the parent directory so we ONLY loop through the 12 patient folders
    all_dirs <- list.dirs(data_dir, recursive = FALSE)
    sample_folders <- all_dirs[grepl("GSM", all_dirs)] 
    
    seurat_list <- list()
    
    for (folder in sample_folders) {
      sample_name <- basename(folder)
      
      # 2. Dynamically grab the exact paths of the 3 files, ignoring whatever Windows named them
      all_files <- list.files(folder, full.names = TRUE)
      mtx_file <- all_files[grepl("matrix", all_files, ignore.case = TRUE)][1]
      barcode_file <- all_files[grepl("barcode", all_files, ignore.case = TRUE)][1]
      feature_file <- all_files[grepl("feature|gene", all_files, ignore.case = TRUE)][1]
      
      # 3. Build the matrix manually (Bypassing Read10X completely)
      counts <- readMM(mtx_file)
      barcodes <- read.delim(barcode_file, header = FALSE)[,1]
      features <- read.delim(feature_file, header = FALSE)
      
      # Safely handle the genes whether the original authors used 1 column or 2 columns
      rownames(counts) <- make.unique(as.character(if(ncol(features) >= 2) features[,2] else features[,1]))
      colnames(counts) <- barcodes
      
      # 4. Create the Seurat Object
      seurat_obj <- CreateSeuratObject(counts = counts, project = sample_name, min.cells = 3, min.features = 200)
      seurat_list[[sample_name]] <- seurat_obj
      message(paste("Successfully Loaded:", sample_name))
    }
    
    # 5. Merge all 12 samples into one master dataset
    pd_merged <- merge(seurat_list[[1]], y = seurat_list[-1], add.cell.ids = names(seurat_list))
    rm(seurat_list, counts, seurat_obj)
    gc() # Clear RAM
    
    # 6. Quality Control
    pd_merged[["percent.mt"]] <- PercentageFeatureSet(pd_merged, pattern = "^MT-")
    pd_filtered <- subset(pd_merged, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 5)
    
    message("Integration Complete. Final dataset:")
    print(pd_filtered)    
    
    
    
    
    library(Seurat)
    library(Matrix) 
    
    data_dir <- "E:/Viral_TLR4_Echo_Project/data/GSE140231_RAW"
    all_dirs <- list.dirs(data_dir, recursive = FALSE)
    sample_folders <- all_dirs[grepl("GSM", all_dirs)] 
    
    seurat_list <- list()
    
    for (folder in sample_folders) {
      sample_name <- basename(folder)
      
      # THE FIX: 'recursive = TRUE' forces R to dig through the hidden Windows sub-folders
      mtx_file <- list.files(folder, pattern = "matrix", ignore.case = TRUE, recursive = TRUE, full.names = TRUE)[1]
      barcode_file <- list.files(folder, pattern = "barcode", ignore.case = TRUE, recursive = TRUE, full.names = TRUE)[1]
      feature_file <- list.files(folder, pattern = "feature|gene", ignore.case = TRUE, recursive = TRUE, full.names = TRUE)[1]
      
      # Build the matrix manually (immune to renaming errors)
      counts <- readMM(mtx_file)
      barcodes <- read.delim(barcode_file, header = FALSE)[,1]
      features <- read.delim(feature_file, header = FALSE)
      
      # Bind rownames and colnames automatically
      rownames(counts) <- make.unique(as.character(if(ncol(features) >= 2) features[,2] else features[,1]))
      colnames(counts) <- barcodes
      
      # Create Seurat Object
      seurat_obj <- CreateSeuratObject(counts = counts, project = sample_name, min.cells = 3, min.features = 200)
      seurat_list[[sample_name]] <- seurat_obj
      message(paste("Successfully Loaded:", sample_name))
    }
    
    # Merge all 12 samples into one master dataset
    pd_merged <- merge(seurat_list[[1]], y = seurat_list[-1], add.cell.ids = names(seurat_list))
    rm(seurat_list, counts, seurat_obj)
    gc() 
    
    # Quality Control
    pd_merged[["percent.mt"]] <- PercentageFeatureSet(pd_merged, pattern = "^MT-")
    pd_filtered <- subset(pd_merged, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 5)
    
    message("Integration Complete. Final human PD dataset:")
    print(pd_filtered)