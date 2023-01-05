suppressPackageStartupMessages({
  library(data.table)
  library(scater)
  library(scran)
  library(dplyr)
  library(ggplot2)
  library(yaml)
})

params <- read_yaml("../../config.yml")
data_path <- params$data_path
local_data_path <- params$local_data_path
samples <- params$samples

sce <- readRDS(paste(local_data_path,"sce_objects/16030X2_clustered.rds",sep="/"))

plotUMAP(sce, colour_by="clusters")

ct <- fread(paste(local_data_path,"celltypist_output/16030X2_predicted_labels.csv", sep="/"))

label_table <- fread(paste(local_data_path,"celltypist_output/simplified_labels.tsv", sep="/"))
setnames(label_table, "Original", "majority_voting")
ct <- left_join(ct, label_table)
sce$CT <- ct$Simplified
sce$granular <- ct$majority_voting

# Switch NK cells to smooth muscle cells
ct[ct$Simplified=="NK cells",]$Simplified <- "Smooth muscle cells"
sce[, sce$CT=="NK cells"]$CT <- "Smooth muscle cells"

ct$keep <- FALSE
ct$clusters <- sce$clusters
ct[ct$clusters %in% c(1, 2, 3, 4, 6) & ct$Simplified %in% c("Epithelial cells"),]$keep <- TRUE
ct[ct$clusters == 5 & ct$Simplified %in% c("T cells", "Macrophages", "Monocytes", "DC"), ]$keep <- TRUE
ct[ct$clusters == 7 & ct$Simplified %in% c("Fibroblasts", "Smooth muscle cells", "Endothelial cells"),]$keep <- TRUE
ct[ct$clusters == 4 & ct$Simplified %in% c("Endothelial cells", "Plasma cells"), ]$keep <- TRUE

sce$keep <- ct$keep
sce <- sce[, sce$keep]
sce$keep <- FALSE

sce$cellType <- sce$CT
sce$cellTypeGranular <- ifelse(sce$cellType=="T cells", sce$granular, sce$cellType)
sce[, sce$cellTypeGranular == "Tem/Trm cytotoxic T cells"]$cellTypeGranular <- "CD8 T cells"

sce$CT <- NULL
sce$granular <- NULL

saveRDS(sce, paste(local_data_path, "sce_objects/16030X2_labeled.rds", sep = "/"))
