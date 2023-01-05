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

sce <- readRDS(paste(local_data_path,"sce_objects/19833X2_clustered.rds",sep="/"))

plotUMAP(sce, colour_by="clusters")

ct <- fread(paste(local_data_path,"celltypist_output/19833X2_predicted_labels.csv", sep="/"))

label_table <- fread(paste(local_data_path,"celltypist_output/simplified_labels.tsv", sep="/"))
setnames(label_table, "Original", "majority_voting")
ct <- left_join(ct, label_table)
sce$CT <- ct$Simplified
sce$granular <- ct$majority_voting

plotUMAP(sce, colour_by="CT")
plotUMAP(sce, colour_by="WFDC2")
plotUMAP(sce, colour_by="JCHAIN")
plotUMAP(sce, colour_by="MS4A1")
plotUMAP(sce, colour_by="LYZ")
plotUMAP(sce, colour_by="CD3D")
plotUMAP(sce, colour_by="VWF")
plotUMAP(sce, colour_by="NKG7")
plotUMAP(sce, colour_by="DCN")
plotUMAP(sce, colour_by="ACTA2")

ct$keep <- FALSE
ct$clusters <- sce$clusters
ct[ct$clusters %in% c(1, 5) & ct$Simplified %in% c("Epithelial cells", "Fibroblasts",
                                                   "Endothelial cells"),]$keep <- TRUE
ct[ct$clusters %in% c(2) & ct$Simplified %in% c("Macrophages", "DC"),]$keep <- TRUE
ct[ct$clusters %in% c(3, 4, 7) & ct$Simplified %in% c("T cells", "NK cells"),]$keep <- TRUE
ct[ct$clusters %in% c(5) & ct$Simplified %in% c("pDC", "DC", "B cells",
                                                "Macrophages"),]$keep <- TRUE
ct[ct$clusters %in% c(6) & ct$Simplified %in% c("Plasma cells"), ]$keep <- TRUE
ct[ct$clusters %in% c(8) & ct$Simplified %in% c("ILC"), ]$keep <- TRUE

sce$keep <- ct$keep
sce <- sce[, sce$keep]
sce$keep <- NULL

sce$cellType <- sce$CT
sce$cellTypeGranular <- ifelse(sce$cellType=="T cells", sce$granular, sce$cellType)
sce[, sce$cellTypeGranular == "Tem/Trm cytotoxic T cells"]$cellTypeGranular <- "CD8 T cells"

sce$CT <- NULL
sce$granular <- NULL

saveRDS(sce, paste(local_data_path, "sce_objects/19833X2_labeled.rds", sep = "/"))
