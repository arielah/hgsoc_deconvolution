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

sce <- readRDS(paste(local_data_path,"sce_objects/01132022_clustered.rds",sep="/"))
sce$Sample <- "01132022"

plotUMAP(sce, colour_by="clusters")

ct <- fread(paste(local_data_path,"celltypist_output/01132022_predicted_labels.csv", sep="/"))

label_table <- fread(paste(local_data_path,"celltypist_output/simplified_labels.tsv", sep="/"))
setnames(label_table, "Original", "majority_voting")
ct <- left_join(ct, label_table)
sce$CT <- ct$Simplified
sce$granular <- ct$majority_voting

plotUMAP(sce, colour_by="CT")
plotUMAP(sce, colour_by="WFDC2")
plotUMAP(sce, colour_by="LYZ")
plotUMAP(sce, colour_by="CD3D")
plotUMAP(sce, colour_by="VWF")
plotUMAP(sce, colour_by="NKG7")
plotUMAP(sce, colour_by="DCN")
plotUMAP(sce, colour_by="ACTA2")
plotUMAP(sce, colour_by="JCHAIN")
plotUMAP(sce, colour_by="MS4A1")

ct$keep <- FALSE
ct$clusters <- sce$clusters

# Switch DCN- fibroblasts/"NK cells" to smooth muscle
ct[ct$clusters==6,]$Simplified <- "Smooth muscle cells"
sce[,sce$clusters==6]$CT <- "Smooth muscle cells"



ct[ct$clusters %in% c(1) & ct$Simplified %in% c("Endothelial cells"),]$keep <- TRUE
ct[ct$clusters %in% c(2) & ct$Simplified %in% c("Epithelial cells"),]$keep <- TRUE
ct[ct$clusters %in% c(3) & ct$Simplified %in% c("Plasma cells"),]$keep <- TRUE
ct[ct$clusters %in% c(4, 5, 7, 10) & ct$Simplified %in% c("T cells", "NK cells",
                                                          "ILC", "pDC",
                                                          "Mast cells")]$keep <- TRUE
ct[ct$clusters %in% c(6) & ct$Simplified %in% c("Smooth muscle cells"),]$keep <- TRUE
ct[ct$clusters %in% c(8) & ct$Simplified %in% c("Monocytes", "Macrophages",
                                                "DC"),]$keep <- TRUE
ct[ct$clusters %in% c(9) & ct$Simplified %in% c("B cells"),]$keep <- TRUE
ct[ct$clusters %in% c(11) & ct$Simplified %in% c("Fibroblasts"),]$keep <- TRUE



sce$keep <- ct$keep
sce <- sce[, sce$keep]
sce$keep <- NULL

sce$cellType <- sce$CT
sce$cellTypeGranular <- ifelse(sce$cellType=="T cells", sce$granular, sce$cellType)
sce[, sce$cellTypeGranular == "Tem/Trm cytotoxic T cells"]$cellTypeGranular <- "CD8 T cells"
sce[, sce$cellTypeGranular == "Tcm/Naive helper T cells"]$cellTypeGranular <- "CD4 T cells"

sce$CT <- NULL
sce$granular <- NULL

# Load in cell labels from genetic demultiplexing, to remove unassigned cells
labels <- fread(paste(data_path,"pooled_tumors/01132022/vireo/chunk_ribo/donor_ids.tsv", sep = "/"))
labels$Sample <- "01132022"

# Swap donor ids for the January pool for merging with December samples
labels$donor_id <- recode(labels$donor_id,
                          "donor0" = "donor4",
                          "donor1" = "donor5",
                          "donor2" = "donor6",
                          "donor3" = "donor7")
setnames(labels, "cell", "Barcode")
labels <- subset(labels, select=c("Barcode", "donor_id"))

colData(sce) <- merge(colData(sce), labels, by = "Barcode", sort = FALSE)

sce <- sce[, sce$donor_id != "unassigned" & sce$donor_id != "doublet"]
sce$cellTypeGranular <- ifelse(sce$cellType=="Epithelial cells",
                               sce$donor_id, sce$cellTypeGranular)

sce$donor_id <- NULL

saveRDS(sce, paste(local_data_path, "sce_objects/01132022_labeled.rds", sep = "/"))
