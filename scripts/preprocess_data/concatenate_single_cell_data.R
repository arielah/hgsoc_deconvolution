suppressPackageStartupMessages({
    library(SingleCellExperiment)
    library(dplyr)
    library(yaml)
})

params <- read_yaml("../../config.yml")
data_path <- params$data_path
local_data_path <- params$local_data_path
sc_samples <- params$sc_samples

# Combine into a mega-object
for (i in 1:length(sc_samples)){
    scepath <- paste(local_data_path, "sce_objects/", sep = "/")
    tmp_sce <- readRDS(paste(scepath, sc_samples[i], "_labeled.rds", sep = ""))
 
    if (i == 1) {
        sce <- tmp_sce
    } else {
        sce <- cbind(sce, tmp_sce)
    }
}

# BayesPrism expects a "cell state" for each cell. For most cell types, this is
# identical to their cell type label. But for cancer (epithelial) cells, cell
# state is the sample of origin to account for inter-patient heterogeneity. The
# pooled samples 12162021 and 01132022 were already labeled by sample of origin
# in the annotation step, so they aren't relabeled here.
sce$cellTypeGranular <- ifelse(sce$cellType != "Epithelial cells" |
                                   sce$Sample %in% c("12162021", "01132022"),
                               sce$cellTypeGranular, sce$Sample)

# Get rid of some metadata for space reasons
colData(sce) <- subset(colData(sce), select = c("Sample", "Barcode", "clusters",
                                              "cellType", "cellTypeGranular"))

# Save as default single cell file
outfile <- paste(local_data_path, "deconvolution_input",
                 "single_cell_data_default.rds", sep = "/")
saveRDS(sce, outfile)

# I'm not confident in my fibroblast vs. smooth muscle cell annotations,
# so I'll try running them as one category with the distinction in Cell State
sce[, sce$cellType == "Smooth muscle cells"]$cellType <- "Fibroblasts"

outfile <- paste(local_data_path, "deconvolution_input",
                 "single_cell_data_fibro_merged.rds", sep = "/")
saveRDS(sce, outfile)
