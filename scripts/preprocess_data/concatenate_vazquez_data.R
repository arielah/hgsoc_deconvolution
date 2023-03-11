# The Vazquez-Garcia data (GSE180661) has over a million cells and is too big
# to fit into memory all at once, so I've broken the mtx file I made out of the
# HDF5 file into chunks. I'm guessing a million cells is far too much for
# BayesPrism either, so we're going to start by taking a random ~100,000 cell
# sample and running deconvolution with that. This file loads all of the chunks
# sequentially and filters them down to only the chosen cells.

suppressPackageStartupMessages({
  library(DropletUtils)
  library(SingleCellExperiment)
  library(data.table)
  library(scater)
  library(dplyr)
  library(yaml)
})

params <- read_yaml("../../config.yml")
data_path <- params$data_path
local_data_path <- params$local_data_path
plot_path <- params$plot_path

# Get list of cells with cell type labels
cells_with_cell_types <- fread(paste(local_data_path, "Vazquez-Garcia", 
                          "GSE180661_GEO_cells.tsv.gz", sep = "/"))

# Take random sample of cells
set.seed(302)
sample_indices <- sample(1:nrow(cells_with_cell_types), size = 100000)
kept_cells <- cells_with_cell_types[sample_indices, ]

# Iterate through all file chunks
parts <- list.files(path = paste(local_data_path, "Vazquez-Garcia", sep = "/"), pattern = "part_*")
for(i in 1:length(parts)){
  chunk_dir <- paste(local_data_path, "Vazquez-Garcia",
                     parts[i], sep = "/")
  tmp_sce <- read10xCounts(chunk_dir); gc()
  
  # Since we had to put the full matrix size for each chunk, each sce object
  # will be the full size of the dataset but with mostly blank cells. One way
  # to remove these cells is to calculate total number of reads in each cell
  # and throw out the cells with 0 reads. We'll also filter down to the cells
  # in our sample here as well.
  tmp_sce <- addPerCellQC(tmp_sce)
  tmp_sce <- tmp_sce[, tmp_sce$Barcode %in% kept_cells$cell_id &
                       tmp_sce$sum > 0]; gc()
  tmp_kept_cells <- kept_cells[kept_cells$cell_id %in% tmp_sce$Barcode, ]
  
  # Add cell type info to new cells
  setnames(tmp_kept_cells, "cell_id", "Barcode")
  colData(tmp_sce) <- merge(colData(tmp_sce), tmp_kept_cells,
                            by = "Barcode", sort = FALSE)
  
  # Merge samples from file chunk into main object
  if(i == 1){
    sce <- tmp_sce
  } else{
    sce <- cbind(sce, tmp_sce)
  }
  gc()
}

# Make UMAP to check that mtx files formatted properly
sce <- logNormCounts(sce); gc()
set.seed(317)
sce <- runUMAP(sce,
               BNPARAM = BiocNeighbors::AnnoyParam(),
               BPPARAM = BiocParallel::MulticoreParam(),
               ## unnecessary options,  only used to make a pretty graph
               min_dist = 0.5,  repulsion_strength = 0.25,
               spread = 0.7,
               n_neighbors = 15)

plotUMAP(sce, colour_by="cell_type_super")

# Rename cell types to match conventions of existing analysis
sce$cellType <- recode(sce$cell_type_super,
                       "B.super" = "B cells",
                       "Endothelial.super" = "Endothelial cells",
                       "Fibroblast.super" = "Fibroblasts",
                       "Myeloid.super" = "Myeloid cells",
                       "Ovarian.cancer.super" = "Epithelial cells",
                       "T.super" = "T cells")

# BayesPrism expects a "cell state" for each cell. For most cell types, this is
# identical to their cell type label. But for cancer (epithelial) cells, cell
# state is the sample of origin to account for inter-patient heterogeneity.
sce$cellTypeGranular <- ifelse(sce$cellType != "Epithelial cells",
                               sce$cell_type, sce$sample)

# Get rid of some metadata for space reasons
colData(sce) <- subset(colData(sce), select = c("Sample", "Barcode",
                                                "cellType", "cellTypeGranular"))

# Save as Vazquez single cell file
outfile <- paste(local_data_path, "deconvolution_input",
                 "single_cell_data_vazquez.rds", sep = "/")
saveRDS(sce, outfile)
