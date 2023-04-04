# BayesPrism has a built in function to run NMF and try and find distinct transcriptional
# programs within the cancer fraction, once the stromal fraction has been computationally
# removed. To do this as with any NMF, you need to select a value of k for the number of
# latent embeddings. This script generates plots to assess the dispersion of the data and
# is based off the BayesPrism tutorial:
# https://github.com/Danko-Lab/BayesPrism/blob/main/tutorial_embedding_learning.html

suppressPackageStartupMessages({
  library(data.table)
  library(SingleCellExperiment)
  library(dplyr)
  library(curatedOvarianData)
  library(yaml)
  library(stringr)
  library(ggplot2)
  library(ggrepel)
  library(BayesPrism)
})

params <- read_yaml("../../config.yml")
data_path <- params$data_path
local_data_path <- params$local_data_path
plot_path <- params$plot_path

# Current options are TCGA (RNA-seq), microarray, and tothill
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 1) {
  print(args)
  dataset <- args[1]
}

# Load bayesprism object
bp <- readRDS(paste(local_data_path, "/deconvolution_output/", dataset,
                    "_default_bayesprism_results_full.rds", sep = ""))

# Load malignant gene expression
Z.tum.norm <- t(bp@reference.update@psi_mal)

# Scan a range of K values for number of malignant programs
estim.Z.tum.norm <- nmf(Z.tum.norm, rank = 2:12, seed = 123456)

filename <- paste(local_data_path, "/embeddings/", dataset, "_malignant_nmf_rank.rds", sep = "")
saveRDS(estim.Z.tum.norm, filename)

png(paste(plot_path, "/embedding_plots/", dataset, "_nmf_rank.png", sep = ""), width = 800)
plot(estim.Z.tum.norm)
dev.off()

png(paste(plot_path, "/embedding_plots/", dataset, "_consensus_map.png", sep = ""), width = 800)
consensusmap(estim.Z.tum.norm, labCol = NA, labRow = NA)
dev.off()
