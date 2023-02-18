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

# Current options for dataset are AACES, TCGA (for RNA-seq), microarray, and tothill
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 2) {
  print(args)
  dataset <- args[1]
  k <- args[2]
} else{
  stop("Need to input a dataset and a number of latent factors k")
}

# Load bayesprism object
bp <- readRDS(paste(local_data_path, "/deconvolution_output/", dataset, 
                    "_default_bayesprism_results_full.rds", sep = ""))

# Run embeddings, this is the part that takes ~1 day
ebd.res <- learn.embedding.nmf(bp = bp, K = k, cycle = 40)

# Save data
embedding_file <- paste(local_data_path, "/embeddings/", dataset,
                        "_k", k, "_embeddings.rds", sep = "")
saveRDS(ebd.res, embedding_file)

weights <- ebd.res$omega
weights_file <- paste(local_data_path, "/embeddings/", dataset,
                      "_k", k, "_program_weights.tsv", sep = "")
write.table(weights, weights_file, sep = "\t", quote = F)

# Load clustering annotations
cluster_file <- paste(local_data_path, "cluster_assignments", "FullClusterMembership.csv", sep = "/")
cluster_list <- fread(cluster_file)

cluster_list$ClusterK4_kmeans <- recode(cluster_list$ClusterK4_kmeans,
                                        "1" = "Mesenchymal",
                                        "2" = "Proliferative",
                                        "3" = "Immunoreactive",
                                        "4" = "Differentiated")
cluster_list <- subset(cluster_list, cluster_list$V1 %in% rownames(weights))

annotation_row = data.frame(
  Full_kmeans = factor(cluster_list$ClusterK4_kmeans)
)
rownames(annotation_row) = cluster_list$V1

# Plot weights
plotfile <- paste(plot_path, "/embedding_plots/", dataset,
                  "_k", k, "_embedding_heatmap.png", sep = "")
png(plotfile)
pheatmap(weights, annotation_row = annotation_row, show_rownames = F)
dev.off()