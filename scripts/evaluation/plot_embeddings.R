# Script to plot data from generate_embeddings.R. We'll try it a couple of
# different ways, using pheatmap's default clustering and also ordering based
# on subtype cluster assignmentsu. We'll also use GSVA to compare expression
# programs to existing (hallmark) gene sets.

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
  library(GSVA)
  library(GSEABase)
  library(pheatmap)
})


params <- read_yaml("../../config.yml")
data_path <- params$data_path
local_data_path <- params$local_data_path
plot_path <- params$plot_path

# Current options for dataset are TCGA (for RNA-seq), microarray, and tothill
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 2) {
  print(args)
  dataset <- args[1]
  k <- as.integer(args[2])
} else {
  stop("Need to input a dataset and a number of latent factors k")
}

# Load data
embedding_file <- paste(local_data_path, "/embeddings/", dataset,
                        "_k", k, "_embeddings.rds", sep = "")
ebd.res <- readRDS(embedding_file)

weights <- ebd.res$omega

# Load clustering annotations
cluster_file <- paste(local_data_path, "cluster_assignments", "FullClusterMembership.csv", sep = "/")
cluster_list <- fread(cluster_file)

cluster_list <- cluster_list[order(cluster_list$ClusterK4_kmeans), ]

cluster_list$ClusterK4_kmeans <- recode(cluster_list$ClusterK4_kmeans,
                                        "1" = "Mesenchymal",
                                        "2" = "Proliferative",
                                        "3" = "Immunoreactive",
                                        "4" = "Differentiated")

if (dataset == "TCGA") {
  rownames(weights) <- str_extract(rownames(weights), "TCGA-\\w\\w-\\w\\w\\w\\w")
  rownames(weights) <- gsub("-", "\\.", rownames(weights))
}

cluster_list <- subset(cluster_list, cluster_list$V1 %in% rownames(weights))
weights <- weights[rownames(weights) %in% cluster_list$V1, ]
weights <- weights[match(cluster_list$V1, rownames(weights)), ]

annotation_row <- data.frame(
  Full_kmeans = factor(cluster_list$ClusterK4_kmeans)
)
rownames(annotation_row) <- cluster_list$V1

# Plot weights
plotfile <- paste(plot_path, "/embedding_plots/", dataset,
                  "_k", k, "_embedding_heatmap.png", sep = "")
png(plotfile)
pheatmap(weights, annotation_row = annotation_row, show_rownames = FALSE,
         cluster_cols = FALSE)
dev.off()

# Plot weights with no clustering
plotfile <- paste(plot_path, "/embedding_plots/", dataset,
                  "_k", k, "_embedding_heatmap_no_clustering.png", sep = "")
png(plotfile)
pheatmap(weights, annotation_row = annotation_row, show_rownames = FALSE,
         cluster_rows = FALSE, cluster_cols = FALSE)
dev.off()

# Plot weights with k3 labels
annotation_row <- data.frame(
  Full_kmeans = factor(cluster_list$ClusterK4_kmeans),
  K3_kmeans = factor(cluster_list$ClusterK3_kmeans)
)
rownames(annotation_row) <- cluster_list$V1

plotfile <- paste(plot_path, "/embedding_plots/", dataset,
                  "_k", k, "_embedding_heatmap_multilabel.png", sep = "")
png(plotfile)
pheatmap(weights, annotation_row = annotation_row, show_rownames = FALSE,
         cluster_cols = FALSE)
dev.off()

gene_expression <- t(ebd.res$eta_post)
c8 <- getGmt("~/Documents/scRNA/deconvolution_pilot/scripts/bulk_de/GSEA_custom_sets/c8.all.v7.5.1.symbols.gmt")
h <- getGmt("~/Downloads/h.all.v2023.1.Hs.symbols.gmt")

gsva.es <- gsva(gene_expression, h)

plotfile <- paste(plot_path, "/embedding_plots/", dataset,
                  "_k", k, "_gsva_hallmark_results.png", sep = "")
png(plotfile)
pheatmap(gsva.es, cluster_cols = FALSE)
dev.off()
