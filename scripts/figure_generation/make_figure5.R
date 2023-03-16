suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(stringr)
  library(pheatmap)
  library(dplyr)
  library(yaml)
})

source("figure_utils.R")

params <- read_yaml("../../config.yml")
data_path <- params$data_path
local_data_path <- params$local_data_path
plot_path <- params$plot_path
figure_path <- params$figure_path

# Load embeddings data
load_data <- function(bulk_set, k) {
  embedding_file <- paste(local_data_path, "/embeddings/", bulk_set,
                          "_k", k, "_embeddings.rds", sep = "")
  ebd.res <- readRDS(embedding_file)
  
  ebd.res
}

tcga <- load_data("TCGA", 3)
microarray <- load_data("microarray", 3)
tothill <- load_data("tothill", 3)

# Load subtype assignments
cluster_file <- paste(local_data_path, "cluster_assignments", "noAACES_FullClusterMembership.csv", sep = "/")
cluster_list <- fread(cluster_file)
cluster_list <- cluster_list[order(cluster_list$ClusterK4_kmeans),]

cluster_list$ClusterK4_kmeans <- recode(cluster_list$ClusterK4_kmeans,
                                        "1" = "Mesenchymal",
                                        "2" = "Proliferative",
                                        "3" = "Immunoreactive",
                                        "4" = "Differentiated")

plot_weights <- function(ebd.res, name) {
  weights <- ebd.res$omega
  
  if(name == "TCGA RNA-seq") {
    rownames(weights) <- str_extract(rownames(weights), "TCGA-\\w\\w-\\w\\w\\w\\w")
    rownames(weights) <- gsub("-","\\.",rownames(weights))
  }
  
  clusters <- subset(cluster_list, cluster_list$V1 %in% rownames(weights))
  weights <- weights[rownames(weights) %in% cluster_list$V1, ]
  weights <- weights[match(clusters$V1, rownames(weights)),]
  
  row_ha <- rowAnnotation(
    k4_kmeans = clusters$ClusterK4_kmeans,
    k3_kmeans = factor(clusters$ClusterK3_kmeans),
    col = list(k3_kmeans = colors_subtypes[5:7],
               k4_kmeans = colors_subtypes[1:4])
  )
  
  Heatmap(weights, left_annotation = row_ha, show_row_names = F,
          cluster_rows = F, cluster_columns = F)
}

pA <- plot_weights(tcga, "TCGA RNA-seq")
pB <- plot_weights(microarray, "TCGA Microarray")
pC <- plot_weights(tothill, "Tothill")

