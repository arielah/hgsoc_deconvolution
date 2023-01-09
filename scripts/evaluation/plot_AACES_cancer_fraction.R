suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(yaml)
  library(BayesPrism)
  library(umap)
  library(ggplot2)
})

params <- read_yaml("../../config.yml")
data_path <- params$data_path
local_data_path <- params$local_data_path
plot_path <- params$plot_path

# Load AACES object
aaces <- readRDS(paste(local_data_path, "deconvolution_output", 
                       "AACES_default_bayesprism_results_full.rds", sep = "/"))

# Get expression from only epithelial cell fraction
epithelium <- get.exp(aaces, state.or.type = "type",
                      cell.name = "Epithelial cells")

# Get subtype annotations
cluster_file <- paste(local_data_path, "cluster_assignments", "FullClusterMembership.csv", sep = "/")
cluster_list <- fread(cluster_file)
setnames(cluster_list, "V1", "ID")

# Calculate UMAP
set.seed(109)
epithelium_umap <- umap(epithelium)
epithelium_umap <- as.data.frame(epithelium_umap[[1]])

epithelium_umap$ID <- rownames(epithelium_umap)

epithelium_umap <- left_join(epithelium_umap, cluster_list)
epithelium_umap$ClusterK4_kmeans <- as.factor(epithelium_umap$ClusterK4_kmeans)
epithelium_umap$ClusterK2_kmeans <- as.factor(epithelium_umap$ClusterK2_kmeans)

# Plot UMAP
ggplot(data = epithelium_umap, mapping = aes(x = V1, y = V2, color = ClusterK2_kmeans)) + 
  geom_point() + ggtitle("AACES Cancer Fraction Expression") + xlab("UMAP 1") + ylab("UMAP 2")
ggplot(data = epithelium_umap, mapping = aes(x = V1, y = V2, color = ClusterK4_kmeans)) + 
  geom_point() + ggtitle("AACES Cancer Fraction Expression") + xlab("UMAP 1") + ylab("UMAP 2")