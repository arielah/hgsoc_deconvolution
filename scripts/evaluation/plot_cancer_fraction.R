suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(yaml)
  library(BayesPrism)
  library(umap)
  library(ggplot2)
  library(stringr)
})

params <- read_yaml("../../config.yml")
data_path <- params$data_path
local_data_path <- params$local_data_path
plot_path <- params$plot_path

# Current options are TCGA (for RNA-seq), microarray, and tothill
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 1) {
	print(args)
	dataset <- args[1]
}

print(dataset)

# Load bayesprism object
bp <- readRDS(paste(local_data_path, "/deconvolution_output/", dataset,
                       "_default_bayesprism_results_full.rds", sep = ""))

# Get expression from only epithelial cell fraction
epithelium <- get.exp(bp, state.or.type = "type",
                      cell.name = "Epithelial cells")

# Get subtype annotations
cluster_file <- paste(local_data_path, "cluster_assignments", "FullClusterMembership.csv", sep = "/")
cluster_list <- fread(cluster_file)

if (dataset == "TCGA") {
	cluster_list$V1 <- gsub("\\.", "-", cluster_list$V1)
}
setnames(cluster_list, "V1", "ID")

# Calculate UMAP
set.seed(109)
epithelium_umap <- umap(epithelium)
epithelium_umap <- as.data.frame(epithelium_umap[[1]])

epithelium_umap$ID <- rownames(epithelium_umap)
if (dataset == "TCGA") {
  epithelium_umap$ID <- str_extract(epithelium_umap$ID, "TCGA-\\w\\w-\\w\\w\\w\\w")
}

epithelium_umap <- inner_join(epithelium_umap, cluster_list)
epithelium_umap$ClusterK4_kmeans <- as.factor(epithelium_umap$ClusterK4_kmeans)
epithelium_umap$ClusterK3_kmeans <- as.factor(epithelium_umap$ClusterK3_kmeans)
epithelium_umap$ClusterK2_kmeans <- as.factor(epithelium_umap$ClusterK2_kmeans)

epithelium_umap$ClusterK4_kmeans <- recode(epithelium_umap$ClusterK4_kmeans,
                                           "1" = "Mesenchymal",
                                           "2" = "Proliferative",
                                           "3" = "Immunoreactive",
                                           "4" = "Differentiated")

# Plot UMAP
png(paste(plot_path, "/evaluation_plots/", dataset, "_cancer_fraction_K2.png", sep = ""))
ggplot(data = epithelium_umap, mapping = aes(x = V1, y = V2, color = ClusterK2_kmeans)) +
  geom_point() + ggtitle("Cancer Fraction Expression") + xlab("UMAP 1") + ylab("UMAP 2")
dev.off()
png(paste(plot_path, "/evaluation_plots/", dataset, "_cancer_fraction_K3.png", sep = ""))
ggplot(data = epithelium_umap, mapping = aes(x = V1, y = V2, color = ClusterK3_kmeans)) +
  geom_point() + ggtitle("Cancer Fraction Expression") + xlab("UMAP 1") + ylab("UMAP 2")
dev.off()
png(paste(plot_path, "/evaluation_plots/", dataset, "_cancer_fraction_K4.png", sep = ""))
ggplot(data = epithelium_umap, mapping = aes(x = V1, y = V2, color = ClusterK4_kmeans)) +
  geom_point() + ggtitle("Cancer Fraction Expression") + xlab("UMAP 1") + ylab("UMAP 2")
dev.off()

# Try with different random seed
set.seed(118)
epithelium_umap2 <- umap(epithelium)
epithelium_umap2 <- as.data.frame(epithelium_umap2[[1]])

setnames(epithelium_umap2, c("UMAP1", "UMAP2"))
epithelium_umap2$ID <- rownames(epithelium_umap2)
if (dataset == "TCGA") {
  epithelium_umap2$ID <- str_extract(epithelium_umap2$ID, "TCGA-\\w\\w-\\w\\w\\w\\w")
}
epithelium_umap <- inner_join(epithelium_umap, epithelium_umap2)

# Plot UMAP
png(paste(plot_path, "/evaluation_plots/", dataset, "_cancer_fraction_seed2_K2.png", sep = ""))
ggplot(data = epithelium_umap, mapping = aes(x = UMAP1, y = UMAP2, color = ClusterK2_kmeans)) +
  geom_point() + ggtitle("Cancer Fraction Expression") + xlab("UMAP 1") + ylab("UMAP 2")
dev.off()
png(paste(plot_path, "/evaluation_plots/", dataset, "_cancer_fraction_seed2_K3.png", sep = ""))
ggplot(data = epithelium_umap, mapping = aes(x = UMAP1, y = UMAP2, color = ClusterK3_kmeans)) +
  geom_point() + ggtitle("Cancer Fraction Expression") + xlab("UMAP 1") + ylab("UMAP 2")
dev.off()
png(paste(plot_path, "/evaluation_plots/", dataset, "_cancer_fraction_seed2_K4.png", sep = ""))
ggplot(data = epithelium_umap, mapping = aes(x = UMAP1, y = UMAP2, color = ClusterK4_kmeans)) +
  geom_point() + ggtitle("Cancer Fraction Expression") + xlab("UMAP 1") + ylab("UMAP 2")
dev.off()
