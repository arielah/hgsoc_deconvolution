suppressPackageStartupMessages({
  library(data.table)
  library(SingleCellExperiment)
  library(dplyr)
  library(cluster)
  library(BayesPrism)
  library(yaml)
  library(ggplot2)
  library(stringr)
  library(pheatmap)
  library(gplots)
})

params <- read_yaml("../../config.yml")
data_path <- params$data_path
local_data_path <- params$local_data_path
plot_path <- params$plot_path

# Current options are AACES, TCGA (for RNA-seq), microarray, and tothill
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 1) {
  print(args)
  dataset <- args[1]
}

# Load bayesprism object
bp <- readRDS(paste(local_data_path, "/deconvolution_output/", dataset, 
                    "_default_bayesprism_results_full.rds", sep = ""))

# Get expression from only epithelial cell fraction
epithelium <- get.exp(bp, state.or.type = "type",
                      cell.name = "Epithelial cells")

# Get full expression data
full <- fread(paste(local_data_path, "/deconvolution_input/bulk_data_", dataset,
                    ".tsv", sep=""))
genes <- full$Gene; full$Gene <- NULL
full <- t(full)
colnames(full) <- genes

if (dataset == "TCGA") {
  rownames(epithelium) <- str_extract(rownames(epithelium), "TCGA-\\w\\w-\\w\\w\\w\\w")
  rownames(epithelium) <- gsub("-","\\.", rownames(epithelium))
}

# Get subtype annotations
cluster_file <- paste(local_data_path, "cluster_assignments", "FullClusterMembership.csv", sep = "/")
cluster_list <- fread(cluster_file)

# Remove genes with low expression
cs <- colSums(epithelium)
some_exp <- which(cs > 200)
epithelium <- epithelium[, some_exp]
full <- subset(full, select=colnames(epithelium))

# Kmeans clustering
set.seed(201)
km <- kmeans(epithelium, centers = 4)

kmeans <- as.data.frame(km$cluster)
kmeans$Sample <- names(km$cluster)
setnames(kmeans, c("epithelial_kmeans", "Sample"))

# Join epithelial kmeans with "full" kmeans
setnames(cluster_list, "V1", "Sample")
cluster_list <- subset(cluster_list, select=c("Sample", "ClusterK4_kmeans"))
kmeans <- inner_join(cluster_list, kmeans)

# Make heatmap
grouped_data <- kmeans %>%
  group_by(ClusterK4_kmeans, epithelial_kmeans) %>%
  summarize(count = n()) %>% as.data.frame()

# Reshape and then melt to replace NAs with 0s
confusion <- reshape(data = grouped_data,
                     idvar = "ClusterK4_kmeans",
                     v.names = "count",
                     timevar = "epithelial_kmeans",
                     direction = "wide")
colnames(confusion) <- gsub("count.", "", colnames(confusion))
confusion[is.na(confusion)] <- 0
grouped_data <- melt(confusion, id.vars = "ClusterK4_kmeans")
setnames(grouped_data, "variable", "epithelial_kmeans")

plotfile <- paste(plot_path, "/evaluation_plots/", dataset, "_kmeans_heatmap.png", sep = "")
png(plotfile)
ggplot(grouped_data, aes(x = ClusterK4_kmeans, y = epithelial_kmeans, fill = value)) +
  geom_raster() + geom_text(aes(label = value), size = 6) +
  xlab("Full data (Way pipeline) kmeans cluster") + ylab("Cancer fraction kmeans cluster")
dev.off()

# Convert expression to log with pseudocounts
epithelium <- epithelium + 1
epithelium <- log(epithelium)
epithelium <- subset(epithelium, rownames(epithelium) %in% kmeans$Sample)

full <- full + 1
full <- log(full)

annotation_row = data.frame(
    Full_kmeans = factor(kmeans$ClusterK4_kmeans),
    Epi_kmeans = factor(kmeans$epithelial_kmeans)
)
rownames(annotation_row) = kmeans$Sample

plotname <- paste(local_data_path, "/evaluation_plots/",
                  dataset, "_expression_heatmap.png")
png(plotname)
pheatmap(epithelium, annotation_row = annotation_row,
         show_colnames = F, show_rownames = F, scale = "column")
dev.off()

