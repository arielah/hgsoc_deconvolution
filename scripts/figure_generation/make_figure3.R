suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(stringr)
  library(ggsankey)
  library(patchwork)
  library(dplyr)
  library(yaml)
})

source("figure_utils.R")

params <- read_yaml("../../config.yml")
data_path <- params$data_path
local_data_path <- params$local_data_path
plot_path <- params$plot_path
figure_path <- params$figure_path

load_datasets <- function() {
  sc_set <- "fibro"
  bulk_sets <- c("TCGA", "microarray", "tothill")
  
  results_all <- data.frame()
  
  for (i in 1:length(bulk_sets)){
    bulk_set <- bulk_sets[i]
    
    results <- fread(paste(local_data_path, "/deconvolution_output/",
                           bulk_set, "_", sc_set, "_bayesprism_results.tsv",
                           sep = ""))
    
    if(bulk_set == "TCGA") {
      colnames(results) <- str_extract(colnames(results), "TCGA-\\w\\w-\\w\\w\\w\\w")
      colnames(results) <- gsub("-","\\.", colnames(results))
    }
    
    cell_types <- results$cell_type
    results$cell_type <- NULL
    results_t <- t(as.matrix(results))
    colnames(results_t) <- cell_types
    results_t <- as.data.frame(results_t)
    results_t$Dataset <- bulk_set
    results_t$Sample <- rownames(results_t)
    
    results_all <- rbind(results_all, results_t)
  }
  
  results_all
}

everything <- load_datasets()

# Get subtype annotations
cluster_file <- paste(local_data_path, "cluster_assignments", "FullClusterMembership.csv", sep = "/")
cluster_list <- fread(cluster_file)
cluster_list$Dataset <- NULL
setnames(cluster_list, "V1", "Sample")

# Merge cell type estimates and subtype annotations
everything <- full_join(everything, cluster_list)
everything <- subset(everything, !is.na(everything$`Epithelial cells`) &
                            !is.na(everything$ClusterK4_kmeans))

everything$Subtype <- recode(everything$ClusterK4_kmeans,
                             "1" = "Mesenchymal",
                             "2" = "Proliferative",
                             "3" = "Immunoreactive",
                             "4" = "Differentiated")
everything$Subtype <- factor(everything$Subtype, levels = c("Mesenchymal", "Proliferative",
                                                            "Immunoreactive", "Differentiated"))

everything$Dataset <- factor(everything$Dataset, levels = c("TCGA", "microarray", "tothill"),
                             labels = c("TCGA RNA-seq", "TCGA Microarray", "Tothill"))

everything$Immune <- everything$`T cells` + everything$Macrophages + everything$Monocytes + everything$`Plasma cells` +
  everything$DC + everything$`NK cells` + everything$pDC + everything$`B cells` + everything$ILC + everything$`Mast cells`

# Make barplots
pA <- ggplot(everything) + 
    geom_boxplot(mapping = aes(x = Subtype, y = Fibroblasts, fill = Dataset)) +
    theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1)) +
    scale_fill_manual(values = colors_bulktypes) +
    labs(x = "Subtype (k = 4)", tag = "A")
pB <- ggplot(everything) + geom_boxplot(mapping = aes(x = Subtype, y = `Epithelial cells`, fill = Dataset)) +
    theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1)) +
    scale_fill_manual(values = colors_bulktypes) +
    labs(x = "Subtype (k = 4)", tag = "B")
pC <- ggplot(everything) + geom_boxplot(mapping = aes(x = Subtype, y = `Endothelial cells`, fill = Dataset)) +
    theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1)) +
    scale_fill_manual(values = colors_bulktypes) +
    labs(x = "Subtype (k = 4)", tag = "C")
pD <- ggplot(everything) + geom_boxplot(mapping = aes(x = Subtype, y = Immune, fill = Dataset)) +
    theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1)) +
    scale_fill_manual(values = colors_bulktypes) +
    labs(x = "Subtype (k = 4)", tag = "D")

top <- pA + pB + pC + pD + plot_layout(ncol = 4, guides = "collect")

pE <- ggplot(everything) + 
    geom_boxplot(mapping = aes(x = as.factor(ClusterK3_kmeans), y = Fibroblasts, fill = Dataset))  +
    scale_fill_manual(values = colors_bulktypes) +
    labs(x = "Cluster assignment (k = 3)", tag = "E")
pF <- ggplot(everything) + 
    geom_boxplot(mapping = aes(x = as.factor(ClusterK3_kmeans), y = `Epithelial cells`, fill = Dataset)) +
    scale_fill_manual(values = colors_bulktypes) +
    labs(x = "Cluster assignment (k = 3)", tag = "F")
pG <- ggplot(everything) +
    geom_boxplot(mapping = aes(x = as.factor(ClusterK3_kmeans), y = `Endothelial cells`, fill = Dataset)) +
    scale_fill_manual(values = colors_bulktypes) +
    labs(x = "Cluster assignment (k = 3)", tag = "G")
pH <- ggplot(everything) +
    geom_boxplot(mapping = aes(x = as.factor(ClusterK3_kmeans), y = Immune, fill = Dataset)) +
    scale_fill_manual(values = colors_bulktypes) +
    labs(x = "Cluster assignment (k = 3)", tag = "H")

middle <- pE + pF + pG + pH + plot_layout(ncol = 4, guides = "collect")


# Make Sankey plots 
sankey_plot <- function(everything, bulk_set) {
  dataset <- subset(everything, everything$Dataset == bulk_set)
  df <- dataset %>% make_long(ClusterK2_kmeans, ClusterK3_kmeans, Subtype) %>%
      dplyr::mutate(
          node = factor(node, levels = c(1, 2, 3, "Mesenchymal","Proliferative",
                                         "Immunoreactive", "Differentiated")),
          next_node = factor(next_node, levels = c(1, 2, 3, "Mesenchymal","Proliferative",
                                                   "Immunoreactive", "Differentiated"))
      )
  df$x <- recode(df$x,
                 "ClusterK2_kmeans" = "k = 2 Clusters",
                 "ClusterK3_kmeans" = "k = 3 Clusters",
                 "Subtype" = "k = 4 Clusters")
  
  ggplot(df, aes(x = x, 
                 next_x = next_x, 
                 node = node, 
                 next_node = next_node,
                 fill = factor(node),
                 label = node)) +
      geom_sankey(flow.alpha = .6,
                  node.color = "gray30") +
      scale_fill_manual(values = colors_subtypes) + 
      geom_sankey_label(size = 3, color = "white", fill = "gray40") +
      theme(axis.text.y=element_blank(),
            axis.ticks=element_blank(),
            axis.title.x = element_blank(),
            panel.grid = element_blank(),
            legend.position = "none") +
      ggtitle(bulk_set)
}

pI <- sankey_plot(everything, "TCGA RNA-seq") + labs(tag="I")
pJ <- sankey_plot(everything, "TCGA Microarray") + labs(tag="J")
pK <- sankey_plot(everything, "Tothill") + labs(tag="K")

bottom <- pI + pJ + pK + plot_layout(guides = "collect")

pdf(paste(figure_path, "figure3.pdf", sep = "/"), width = 16, height = 12, family = "sans")
top / middle / bottom + plot_layout(nrow=3, heights = c(4.5, 4.5, 3))
dev.off()
