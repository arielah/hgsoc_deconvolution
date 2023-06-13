suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(stringr)
  library(pheatmap)
  library(dplyr)
    library(patchwork)
  library(ComplexHeatmap)
  library(yaml)
    library(ggrepel)
    library(gridExtra)
    library(ggplotify)
    library(GSEABase)
    library(GSVA)
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
cluster_file <- paste(local_data_path, "cluster_assignments", "FullClusterMembership.csv", sep = "/")
cluster_list <- fread(cluster_file)
cluster_list <- cluster_list[order(cluster_list$ClusterK4_kmeans), ]

cluster_list$ClusterK4_kmeans <- recode(cluster_list$ClusterK4_kmeans,
                                        "1" = "Mesenchymal",
                                        "2" = "Proliferative",
                                        "3" = "Immunoreactive",
                                        "4" = "Differentiated")

plot_weights <- function(ebd.res, name, tag) {
  weights <- ebd.res$omega
  
  if(name == "TCGA RNA-seq") {
    rownames(weights) <- str_extract(rownames(weights), "TCGA-\\w\\w-\\w\\w\\w\\w")
    rownames(weights) <- gsub("-", "\\.", rownames(weights))
  }
  
  clusters <- subset(cluster_list, cluster_list$V1 %in% rownames(weights))
  weights <- weights[rownames(weights) %in% cluster_list$V1, ]
  weights <- weights[match(clusters$V1, rownames(weights)), ]
  
  row_ha <- rowAnnotation(
    k4_kmeans = clusters$ClusterK4_kmeans,
    k3_kmeans = factor(clusters$ClusterK3_kmeans),
    col = list(k3_kmeans = colors_subtypes[5:7],
               k4_kmeans = colors_subtypes[1:4]),
    annotation_label = list(k3_kmeans = "Cluster (k = 3)",
                            k4_kmeans = "Subtype (k = 4)")
  )
  
  as.ggplot(Heatmap(weights, left_annotation = row_ha, show_row_names = F,
          cluster_rows = F, cluster_columns = F, col = heatmap_scale_2d,
<<<<<<< HEAD
          heatmap_legend_param = list(title = " ", at = c(0, 0.2, 0.4, 0.6, 0.8, 1)),
          name = name,
          column_names_rot = 45
          )) +
=======
          heatmap_legend_param = list(title = " ", at = c(0,0.2,0.4,0.6,0.8,1)),
          name = name)) +
>>>>>>> 91ac905c3e97e2a63a936b95a0d22b9eba81ba61
      labs(tag = tag) + ggtitle(name) + 
      theme(plot.title = element_text(hjust = 0.43))
}


hallmark_sets <- getGmt("~/Downloads/h.all.v2023.1.Hs.symbols.gmt")

plot_gsva <- function(ebd.res, name, tag) {
    gene_expression <- t(ebd.res$eta_post)
    gsva.es <- gsva(gene_expression, hallmark_sets)
    rownames(gsva.es) <- gsub("HALLMARK_", "", rownames(gsva.es))
    rownames(gsva.es)[rownames(gsva.es) == 
                          "REACTIVE_OXYGEN_SPECIES_PATHWAY"] <- "REACTIVE_OXYGEN_SPECIES"
    
    as.ggplot(Heatmap(gsva.es, col = heatmap_scale_2d, cluster_columns = F,
                      heatmap_legend_param = list(title = " ",
<<<<<<< HEAD
                                                  at = c(-1, -0.5, 0, 0.5, 1)), 
                      name = name, column_names_rot = 45,
=======
                                                  at = c(-1,-0.5,0,0.5,1)), 
                      name = name,
>>>>>>> 91ac905c3e97e2a63a936b95a0d22b9eba81ba61
                      row_names_gp = grid::gpar(fontsize = 10))) +
        labs(tag = tag) + ggtitle(name) + 
        theme(plot.title = element_text(hjust = 0.3))
}

pA <- plot_weights(tcga, "TCGA RNA-seq", "A") 
pB <- plot_weights(microarray, "TCGA Microarray", "B")
pC <- plot_weights(tothill, "Tothill", "C")

pD <- plot_gsva(tcga, "TCGA RNA-seq", "D")
pE <- plot_gsva(microarray, "TCGA Microarray", "E")
pF <- plot_gsva(tothill, "Tothill", "F")

pdf(paste(figure_path, "figure5.pdf", sep = "/"), width = 16, height = 24, family = "sans")
pA + pD + pB + pE + pC + pF + plot_layout(ncol = 2)
dev.off()


tcga_4 <- load_data("TCGA", 4)
microarray_4 <- load_data("microarray", 4)
tothill_4 <- load_data("tothill", 4)

qA <- plot_weights(tcga_4, "TCGA RNA-seq", "A") 
qB <- plot_weights(microarray_4, "TCGA Microarray", "B")
qC <- plot_weights(tothill_4, "Tothill", "C")

qD <- plot_gsva(tcga_4, "TCGA RNA-seq", "D")
qE <- plot_gsva(microarray_4, "TCGA Microarray", "E")
qF <- plot_gsva(tothill_4, "Tothill", "F")

pdf(paste(figure_path, "suppfig3.pdf", sep = "/"), width = 16, height = 24, family = "sans")
<<<<<<< HEAD
qA + qD + qB + qE + qC + qF + plot_layout(ncol = 2)
=======
qA + qD + qB + qE + qC + qF + plot_layout(ncol=2)
>>>>>>> 91ac905c3e97e2a63a936b95a0d22b9eba81ba61
dev.off()
