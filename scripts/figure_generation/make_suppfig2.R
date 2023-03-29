suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(stringr)
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

tcga <- readRDS(paste(local_data_path, 
                      "/embeddings/TCGA_malignant_nmf_rank.rds", sep = ""))
microarray <- readRDS(paste(local_data_path,
                            "/embeddings/microarray_malignant_nmf_rank.rds", sep = ""))
tothill <- readRDS(paste(local_data_path,
                         "/embeddings/tothill_malignant_nmf_rank.rds", sep = ""))

make_heatmap <- function(dataset, k) {
  thing <- dataset$consensus[[k-1]]
  as.ggplot(Heatmap(thing, show_column_names = F, show_row_names = F,
                    heatmap_legend_param = list(title = " "),
                    show_column_dend = F, show_row_dend = F))
}

make_cophenetic <- function(dataset) {
  thing <- dataset$measures
  ggplot(thing, mapping = aes(x = rank, y = cophenetic)) + geom_line() +
    scale_x_continuous(breaks=2:12) + xlab("k") + ylab("Cophenetic correlation")
}

pA <- make_heatmap(tcga, 2) + labs(tag = "A")
pB <- make_heatmap(tcga, 3)
pC <- make_heatmap(tcga, 4)
pD <- make_heatmap(tcga, 5)
pE <- make_heatmap(tcga, 6)
pF <- make_cophenetic(tcga)

top <- pA + pB + pC + pD + pE + pF

pG <- make_heatmap(microarray, 2) + labs(tag = "B")
pH <- make_heatmap(microarray, 3)
pI <- make_heatmap(microarray, 4)
pJ <- make_heatmap(microarray, 5)
pK <- make_heatmap(microarray, 6)
pL <- make_cophenetic(microarray)

middle <- pG + pH + pI + pJ + pK + pL

pM <- make_heatmap(tothill, 2) + labs(tag = "C")
pN <- make_heatmap(tothill, 3)
pO<- make_heatmap(tothill, 4)
pP <- make_heatmap(tothill, 5)
pQ <- make_heatmap(tothill, 6)
pR <- make_cophenetic(tothill)

bottom <- pM + pN + pO + pP + pQ + pR

pdf(paste(figure_path, "suppfig2.pdf", sep = "/"), width = 7, height = 12, family = "sans")
top / middle / bottom
dev.off()
