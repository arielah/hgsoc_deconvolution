suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(stringr)
  library(dplyr)
  library(yaml)
})

source("figure_utils.R")

params <- read_yaml("../../config.yml")
data_path <- params$data_path
local_data_path <- params$local_data_path
plot_path <- params$plot_path
figure_path <- params$figure_path

load_dataset <- function(bulk_set, sc_set) {
  results <- fread(paste(local_data_path, "/deconvolution_output/",
                         bulk_set, "_", sc_set, "_bayesprism_results.tsv",
                         sep = ""))
  
  results_melt <- melt(results)
  setnames(results_melt, "value", sc_set)
  
  if (sc_set == "vazquez") {
      results_melt$cell_type <- recode(results_melt$cell_type,
                                       "Myeloid cells" = "Macrophages")
  }
  
  results_melt
}

tcga_ours <- load_dataset("TCGA", "fibro")
microarray_ours <- load_dataset("microarray", "fibro")
tothill_ours <- load_dataset("tothill", "fibro")
tcga_vazquez <- load_dataset("TCGA", "vazquez")
microarray_vazquez <- load_dataset("microarray", "vazquez")
tothill_vazquez <- load_dataset("tothill", "vazquez")

tcga <- full_join(tcga_ours, tcga_vazquez)
pA <- ggplot(tcga, mapping = aes(x = fibro, y = vazquez, color = cell_type)) + 
  geom_point(size = 0.8) + geom_abline() +
  scale_color_manual(values = colors_celltypes) +
  labs(x = "Proportions with HGSOC Penn/Utah reference", 
       y = "Proportions with Vázquez-García reference",
       title = "TCGA RNA-seq",
       color = "Cell type", tag = "A")


microarray <- full_join(microarray_ours, microarray_vazquez)
pB <- ggplot(microarray, mapping = aes(x = fibro, y = vazquez, color = cell_type)) +
  geom_point(size = 0.8) + geom_abline() +
  scale_color_manual(values = colors_celltypes) +
  labs(x = "Proportions with HGSOC Penn/Utah reference", 
       y = "Proportions with Vázquez-García reference",
       title = "TCGA Microarray",
       color = "Cell type", tag = "B")

tothill <- full_join(tothill_ours, tothill_vazquez)
pC <- ggplot(tothill, mapping = aes(x = fibro, y = vazquez, color = cell_type)) +
  geom_point(size = 0.8) + geom_abline() +
  scale_color_manual(values = colors_celltypes) +
  labs(x = "Proportions with HGSOC Penn/Utah reference", 
       y = "Proportions with Vázquez-García reference",
       title = "Tothill",
       color = "Cell type", tag = "C")

# Compare TCGA and microarray

tcga_ours$variable <- gsub("-", "\\.", tcga_ours$variable)
tcga_ours$variable <- str_extract(tcga_ours$variable, "TCGA.\\w\\w.\\w\\w\\w\\w")

setnames(tcga_ours, "fibro", "rnaseq")
setnames(microarray_ours, "fibro", "microarray")

tcga_both <- inner_join(tcga_ours, microarray_ours)
pD <- ggplot(tcga_both, mapping = aes(x = microarray, y = rnaseq, color = cell_type)) +
  geom_point(size = 0.8) + geom_abline() +
  scale_color_manual(values = colors_celltypes) +
  labs(x = "TCGA Microarray Proportions", y = "TCGA RNA-seq Proportions",
       color = "Cell type", tag = "D",
       title = "Cross-plaform comparison")

pdf(paste(figure_path, "figure2.pdf", sep = "/"), width = 16, height = 12, family = "sans")
pA + pB + pC + pD
dev.off()
