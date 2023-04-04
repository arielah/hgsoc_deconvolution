# Surgical debulking is one of the main treatments for ovarian cancer, and the
# outcome of surgery is a big determinant of survial. Surgeries are classified
# as optimal or suboptimal based on how much residal tumor is left afterward.
# Here, we're going to compare the percentage of fibroblasts in tumors with each
# outcome and see if there's a significant difference.

suppressPackageStartupMessages({
  library(data.table)
  library(SingleCellExperiment)
  library(dplyr)
  library(yaml)
  library(stringr)
  library(ggplot2)
  library(survival)
  library(ggfortify)
})

params <- read_yaml("../../config.yml")
data_path <- params$data_path
local_data_path <- params$local_data_path
plot_path <- params$plot_path

debulk_set <- fread(paste(local_data_path, "cluster_assignments", "AnalSet.csv", sep = "/"))

tcga <- fread(paste(local_data_path, "deconvolution_output",
                    "TCGA_default_bayesprism_results.tsv", sep = "/"))

cell_types <- tcga$cell_type

tcga$cell_type <- NULL
tcga_t <- t(as.matrix(tcga))
colnames(tcga_t) <- cell_types
tcga_t <- as.data.frame(tcga_t)

tcga_t$sampleid <- str_extract(rownames(tcga_t), "TCGA-\\w\\w-\\w\\w\\w\\w")
tcga_t$sampleid <- gsub("-", "\\.", tcga_t$sampleid)

tcga_t <- inner_join(tcga_t, debulk_set)

table(tcga_t$debulking)

ggplot(tcga_t, mapping = aes(y = Fibroblasts, x = debulking, group = debulking)) + geom_boxplot(notch = TRUE)

# Tothill
tothill <- fread(paste(local_data_path, "deconvolution_output",
                       "tothill_default_bayesprism_results.tsv", sep = "/"))
tothill$cell_type <- NULL
tothill_t <- t(as.matrix(tothill))
colnames(tothill_t) <- cell_types
tothill_t <- as.data.frame(tothill_t)

tothill_t$sampleid <- rownames(tothill_t)
tothill_t <- inner_join(tothill_t, debulk_set)

ggplot(tothill_t, mapping = aes(y = Fibroblasts, x = debulking, group = debulking)) + geom_boxplot(notch = TRUE)

# TCGA Microarray
microarray <- fread(paste(local_data_path, "deconvolution_output",
                       "microarray_default_bayesprism_results.tsv", sep = "/"))
microarray$cell_type <- NULL
microarray_t <- t(as.matrix(microarray))
colnames(microarray_t) <- cell_types
microarray_t <- as.data.frame(microarray_t)

microarray_t$sampleid <- rownames(microarray_t)
microarray_t <- inner_join(microarray_t, debulk_set)

ggplot(microarray_t, mapping = aes(y = Fibroblasts, x = debulking, group = debulking)) + geom_boxplot(notch = TRUE)

# All together
total <- rbind(tcga_t, tothill_t, microarray_t)
total$debulking <- recode(total$debulking,
       "suboptimal" = "Suboptimal",
       "optimal" = "Optimal")
total[is.na(total$debulking), ]$debulking <- "Missing"
ggplot(total, mapping = aes(y = Fibroblasts, x = debulking, group = debulking)) + geom_boxplot(notch = TRUE)
