# This is the first-look analysis to see if there's a difference in estimated
# cell type proportions for TCGA RNA-seq and TCGA microarray data. The answer:
# mostly, but also the Spearman correlation is high, so it may just be a
# platform difference.

suppressPackageStartupMessages({
  library(data.table)
  library(SingleCellExperiment)
  library(dplyr)
  library(curatedOvarianData)
  library(yaml)
  library(stringr)
  library(ggplot2)
  library(ggrepel)
})

params <- read_yaml("../../config.yml")
data_path <- params$data_path
plot_path <- params$plot_path
local_data_path <- params$local_data_path

microarray <- fread(paste(local_data_path, "deconvolution_output",
                    "microarray_default_bayesprism_results.tsv", sep = "/"))
rnaseq <- fread(paste(local_data_path, "deconvolution_output",
                      "TCGA_default_bayesprism_results.tsv", sep = "/"))

colnames(rnaseq) <- gsub("-", "\\.", colnames(rnaseq))
colnames(rnaseq) <- str_extract(colnames(rnaseq), "TCGA.\\w\\w.\\w\\w\\w\\w")

# Filter datasets down to only shared IDs and order them the same
common_samples <- intersect(colnames(rnaseq), colnames(microarray))
rnaseq <- subset(rnaseq, select = c(common_samples))
microarray <- subset(microarray, select = c(common_samples))

# Calculate the pearson correlation between cell type estimates
correlations <- numeric()
for (i in 2:length(common_samples)){
  correlations[i - 1] <- cor(rnaseq[, ..i], microarray[, ..i])
}
hist(correlations, main = "Correlation between TCGA RNAseq proportions and microarray proportions")

# Plot the low correlation values
low_corr <- which(correlations < 0.99)
low_corr <- low_corr + 1
low_corr <- c(1, low_corr)

low_rna <- rnaseq[, ..low_corr]
low_micro <- microarray[, ..low_corr]

melted_rna <- melt(low_rna, id.vars = "cell_type"); setnames(melted_rna, "value", "rnaseq")
melted_micro <- melt(low_micro, id.vars = "cell_type"); setnames(melted_micro, "value", "microarray")
melted_both <- left_join(melted_rna, melted_micro)

ggplot(melted_both, mapping = aes(x = microarray, y = rnaseq, color = cell_type)) +
  geom_point() + geom_label_repel(aes(label = cell_type, size = NULL)) + geom_abline()

# Plot all correlation values
melted_rna <- melt(rnaseq, id.vars = "cell_type"); setnames(melted_rna, "value", "rnaseq")
melted_micro <- melt(microarray, id.vars = "cell_type"); setnames(melted_micro, "value", "microarray")
melted_both <- left_join(melted_rna, melted_micro)

ggplot(melted_both, mapping = aes(x = microarray, y = rnaseq, color = cell_type)) +
  geom_point() + geom_abline()

# Calculate spearman correlation of individual cell types
rownames(rnaseq) <- rnaseq$cell_type; rnaseq$cell_type <- NULL
rownames(microarray) <- microarray$cell_type; microarray$cell_type <- NULL

tmp3 <- transpose(rnaseq)
rownames(tmp3) <- colnames(rnaseq)
colnames(tmp3) <- rownames(rnaseq)

tmp4 <- transpose(microarray)
rownames(tmp4) <- colnames(microarray)
colnames(tmp4) <- rownames(microarray)

cor(tmp3$`Epithelial cells`, tmp4$`Epithelial cells`, method = "spearman")
cor(tmp3$Fibroblasts, tmp4$Fibroblasts, method = "spearman")
cor(tmp3$Macrophages, tmp4$Macrophages, method = "spearman")
cor(tmp3$`Endothelial cells`, tmp4$`Endothelial cells`, method = "spearman")
