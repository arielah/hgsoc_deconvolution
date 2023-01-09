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


# Plot early AACES results
aaces <- fread(paste(local_data_path, "deconvolution_output",
                     "AACES_default_bayesprism_results.tsv", sep = "/"))
aaces_melt <- melt(aaces)

plotfile <- paste(plot_path, "deconvolution_plots",
                  "AACES_proportion_barchart.png", sep = "/")
png(filename = plotfile, width = 1200)
ggplot(aaces_melt, mapping = aes(x=variable, y=value, fill=cell_type)) +
    geom_bar(stat = "identity") + 
    theme(axis.text.x=element_blank(), #remove x axis labels
          axis.ticks.x=element_blank(), #remove x axis ticks
          axis.text.y=element_blank(),  #remove y axis labels
          axis.ticks.y=element_blank()) #remove y axis ticks
dev.off()

# Switch so cell types are columns and samples are rows for easier analysis
cell_types <- aaces$cell_type

aaces$cell_type <- NULL
aaces_t <- t(as.matrix(aaces))
colnames(aaces_t) <- cell_types
aaces_t <- as.data.frame(aaces_t)
aaces_t <- cbind(rownames(aaces_t), aaces_t)
setnames(aaces_t, "rownames(aaces_t)", "ID")

# Load survival data
aaces_survival <- fread(paste(local_data_path, "AACES", "AACES_survival.tsv",
                              sep = "/"))
aaces_t <- left_join(aaces_t, aaces_survival)

plot(aaces_t$survival_days, aaces_t$Fibroblasts, xlab="OS (days)", ylab="% Fibroblasts")
plot(aaces_t$survival_days, aaces_t$`Endothelial cells`, xlab="OS (days)", ylab="% Endothelial")
plot(aaces_t$survival_days, aaces_t$`T cells`, xlab="OS (days)", ylab="% T cells")
plot(aaces_t$survival_days, aaces_t$`B cells`, xlab="OS (days)", ylab="% B cells")

# Put the samples into quartiles based on fibroblast content
quantiles <- quantile(aaces_t$Fibroblasts)
q1 <- quantiles[2]
q3 <- quantiles[4]
aaces_t$high_fibro <- ifelse(aaces_t$Fibroblasts > q3, 1, 0)

# Get Kaplan-Meier curves
aaces_t$vitalstatus <- ifelse(aaces_t$vitalstatus=="Alive", 0, 1)
km <- Surv(time = aaces_t$survival_days,  event = aaces_t$vitalstatus)
km_treatment<-survfit(km~high_fibro,data=aaces_t,type='kaplan-meier',conf.type='log')

autoplot(km_treatment)

# Compare cell type proportions of subtypes
ggplot(aaces_t, mapping = aes(x=ClusterK4_kmeans_TCGA_names, y=Fibroblasts)) + geom_boxplot()
ggplot(aaces_t, mapping = aes(x=ClusterK4_kmeans_TCGA_names, y=`Epithelial cells`)) + geom_boxplot()
ggplot(aaces_t, mapping = aes(x=ClusterK4_kmeans_TCGA_names, y=`T cells`)) + geom_boxplot()
ggplot(aaces_t, mapping = aes(x=ClusterK4_kmeans_TCGA_names, y=`Endothelial cells`)) + geom_boxplot()

# Check for pan-immune fraction
aaces_t$Immune <- aaces_t$`T cells` + aaces_t$Macrophages + aaces_t$Monocytes + aaces_t$`Plasma cells` +
  aaces_t$DC + aaces_t$`NK cells` + aaces_t$pDC + aaces_t$`B cells` + aaces_t$ILC + aaces_t$`Mast cells`

ggplot(aaces_t, mapping = aes(x=ClusterK4_kmeans_TCGA_names, y=Immune)) + geom_boxplot()
