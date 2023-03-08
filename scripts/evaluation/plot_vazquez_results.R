
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

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 1) {
    print(args)
    dataset <- args[1]
}

data <- fread(paste(local_data_path, "/deconvolution_output/", dataset,
                     "_vazquez_bayesprism_results.tsv", sep = ""))
data_melt <- melt(data)

g <- ggplot(data_melt, mapping = aes(x=variable, y=value, fill=cell_type, color=cell_type)) +
    geom_bar(stat = "identity") + 
    theme(axis.text.x=element_blank(), #remove x axis labels
          axis.ticks.x=element_blank(), #remove x axis ticks
          axis.text.y=element_blank(),  #remove y axis labels
          axis.ticks.y=element_blank()) #remove y axis ticks

plotfile <- paste(plot_path, "/evaluation_plots/", dataset,
                  "_proportion_barchart_vazquez.png", sep = "")
png(filename = plotfile, width = 1200)
g
dev.off()



# Switch so cell types are columns and samples are rows for easier analysis
cell_types <- data$cell_type

data$cell_type <- NULL
data_t <- t(as.matrix(data))
colnames(data_t) <- cell_types
data_t <- as.data.frame(data_t)
data_t <- cbind(rownames(data_t), data_t)
setnames(data_t, "rownames(data_t)", "ID")

if(dataset == "TCGA") {
    data_t$ID <- str_extract(data_t$ID, "TCGA-\\w\\w-\\w\\w\\w\\w")
    data_t$ID <- gsub("-","\\.", data_t$ID)
}

survival <- fread(paste(local_data_path, "cluster_assignments",
                        "AnalSet.csv", sep = "/"))
setnames(survival, "sampleid", "ID")

data_t <- left_join(data_t, survival)

data_t$Immune <- data_t$`T cells` + data_t$`B cells` + data_t$`Myeloid cells`

# Put the samples into quartiles based on fibroblast content
hist(data_t$Fibroblasts)
quantiles <- quantile(data_t$Fibroblasts)
q1 <- quantiles[2]
q3 <- quantiles[4]
data_t$high_fibro <- ifelse(data_t$Fibroblasts > q3, 1, 0)


# Get Kaplan-Meier curves
km <- Surv(time = data_t$months,  event = data_t$vital)
km_treatment<-survfit(km~high_fibro,data=data_t,type='kaplan-meier',conf.type='log')

plotfile <- paste(plot_path, "/evaluation_plots/", dataset, "_KaplanMeier_fibroblasts_vazquez.png", sep = "")
png(filename = plotfile)
autoplot(km_treatment)
dev.off()


cluster_file <- paste(local_data_path, "cluster_assignments", "FullClusterMembership.csv", sep = "/")
cluster_list <- fread(cluster_file)

cluster_list$V1 <- gsub("\\.", "-", cluster_list$V1)
setnames(cluster_list, "V1", "ID")
cluster_list$Dataset <- NULL

data_t <- left_join(data_t, cluster_list)

data_t <- subset(data_t, !is.na(data_t$ClusterK4_kmeans))

data_t$ClusterK4_kmeans <- recode(data_t$ClusterK4_kmeans,
                                  "1" = "Mesenchymal",
                                  "2" = "Proliferative",
                                  "3" = "Immunoreactive",
                                  "4" = "Differentiated")

# Compare cell type proportions of subtypes
g <- ggplot(data_t, mapping = aes(x=ClusterK4_kmeans, y=Fibroblasts)) + geom_boxplot()
plotfile <- paste(plot_path, "/evaluation_plots/", dataset, "_fibroblasts_by_subtype_vazquez.png", sep = "")
png(filename = plotfile); g; dev.off()
g

g <- ggplot(data_t, mapping = aes(x=ClusterK4_kmeans, y=`Epithelial cells`)) + geom_boxplot()
plotfile <- paste(plot_path, "/evaluation_plots/", dataset, "_epithelial_by_subtype_vazquez.png", sep = "")
png(filename = plotfile); g; dev.off()
g

g <- ggplot(data_t, mapping = aes(x=ClusterK4_kmeans, y=Immune)) + geom_boxplot()
plotfile <- paste(plot_path, "/evaluation_plots/", dataset, "_immune_by_subtype_vazquez.png", sep = "")
png(filename = plotfile); g; dev.off()
g

g <- ggplot(data_t, mapping = aes(x=ClusterK4_kmeans, y=`Myeloid cells`)) + geom_boxplot()
plotfile <- paste(plot_path, "/evaluation_plots/", dataset, "_myeloid_by_subtype_vazquez.png", sep = "")
png(filename = plotfile); g; dev.off()
g

g <- ggplot(data_t, mapping = aes(x=ClusterK4_kmeans, y=`T cells`)) + geom_boxplot()
plotfile <- paste(plot_path, "/evaluation_plots/", dataset, "_tcells_by_subtype_vazquez.png", sep = "")
png(filename = plotfile); g; dev.off()
g

g <- ggplot(data_t, mapping = aes(x=ClusterK4_kmeans, y=`Endothelial cells`)) + geom_boxplot()
plotfile <- paste(plot_path, "/evaluation_plots/", dataset, "_endothelial_by_subtype_vazquez.png", sep = "")
png(filename = plotfile); g; dev.off()
g

g <- ggplot(data_t, mapping = aes(x=as.factor(ClusterK3_kmeans), y=Fibroblasts)) + geom_boxplot()
plotfile <- paste(plot_path, "/evaluation_plots/", dataset, "_fibroblasts_by_subtype_vazquez_k3.png", sep = "")
png(filename = plotfile); g; dev.off()
g

g <- ggplot(data_t, mapping = aes(x=factor(ClusterK3_kmeans), y=`Epithelial cells`)) + geom_boxplot()
plotfile <- paste(plot_path, "/evaluation_plots/", dataset, "_epithelial_by_subtype_vazquez_k3.png", sep = "")
png(filename = plotfile); g; dev.off()
g 

g <- ggplot(data_t, mapping = aes(x=factor(ClusterK3_kmeans), y=Immune)) + geom_boxplot()
plotfile <- paste(plot_path, "/evaluation_plots/", dataset, "_immune_by_subtype_vazquez_k3.png", sep = "")
png(filename = plotfile); g; dev.off()
g

g <- ggplot(data_t, mapping = aes(x=factor(ClusterK3_kmeans), y=`Myeloid cells`)) + geom_boxplot()
plotfile <- paste(plot_path, "/evaluation_plots/", dataset, "_myeloid_by_subtype_vazquez_k3.png", sep = "")
png(filename = plotfile); g; dev.off()
g

g <- ggplot(data_t, mapping = aes(x=factor(ClusterK3_kmeans), y=`T cells`)) + geom_boxplot()
plotfile <- paste(plot_path, "/evaluation_plots/", dataset, "_tcells_by_subtype_vazquez_k3.png", sep = "")
png(filename = plotfile); g; dev.off()
g

g <- ggplot(data_t, mapping = aes(x=factor(ClusterK3_kmeans), y=`Endothelial cells`)) + geom_boxplot()
plotfile <- paste(plot_path, "/evaluation_plots/", dataset, "_endothelial_by_subtype_vazquez_k3.png", sep = "")
png(filename = plotfile); g; dev.off()
g