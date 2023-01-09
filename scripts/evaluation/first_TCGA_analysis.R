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

# Plot early TCGA results
tcga <- fread(paste(local_data_path, "deconvolution_output", 
                    "TCGA_default_bayesprism_results.tsv", sep = "/"))
tcga_melt <- melt(tcga)

plotfile <- paste(plot_path, "deconvolution_plots",
                  "TCGA_proportion_barchart.png", sep = "/")
png(filename = plotfile, width = 1200)
ggplot(tcga_melt, mapping = aes(x=variable, y=value, fill=cell_type)) +
    geom_bar(stat = "identity") + 
    theme(axis.text.x=element_blank(), #remove x axis labels
          axis.ticks.x=element_blank(), #remove x axis ticks
          axis.text.y=element_blank(),  #remove y axis labels
          axis.ticks.y=element_blank()) #remove y axis ticks
dev.off()
    
# Switch so cell types are columns and samples are rows for easier analysis
cell_types <- tcga$cell_type

tcga$cell_type <- NULL
tcga_t <- t(as.matrix(tcga))
colnames(tcga_t) <- cell_types
tcga_t <- as.data.frame(tcga_t)

# Load survival data
tcga_survival <- fread(paste(local_data_path, "TCGA", "TCGA_OV_survival.tsv", 
                             sep = "/"))
tcga_patients <- str_extract(colnames(tcga), "TCGA-\\w\\w-\\w\\w\\w\\w")
tcga_survival <- subset(tcga_survival, tcga_survival$bcr_patient_barcode %in%
                            tcga_patients)

# Combine survival data with %
tcga_t$bcr_patient_barcode <- tcga_patients
tcga_master <- full_join(tcga_survival, tcga_t)

plot(tcga_master$OS.time, tcga_master$Fibroblasts, xlab="OS (days)", ylab="% Fibroblasts")
plot(tcga_master$OS.time, tcga_master$Macrophages, xlab="OS (days)", ylab="% Macrophages")

# Put the samples into quartiles based on fibroblast content
quantiles <- quantile(tcga_master$Fibroblasts)
q1 <- quantiles[2]
q3 <- quantiles[4]
tcga_master$high_fibro <- ifelse(tcga_master$Fibroblasts > q3, 1, 0)

# Get Kaplan-Meier curves
km <- Surv(time = tcga_master$PFI.time,  event = tcga_master$PFI)
km_treatment<-survfit(km~high_fibro,data=tcga_master,type='kaplan-meier',conf.type='log')

autoplot(km_treatment)

# Get subtype annotations
cluster_file <- paste(local_data_path, "cluster_assignments", "FullClusterMembership.csv", sep = "/")
cluster_list <- fread(cluster_file)

cluster_list$V1 <- gsub("\\.", "-", cluster_list$V1)
setnames(cluster_list, "V1", "ID")

tcga_t$ID <- str_extract(rownames(tcga_t), "TCGA-\\w\\w-\\w\\w\\w\\w")

tcga_t <- left_join(tcga_t, cluster_list)

tcga_t$Subtype <- recode(tcga_t$ClusterK4_kmeans,
                           "1" = "Mesenchymal",
                           "2" = "Proliferative",
                           "3" = "Immunoreactive",
                           "4" = "Differentiated")

# Compare cell type proportions of subtypes
ggplot(tcga_t, mapping = aes(x=Subtype, y=Fibroblasts)) + geom_boxplot()
ggplot(tcga_t, mapping = aes(x=Subtype, y=`Epithelial cells`)) + geom_boxplot()
ggplot(tcga_t, mapping = aes(x=Subtype, y=Macrophages)) + geom_boxplot()
ggplot(tcga_t, mapping = aes(x=Subtype, y=`Endothelial cells`)) + geom_boxplot()

# Check for pan-immune fraction
tcga_t$Immune <- tcga_t$`T cells` + tcga_t$Macrophages + tcga_t$Monocytes + tcga_t$`Plasma cells` +
  tcga_t$DC + tcga_t$`NK cells` + tcga_t$pDC + tcga_t$`B cells` + tcga_t$ILC + tcga_t$`Mast cells`

ggplot(tcga_t, mapping = aes(x=Subtype, y=Immune)) + geom_boxplot()