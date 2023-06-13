suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(stringr)
    library(patchwork)
    library(dplyr)
    library(survival)
    library(ggfortify)
    library(yaml)
})

source("figure_utils.R")

params <- read_yaml("../../config.yml")
data_path <- params$data_path
local_data_path <- params$local_data_path
plot_path <- params$plot_path
figure_path <- params$figure_path

load_datasets <- function(bulk_set, sc_set) {
    results_t <- fread(paste(local_data_path, "/deconvolution_output/", bulk_set,
                           "_", sc_set, "_bayesprism_results.tsv", sep = ""))
    
    cell_types <- results_t$cell_type
    
    results_t$cell_type <- NULL
    results <- t(as.matrix(results_t))
    colnames(results) <- cell_types
    results <- as.data.frame(results)
    results <- cbind(rownames(results), results)
    setnames(results, "rownames(results)", "ID")
    results$dataset <- bulk_set
    
    if (bulk_set == "TCGA"){
        results$ID <- str_extract(results$ID, "TCGA-\\w\\w-\\w\\w\\w\\w")
        results$ID <- gsub("-", "\\.", results$ID)
    }
    
    results$Immune <- results$`T cells` + results$Macrophages + results$Monocytes +
        results$`Plasma cells` + results$`Mast cells` + results$`NK cells` +
        results$pDC + results$`B cells` + results$ILC + results$DC 
    
    results
}

tcga <- load_datasets("TCGA", "fibro")
microarray <- load_datasets("microarray", "fibro")
tothill <- load_datasets("tothill", "fibro")

# Load survival data
covariates <- fread(paste(local_data_path, "cluster_assignments",
                          "AnalSet.csv", sep = "/"))
covariates$V1 <- NULL
setnames(covariates, "sampleid", "ID")
covariates$debulking <- recode(covariates$debulking,
                               "optimal" = "Optimal",
                               "suboptimal" = "Suboptimal")

# Merge survival and composition data
tcga <- inner_join(tcga, covariates)
microarray <- inner_join(microarray, covariates)
tothill <- inner_join(tothill, covariates)

composition <- rbind(tcga, microarray, tothill)

kaplan_meier <- function(data, title, tag) {
    quantiles <- quantile(data$Fibroblasts)
    q1 <- quantiles[2]
    median <- quantiles[3]
    q3 <- quantiles[4]
    data$high_fibro <- ifelse(data$Fibroblasts > q3, "High fibroblast", "Other")
    
    km <- Surv(data$months, data$vital)
    km_treatment<-survfit(km~high_fibro,data=data,type='kaplan-meier',conf.type='log')
    
    autoplot(km_treatment, conf.int = F) +
        labs(x="Time since diagnosis (months)", y = "Survival",
             title = title, tag = tag) +
        theme(legend.title=element_blank()) +
        scale_color_manual(values=c("#F8766D","#B79F00"))
    }

pA <- kaplan_meier(tcga, "TCGA RNA-seq", "A")
pB <- kaplan_meier(microarray, "TCGA Microarray", "B")
pC <- kaplan_meier(tothill, "Tothill", "C")

pD <- kaplan_meier(composition, "All Datasets", "D")

top <- pA + pB + pC + pD + plot_layout(nrow = 2, guides = "collect")

surgery <- subset(composition, !is.na(composition$debulking))

surgery$dataset <- recode(surgery$dataset,
                          "microarray" = "TCGA Microarray",
                          "TCGA" = "TCGA RNA-seq",
                          "tothill" = "Tothill")

pE <- ggplot(surgery, mapping = aes(y=Fibroblasts, x=dataset, fill=debulking)) +
    geom_boxplot(notch=TRUE) +
    scale_fill_manual(values = c("#39B600", "#00B0F6")) +
    labs(x="Dataset", fill="Debulking status", tag="E")

pF <- ggplot(surgery, mapping = aes(y=Immune, x=dataset, fill=debulking)) +
    geom_boxplot(notch=TRUE) +
    scale_fill_manual(values = c("#39B600", "#00B0F6")) +
    labs(x="Dataset", fill="Debulking status", tag="F")

bottom <- pE + pF + plot_layout(guides="collect")

pdf(paste(figure_path, "figure4.pdf", sep = "/"), width = 12, height = 14, family = "sans")
top / bottom + plot_layout(heights = c(2, 1))
dev.off()


plot_histogram <- function(data, bulk_set, cell_type, lab) {
    if (cell_type == "Fibroblasts") {
        data$celltype <- data$Fibroblasts
    } else if (cell_type == "Immune") {
        data$celltype <- data$Immune
    }
    quantiles <- quantile(data$celltype)
    q3 <- quantiles[4]
    print(q3)
    
    ggplot(data, mapping = aes(x=celltype))  + geom_histogram() +
        geom_vline(xintercept = q3, linetype = "dashed", color = "red") +
        labs(title = bulk_set, tag = lab, x = cell_type, y = "Count")
}

qA <- plot_histogram(tcga, "TCGA RNA-seq", "Fibroblasts", "A")
qB <- plot_histogram(microarray, "TCGA Microarray", "Fibroblasts", "B")
qC <- plot_histogram(tothill, "Tothill", "Fibroblasts", "C")

qD <- plot_histogram(tcga, "TCGA RNA-seq", "Immune", "D")
qE <- plot_histogram(microarray, "TCGA Microarray", "Immune", "E")
qF <- plot_histogram(tothill, "Tothill", "Immune", "F")

pdf(paste(figure_path, "suppfig1.pdf", sep = "/"), width = 18, height = 12, family = "sans")
qA + qB + qC + qD + qE + qF + plot_layout(nrow = 2)
dev.off()


kaplan_meier_immune <- function(data, title, tag) {
    quantiles <- quantile(data$Immune)
    q1 <- quantiles[2]
    median <- quantiles[3]
    q3 <- quantiles[4]
    data$high_immune <- ifelse(data$Immune > q3, "High immune", "Other")

    km <- Surv(data$months, data$vital)
    km_treatment<-survfit(km~high_immune,data=data,type='kaplan-meier',conf.type='log')

    autoplot(km_treatment, conf.int = F) +
        labs(x="Time since diagnosis (months)", y = "Survival",
             title = title, tag = tag) +
        theme(legend.title=element_blank()) +
        scale_color_manual(values=c("#F8766D","#B79F00"))
}

rA <- kaplan_meier_immune(tcga, "TCGA RNA-seq", "A")
rB <- kaplan_meier_immune(microarray, "TCGA Microarray", "B")
rC <- kaplan_meier_immune(tothill, "Tothill", "C")
rD <- kaplan_meier_immune(composition, "All Datasets", "D")

pdf(paste(figure_path, "suppfig2.pdf", sep = "/"), width = 12, height = 9.3, family = "sans")
rA + rB + rC + rD + plot_layout(nrow = 2, guides = "collect")
dev.off()
