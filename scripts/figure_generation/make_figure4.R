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

params <- read_yaml("../../config.yml")
data_path <- params$data_path
local_data_path <- params$local_data_path
plot_path <- params$plot_path
figure_path <- params$figure_path

theme_set(theme_bw() +
              theme(text = element_text(size = 14),
                    strip.background = element_rect(colour = NA,
                                                    fill = "white"),
                    plot.title = element_text(hjust = 0.5)
              )
)

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
    quantiles <- quantile(data$Immune)
    q1 <- quantiles[2]
    median <- quantiles[3]
    q3 <- quantiles[4]
    data$high_fibro <- ifelse(data$Immune > q3, "High immune", "Other")
    
    km <- Surv(data$months, data$vital)
    km_treatment<-survfit(km~high_fibro,data=data,type='kaplan-meier',conf.type='log')
    
    autoplot(km_treatment, conf.int = F) +
        labs(x="Time since diagnosis (months)", y = "Survival",
             title = title, tag = tag) +
        theme(legend.title=element_blank())
    }

pA <- kaplan_meier(tcga, "TCGA RNA-seq", "A")
pB <- kaplan_meier(microarray, "TCGA Microarray", "B")
pC <- kaplan_meier(tothill, "Tothill", "C")

pD <- kaplan_meier(composition, "All Datasets", "D")

top <- pA + pB + pC + pD + plot_layout(nrow = 2, guides = "collect")

surgery <- subset(composition, !is.na(composition$debulking))

pE <- ggplot(surgery, mapping = aes(y=Fibroblasts, x=dataset, fill=debulking)) +
    geom_boxplot(notch=TRUE) +
    labs(x="Dataset", fill="Debulking status", tag="E")

pF <- ggplot(surgery, mapping = aes(y=Immune, x=dataset, fill=debulking)) +
    geom_boxplot(notch=TRUE) +
    labs(x="Dataset", fill="Debulking status", tag="F")

bottom <- pE + pF + plot_layout(guides="collect")

pdf(paste(figure_path, "figure4.pdf", sep = "/"), width = 12, height = 14, family = "sans")
top / bottom + plot_layout(heights = c(2, 1))
dev.off()

# Delete this later, checking the notch
#optimal <- subset(surgery, surgery$debulking=="optimal")
#suboptimal <- subset(surgery, surgery$debulking=="suboptimal")

#median(optimal$Fibroblasts) + 1.5 * IQR(optimal$Fibroblasts) / sqrt(nrow(optimal))
#median(suboptimal$Fibroblasts) - 1.5 * IQR(suboptimal$Fibroblasts) / sqrt(nrow(suboptimal))
