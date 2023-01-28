suppressPackageStartupMessages({
  library(data.table)
  library(SingleCellExperiment)
  library(dplyr)
  library(curatedOvarianData)
  library(yaml)
})

params <- read_yaml("../../config.yml")
data_path <- params$data_path
local_data_path <- params$local_data_path

# Load in microarray data
data("TCGA_eset")
tcga <- exprs(TCGA_eset)
pheno <- pData(TCGA_eset)

# Filter out non-HGSOC samples
hgsoc_samples <- pheno$histological_type=="ser" & pheno$sample_type=="tumor" & pheno$summarygrade=="high"
hgsoc_samples[is.na(hgsoc_samples)] <- FALSE
tcga <- tcga[, hgsoc_samples]

# Load a small SCE object to merge gene names
sce <- readRDS(paste(local_data_path, "sce_objects/19595X1_labeled.rds", sep = "/"))
common_genes <- intersect(rownames(tcga), rowData(sce)$Symbol)

tcga <- tcga[rownames(tcga) %in% common_genes,]

# Convert from log intensity
tcga <- 2^tcga

tcga <- cbind(rownames(tcga), tcga)
colnames(tcga)[1] <- "Gene"

# Save expression file
outfile <- paste(local_data_path, "deconvolution_input",
                 "bulk_data_microarray.tsv", sep= "/")
write.table(tcga, outfile, sep = "\t", row.names = F, quote = F)

# Make survival file
pheno <- pheno[hgsoc_samples, ]
pheno <- subset(pheno, select = c("unique_patient_ID", "grade", "age_at_initial_pathologic_diagnosis",
                                  "days_to_tumor_recurrence", "recurrence_status", "days_to_death",
                                  "vital_status"))

outfile <- paste(local_data_path, "TCGA",
                 "microarray_survival.tsv", sep = "/")
write.table(pheno, outfile, row.names = F, sep = "\t", quote = F)
