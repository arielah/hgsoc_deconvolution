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
data("GSE9891_eset")
tothill <- exprs(GSE9891_eset)
pheno <- pData(GSE9891_eset)

# Filter out non-HGSOC samples
hgsoc_samples <- pheno$histological_type=="ser" & pheno$sample_type=="tumor" & pheno$summarygrade=="high"
hgsoc_samples[is.na(hgsoc_samples)] <- FALSE
tothill <- tothill[, hgsoc_samples]

# Load a small SCE object to merge gene names
sce <- readRDS(paste(local_data_path, "sce_objects/19595X1_labeled.rds", sep = "/"))
common_genes <- intersect(rownames(tothill), rowData(sce)$Symbol)

tothill <- tothill[rownames(tothill) %in% common_genes,]

# Convert from log intensity
tothill <- 2^tothill

tothill <- cbind(rownames(tothill), tothill)
colnames(tothill)[1] <- "Gene"

# Save expression file
outfile <- paste(local_data_path, "deconvolution_input",
                 "bulk_data_tothill.tsv", sep= "/")
write.table(tothill, outfile, sep = "\t", row.names = F, quote = F)

# Make survival file
pheno <- pheno[hgsoc_samples, ]
pheno <- subset(pheno, select = c("alt_sample_name", "grade", "age_at_initial_pathologic_diagnosis",
                                  "days_to_tumor_recurrence", "recurrence_status", "days_to_death",
                                  "vital_status"))

outfile <- paste(local_data_path, "Tothill",
                 "tothill_survival.tsv", sep = "/")
write.table(pheno, outfile, row.names = F, sep = "\t", quote = F)
