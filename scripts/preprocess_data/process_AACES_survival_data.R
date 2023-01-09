suppressPackageStartupMessages({
    library(data.table)
    library(SingleCellExperiment)
    library(dplyr)
    library(yaml)
})

params <- read_yaml("../../config.yml")
data_path <- params$data_path
local_data_path <- params$local_data_path

aaces <- fread(paste(local_data_path, "AACES",
                     "tissuegrant_epidata_08102022.csv", sep = "/"))

id_map <- fread(paste(local_data_path, "AACES",
                      "main_metadata_table.tsv", sep = "/"))

aaces <- left_join(aaces, id_map)
aaces <- subset(aaces, select = c("ID", "study", "refage", "diagyear", 
                                  "vitalstatus", "survival_days", "neoadj", 
                                  "dblkstat", "seer_stage", "tissue_source",
                                  "ClusterK4_kmeans_TCGA_names", 
                                  "external_HGSCsubtype_estimate"))

aaces$ID <- paste("Sample", aaces$ID, sep = "_")

outfile <- paste(local_data_path, "AACES", "AACES_survival.tsv", sep = "/")
write.table(aaces, file = outfile, sep = "\t", row.names = F, quote = F)
