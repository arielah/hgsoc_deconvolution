suppressPackageStartupMessages({
    library(data.table)
    library(SingleCellExperiment)
    library(dplyr)
    library(yaml)
})

params <- read_yaml("../../config.yml")
data_path <- params$data_path
local_data_path <- params$local_data_path

# Load AACES data, copied from hgsc_characterization repo
aaces <- fread(paste(local_data_path, "AACES", 
                     "salmon_raw_counts_for_way_pipeline.tsv", sep = "/"))
aaces[1:5,1:5]
setnames(aaces, "V1", "gene_name")

# Load a small SCE object to merge gene names
sce <- readRDS(paste(local_data_path, "sce_objects/19595X1_labeled.rds", sep = "/"))
common_genes <- intersect(aaces$gene_name, rowData(sce)$Symbol)

aaces <- subset(aaces, aaces$gene_name %in% rowData(sce)$Symbol)
setnames(aaces, "gene_name", "Gene")

outfile <- paste(local_data_path, "deconvolution_input",
                 "bulk_data_AACES.tsv", sep= "/")
write.table(aaces, outfile, sep = "\t", row.names = F, quote = F)
