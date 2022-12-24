suppressPackageStartupMessages({
    library(data.table)
    library(stringr)
    library(dplyr)
    library(SingleCellExperiment)
    library(yaml)
})

params <- read_yaml("../../config.yml")
data_path <- params$data_path
local_data_path <- params$local_data_path

# Load tcga data, downloaded from GDC
tcga_rna <- fread(paste(local_data_path, "TCGA",
                        "EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2-v2.geneExp.tsv",
                        sep = "/"))

# Trim down names to patient only
old_names <- colnames(tcga_rna)
rna_samples <- str_extract(old_names, "TCGA-\\w\\w-\\w\\w\\w\\w")
rna_samples[1] <- "gene_id"
head(rna_samples)

# Load in survival data to filter down to just HGSOC cases
tcga_survival <- fread(paste(local_data_path, "TCGA", "TCGA_OV_survival.tsv", sep = "/"))
columns_to_keep <- c(1, which(rna_samples %in% tcga_survival$bcr_patient_barcode))
tcga_rna <- tcga_rna %>% select(columns_to_keep)

# Relabel columns with patient and sample id
new_names <- str_extract(colnames(tcga_rna), "TCGA-\\w\\w-\\w\\w\\w\\w-\\w\\w\\w")
new_names[1] <- "gene_id"
colnames(tcga_rna) <- new_names

# The gene names in the data directly from TCGA are formatted like "ABCF1|23",
# probably a mapping to an index I can't find on their website. Also, 29 of
# them have question marks for names, e.g. "?|100133144". This separates out
# only the HGNC symbols and removes the ? genes.
gene_names <- strsplit(tcga_rna$gene_id, split = "\\|")
gene_names <- sapply(gene_names, "[[", 1)
tcga_rna$gene_id <- gene_names

tcga_rna <- subset(tcga_rna, tcga_rna$gene_id != "?")
tcga_rna[1:5, 1:5]

# Duplicate gene to remove
gene_duplicates <- which(duplicated(tcga_rna$gene_id))
tcga_rna <- tcga_rna[-gene_duplicates, ]

# Filter to genes in master list (in both TCGA and single-cell annotations)
master_gene_list <- fread(paste(local_data_path, "TCGA",
                                "master_gene_list.tsv", sep = "/"))
tcga_rna <- subset(tcga_rna, tcga_rna$gene_id %in% master_gene_list$bulk_name)

# Switch TCGA gene names to single-cell names for clarity
setnames(tcga_rna, "gene_id", "bulk_name")
tcga_rna <- inner_join(master_gene_list, tcga_rna)
tcga_rna$bulk_name <- NULL
setnames(tcga_rna, "sc_name", "Gene")

# Remove genes with NAs, these break BayesPrism
nas_per_gene <- rowSums(is.na(tcga_rna))
tcga_rna <- tcga_rna[which(nas_per_gene == 0), ]

# Switch any negative values to 0, smallest value is approx. -1
tcga_rna[tcga_rna < 0] <- 0

outfile <- paste(local_data_path, "deconvolution_input",
                 "bulk_data_TCGA.tsv", sep= "/")
write.table(tcga_rna, outfile, sep = "\t", row.names = F, quote = F)
