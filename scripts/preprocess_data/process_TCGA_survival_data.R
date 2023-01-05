suppressPackageStartupMessages({
    library(data.table)
    library(SingleCellExperiment)
    library(dplyr)
    library(yaml)
})

params <- read_yaml("../../config.yml")
data_path <- params$data_path
local_data_path <- params$local_data_path

# I downloaded an xlsv file from GDC and then converted each
# tab of the spreadsheet into a tsv file.
tcga <- fread(paste(local_data_path, "TCGA",
                    "TCGA-CDR-SupplementalTableS1/TCGA-CDR-Table\ 1.tsv",
                    sep = "/"))
tcga$V1 <- NULL

table(tcga$type)

# Subset to only ovarian tumors
tcga <- subset(tcga, tcga$type=="OV")
table(tcga$histological_grade)

# Subset to only high-grade tumors
tcga <- subset(tcga, tcga$histological_grade %in% c("G2", "G3", "G4"))

outfile <- paste(local_data_path, "TCGA",
                 "TCGA_OV_survival.tsv", sep = "/")
write.table(tcga, outfile, row.names = F, sep = "\t", quote = F)
