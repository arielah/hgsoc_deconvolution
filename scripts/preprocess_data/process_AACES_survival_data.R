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