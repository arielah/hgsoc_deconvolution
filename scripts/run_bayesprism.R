suppressPackageStartupMessages({
    library(data.table)
    library(SingleCellExperiment)
    library(dplyr)
    library(BayesPrism)
    library(yaml)
})

bulk_type <- snakemake@wildcards[['bulk_type']]
sc_type <- snakemake@wildcards[['sc_type']]
params <- read_yaml("../config.yml")
data_path <- params$data_path
local_data_path <- params$local_data_path

outdir <- paste(local_data_path, "/deconvolution_output/", bulk_type, "_", sc_type, sep = "")
print(outdir)

# Get bulk counts matrix
bulk_matrix <- fread(paste(local_data_path, "/deconvolution_input/",
                           "bulk_data_", bulk_type, ".tsv", sep = ""),
                     header = TRUE)
genes <- bulk_matrix$Gene; bulk_matrix$Gene <- NULL
sample_names <- colnames(bulk_matrix)
bulk_matrix <- t(bulk_matrix)
colnames(bulk_matrix) <- genes

# Get single cell data

scefile <- paste(local_data_path, "/deconvolution_input/",
                 "single_cell_data_", sc_type,".rds", sep = "")
sce <- readRDS(scefile)
rownames(sce) <- rowData(sce)$Symbol

gene_duplicates <- which(duplicated(rowData(sce)$Symbol))
sce <- sce[-gene_duplicates, ]

# Subset down to genes in bulk matrix
sce <- sce[which(rownames(sce) %in% genes),]
single_cell_matrix <- t(as.matrix(assay(sce)))

# Get cell type and cell state labels
cell_types <- sce$cellType
cell_states <- sce$cellTypeGranular


# Cleanup genes
sc.dat.filtered <- cleanup.genes(input = single_cell_matrix, input.type = "count.matrix",
                                 species="hs", gene.group = c("Rb","Mrp","other_Rb","chrM","MALAT1"), exp.cells = 5)

# Filter to protein coding genes for faster computation
sc.dat.filtered.pc <- select.gene.type(sc.dat.filtered, gene.type = "protein_coding")

rm(single_cell_matrix, sc.dat.filtered, sce, genes); gc()

# Make prism object
myPrism <- new.prism(reference = sc.dat.filtered.pc,
                     mixture = bulk_matrix,
                     input.type = "count.matrix",
                     cell.type.labels = cell_types,
                     cell.state.labels = cell_states,
                     key="Epithelial cells")

rm(sc.dat.filtered.pc); gc()

bp.res <- run.prism(prism = myPrism, n.cores=6)
theta <- get.fraction(bp = bp.res,
                      which.theta = "final",
                      state.or.type = "type")

# Save BayesPrism object for later perusal
object_file <- paste(outdir, "_bayesprism_results_full.rds", sep = "")
saveRDS(bp.res, file = object_file)

# Format text version of proportion estimates
theta <- as.data.frame(t(theta))
theta <- cbind(rownames(theta), theta)
colnames(theta) <- c("cell_type", sample_names)

# Save proportion estimates
text_file <- paste(outdir, "_bayesprism_results.tsv", sep = "")
write.table(theta, file = text_file, sep = "\t", quote = F, row.names = F)
