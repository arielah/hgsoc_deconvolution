# Scripts: cell type deconvolution analysis of HGSOC subtypes

## single_cell_annotation
This directory contains R scripts for manually annotating our single cell data using a combination of unsupervised clustering and [CellTypist](https://www.celltypist.org/). 

## preprocess_data
This directory contains for preprocessing single-cell and bulk data for deconvolution, as well as preparing clinical covariate data for later evaluations. It includes scripts to:

- Merge all HGSOC Penn/Utah datasets into one file
- Randomly downsample the Vázquez-Garcia data into a workable size for BayesPrism
- Process TCGA and Tothill expression and survival data

## Snakefile
A workflow definition file for Snakemake. It defines the order in which the main deconvolution scripts should be run.

## run_bayesprism.R
An R script to run BayesPrism package on each combination of bulk dataset and single-cell reference profile.

## evaluation
This directory contains R scripts for various evaluation tasks such as:

- Comparing deconvolution results for TCGA RNA-seq and TCGA microarray datasets
- Comparing cell composition data to optimal/suboptimal surgical outcomes
- Selecting a good value of k for the cancer fraction embedding (NMF) analysis
- Generating embeddings from the cancer fraction and visualizing the results
- Performing Cox proportional hazards tests on survival based on fibroblast and immune proportions.

## figure_generation
This directory contains R scripts for generation figures for the final paper.
