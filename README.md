# Cell type deconvolution analysis of HGSOC subtypes

## Overview

Ever since the first TCGA paper on the subject came out in 2011, high-grade serous ovarian cancer (HGSOC) has been described as having four transcriptional subtypes: differentiated, immmunoreactive, mesenchymal, and proliferative. These subtypes seemed to have differential survival prognoses. However, the precise number of subtypes and the most appropriate way to describe them has been questioned. We suspected that the transcriptional subtypes represented different cell type compositions in the tumor, meaning that subtypes had more to do with stromal admixture than the actual cancer cells themselves.

To see if there are differences in the cell type composition of the HGSOC subtypes, we estimated cell type profiles for several large bulk transcriptomic datasets using a deconvolution method called BayesPrism. We chose to use BayesPrism based on results from our [pilot study](https://github.com/greenelab/deconvolution_pilot). As an example, we found the following cell type proportions in the RNA-seq data from TCGA:

![](https://github.com/arielah/hgsoc_deconvolution/blob/prepublic/figures/readme_figure.png)

We then compared cell type info across subtype assignments at k = 3 and k = 4. We looked for a relationship between cell type proportions and survival/other clinical factors. Finally, we ran nonnegative matrix factorization (NMF) on the cancer fraction, to see if there were distinct cancer expression profiles across subtypes.

## Data

Our main reference profile, called HGSOC Penn/Utah, is made of single-cell RNA-sequencing data originating from 15 HGSOC tumors that were collected at two sites, the University of Pennsylvania (Penn) and the University of Utah (Utah). The raw data (FASTQ files) for these tumors are available at dbGaP: https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs002262.v3.p3 (note: the dbGaP submission process is still underway. In the meantime refer to https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs002262.v1.p1 for the first few samples.) The processed data (gene count matrices) are available at GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE217517.

Our validation reference profile, called Vázquez-Garcia, is a subset of the single-cell RNA-sequencing data from [this paper](https://www.nature.com/articles/s41586-022-05496-1). The full dataset is available in GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE180661.

The bulk data and associated covariates are all publicly available: we recommend using the Bioconductor package `curatedOvarianData` to obtain them, it's faster than finding them on GEO.

## Repository layout

A description of all scripts in this repository is available [here](https://github.com/arielah/hgsoc_deconvolution/tree/prepublic/scripts).

## Citation

%Placeholder until preprint is approved on biorXiv
