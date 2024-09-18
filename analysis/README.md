# HBV Genotype Analysis Pipeline

This pipeline is used to analyze the genotype of the Hepatitis B Virus (HBV). It consists of three main scripts: `anc_reconstruct.r`, `get_allele_homoplasy.r`, and `lasso_analysis_nofused.r`.

## anc_reconstruct.r

This script is the first step in the pipeline. It is used to reconstruct the ancestral state of the HBV genotype. The output of this script is a file containing the ancestral restructure data, which is used as input for the next script in the pipeline.

parameters in the file that you might want to change:


## get_allele_homoplasy.r

This script takes the ancestral restructure data from `anc_reconstruct.r` and analyzes allele homoplasy. It uses several R packages including `tidyverse`, `ape`, `phangorn`, `foreach`, and `doParallel`. The output of this script is a file containing the allele homoplasy data, which is used as input for the final script in the pipeline.

## lasso_analysis_nofused.r

This script takes the allele homoplasy data from `get_allele_homoplasy.r` and performs a LASSO analysis. The LASSO (Least Absolute Shrinkage and Selection Operator) is a regression analysis method that performs both variable selection and regularization in order to enhance the prediction accuracy and interpretability of the statistical model it produces.

## Pipeline Usage

To run the pipeline, use the following commands in the terminal:

```bash
Rscript anc_reconstruct.r 
Rscript get_allele_homoplasy.r <genotype> <anc_restruct_file> <homoplasy_file>
Rscript lasso_analysis_nofused.r <lasso_output_file>
```