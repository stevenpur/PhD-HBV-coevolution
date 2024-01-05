#!/bin/bash
# set the env needed for the diong the analysis
# check if conda env gwas already exists. If not, create it
echo "checking if conda env gwas exists..."
    if [ -z "$(conda env list | grep hbv_covar)" ]
    then
        # create the env and install the required packages]
        echo "conda env gwas does not exist, creating it..."
        conda env create -f requirements.yml
    fi
    # activate the env
    echo "activating conda env gwas..."
    conda activate hbv_covar