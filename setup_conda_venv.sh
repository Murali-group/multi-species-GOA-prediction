#!/bin/bash

# Set-up Anaconda virtual environment
echo "Setting up Anaconda Python virtual environment..."

conda create -y --name mult-sp-pred python=3.7 r=3.6 
conda activate mult-sp-pred
pip install -r requirements.txt

# Install the PRROC package for computing area under PR curve
# TODO: Write the PRROC AUC function without using rpy2
R -e "install.packages('https://cran.r-project.org/src/contrib/PRROC_1.3.1.tar.gz', type = 'source')"
