# multi-species-GOA-prediction

This repository contains the main LOSO validation pipeline, algorithms, and plotting functions for the paper 
[Accurate and Efficient Gene Function Prediction using a Multi-Bacterial Network](http://dx.doi.org/10.1093/bioinformatics/btaa885).

## Installation
These scripts require Python 3 due to the use of obonet to build the GO DAG.

- Required Python packages: `networkx`, `numpy`, `scipy`, `pandas`, `sklearn`, `obonet`, `pyyaml`, `tqdm`, `rpy2`
- Required R packages: PPROC

We recommend using [Anaconda](https://www.anaconda.com/) for Python, especially to access the needed R packages. 
To setup your environment, use the following commands:

```
conda create -n fastsinksource python=3.7 r=3.6
conda activate fastsinksource
pip install -r requirements.txt
```
To install the R packages:
```
R -e "install.packages('https://cran.r-project.org/src/contrib/PRROC_1.3.1.tar.gz', type = 'source')"
conda install -c bioconda bioconductor-clusterprofiler
```

## Download Datasets
The networks and GO term annotations for the 200 bacterial species with the most EXPC and COMP annotations are available here: http://bioinformatics.cs.vt.edu/~jeffl/supplements/2020-fastsinksource/

## Run the Pipeline
This script requires a config file specific to the [annotation_prediction pipeline](https://github.com/Murali-group/annotation-prediction), with a few additional options set.

> Note that the annotation_prediction code was added as a [git-subrepo](https://github.com/ingydotnet/git-subrepo), so to make changes to that code, please commit them to that repository directly, and then pull them with `git subrepo pull src/annotation_prediction/`, or follow the suggestions in the git-subrepo documentation.


## Cite
If you use FastSinkSource or other methods in this package, please cite:

Jeffrey N. Law, Shiv D. Kale, and T. M. Murali. [Accurate and Efficient Gene Function Prediction using a Multi-Bacterial Network](http://dx.doi.org/10.1093/bioinformatics/btaa885), _Bioinformatics_ (2020). doi.org/10.1093/bioinformatics/btaa885

