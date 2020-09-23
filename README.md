# multi-species-GOA-prediction

This repository contains the main LOSO validation pipeline, algorithms, and plotting functions for the paper 
[Accurate and Efficient Gene Function Prediction using a Multi-Bacterial Network](https://doi.org/10.1101/646687).

## Installation
These scripts requires Python 3 due to the use of obonet to build the GO DAG.

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


#### Plot
After CV has finished, to visualize the results, use the `plot.py` script. For example:
```
python plot.py --config config.yaml --box --measure fmax
```

## Cite
If you use FastSinkSource or other methods in this package, please cite:

Jeffrey N. Law, Shiv D. Kale, and T. M. Murali. [Accurate and Efficient Gene Function Prediction using a Multi-Bacterial Network](https://doi.org/10.1101/646687), _bioRxiv_ (2020). doi.org/10.1101/646687

