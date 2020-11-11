# annotation-prediction
Pipeline for running and evaluating algorithms for gene/protein annotation prediction, 
e.g., to the Gene Ontology (GO) and to the Human Phenotype Ontology(HPO).


## Installation

- Required Python packages: `networkx`, `numpy`, `scipy`, `pandas`, `sklearn`, `pyyaml`, `tqdm`, `rpy2`
- Required R packages: PPROC

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

## Usage 

### Cross Validation
The relevant options are below. See `python run_eval_algs.py --help` for more details.
  - `cross_validation_folds`
    - Number of folds to use for cross validation. Specifying this parameter will also run CV
  - `cv_seed`
    - Can be used to specify the seed to use when generating the CV splits. 
    
Example:
```
python run_eval_algs.py  --config config.yaml --cross-validation-folds 5 --only-eval
```

#### Plot
After CV has finished, to visualize the results, use the `plot.py` script. For example:
```
python plot.py --config config.yaml --box --measure fmax
```

## Cite
If you use FastSinkSource or other methods in this package, please cite:

Jeffrey Law, Shiv D. Kale, and T. M. Murali. [Accurate and Efficient Gene Function Prediction using a Multi-Bacterial Network](https://doi.org/10.1101/646687), _bioRxiv_ (2020). doi.org/10.1101/646687
