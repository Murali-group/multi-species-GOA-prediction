# multi-species-GOA-prediction

This repository contains the main LOSO validation pipeline, algorithms, and plotting functions for the paper 
[Accurate and Efficient Gene Function Prediction using a Multi-Bacterial Network](http://dx.doi.org/10.1093/bioinformatics/btaa885).

## Download Datasets
The networks and GO term annotations for the 200 bacterial species with the most EXPC and COMP annotations are available here: https://doi.org/10.5281/zenodo.4280990

## Setup Environment

### Datasets
Run the following command (from the base directory) to download the datasets and setup the environment in which the datasets are used. 

```
bash setup.sh
```

### Python
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

## Run the Pipeline

### Config file
This script requires a config ([YAML](https://yaml.org/)) file specific to the [annotation_prediction pipeline](https://github.com/Murali-group/annotation_prediction), with a few additional options set.

The config files contain the documentation for the different options. See `config_files/expc/expc-loso-net-combinations-bp.yaml` for an example.
To change which algorithms are run, change the `should_run` flag to `True/False`. 

> Note that the annotation_prediction code was added as a [git-subrepo](https://github.com/ingydotnet/git-subrepo), so to make changes to that code, please commit them to that repository directly, and then pull them with `git subrepo pull src/annotation_prediction/`, or follow the suggestions in the git-subrepo documentation.

### Generate Predictions
To generate predictions, the following command extracts the core network from the network files and generates predictions.

```
python make_predictions.py --config config-files/expc-core-bp.yaml
```

If you do not need to specify the core and target species meaning you have already built a network and wish to make predictions for all nodes in the network, you can use the standard `run_eval_algs.py` script in the [annotation_prediction pipeline](https://github.com/Murali-group/annotation_prediction).

```
python src/annotation_prediction/run_eval_algs.py --config <config_file>
```

> TODO make an example config_file and datasets.

### Run the Evaluations
To run the leave-one-species-out (LOSO) validation and the COMP and ELEC validations in the paper, use the following commands:

EXPC LOSO:
```
python -u run_experiments.py --config config-files/expc/expc-loso-net-combinations-bp.yaml
```

EXPC eval COMP:
```
python -u run_experiments.py --config config-files/expc-eval-comp/expc-eval-comp-net-comb-bp.yaml --skip-core
```

EXPC eval ELEC:
```
python -u run_experiments.py --config config-files/expc-eval-elec/expc-eval-elec-net-comb-bp.yaml --skip-core
```

### Use Screen or Job Scheduler
The script `start_jobs.py` allows you to run multiple evaluations/algorithms in parallel using either a screen sessions on a single machine, or by submitting jobs to a job scheduler, and also automatically writes a job-specific config file and log file. This script ignores the `should_run` flag set in the config file, allowing you to specify which algorithms should be run directly.

If using an anaconda environment, include conda's `activate` command in the `--python` option to load the correct environment when running the jobs.

Example call:
```
python start_jobs.py \
  --config config-files/expc-eval-comp/expc-eval-comp-net-comb-bp.yaml  \
  --alg fastsinksource --alg genemania \
  --alg birgrank --alg sinksource --alg localplus   \
  --job-per-dataset --job-per-param \
  --qsub --cores 3 \
  --script-to-run run_experiments.py \
  --python="source /data/jeff-law/tools/anaconda3/bin/activate fungcat; python" \
  --pass-to-script --skip-core
```

### Plot the results
TODO.

## Cite
If you use FastSinkSource or other methods in this package, please cite:

Jeffrey N. Law, Shiv D. Kale, and T. M. Murali. [Accurate and Efficient Gene Function Prediction using a Multi-Bacterial Network](http://dx.doi.org/10.1093/bioinformatics/btaa885), _Bioinformatics_ (2020). doi.org/10.1093/bioinformatics/btaa885

