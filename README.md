# multi-species-GOA-prediction

This repository contains the main LOSO validation pipeline, algorithms, and plotting functions for the paper 
[Accurate and Efficient Gene Function Prediction using a Multi-Bacterial Network](http://dx.doi.org/10.1093/bioinformatics/btaa885).

## Download Datasets
The networks and GO term annotations for the 200 bacterial species with the most EXPC and COMP annotations are available here: https://zenodo.org/record/4280990

## Setup Environment

### Datasets
More information about the input datasets is available in the [inputs README](https://github.com/Murali-group/multi-species-GOA-prediction/blob/master/inputs/README.md).

Run the `setup.sh` script (from the base directory) to download the datasets and setup the environment in which the datasets are used. 

```
bash setup.sh
```

### Python
Run the `setup_conda_venv.sh` script to install the required packages. The commands are also listed below.

```
bash setup_conda_venv.sh
```

- Required Python packages: `networkx`, `numpy`, `scipy`, `pandas`, `sklearn`, `obonet`, `pyyaml`, `tqdm`, `rpy2`
- Required R packages: PPROC

We recommend using [Anaconda](https://www.anaconda.com/) for Python, especially to access the needed R packages. 
To setup your environment, use the following commands:

```
conda create -y --name mult-sp-pred python=3.7 r=3.6 
conda activate mult-sp-pred
pip install -r requirements.txt
```
To install the R packages:
```
R -e "install.packages('https://cran.r-project.org/src/contrib/PRROC_1.3.1.tar.gz', type = 'source')"
```

## Run the Pipeline

### Config file
This script requires a config ([YAML](https://yaml.org/)) file specific to the [annotation_prediction pipeline](https://github.com/Murali-group/annotation_prediction), with a few additional options set.

The config files contain the documentation for the different options. See `config_files/expc/expc-loso-net-combinations-bp.yaml` for an example.
To change which algorithms are run, change the `should_run` flag to `True/False`. 

> Note that the annotation_prediction code was added as a [git-subrepo](https://github.com/ingydotnet/git-subrepo), so to make changes to that code, please commit them to that repository directly, and then pull them with `git subrepo pull src/annotation_prediction/`, or follow the suggestions in the git-subrepo documentation.

### Generate Predictions
To generate predictions, `make_predictions.py` extracts the core network from the network files and generates predictions. The default number of predictions stored is 10. To write more, use either the `--num-pred-to-write` or `--factor-pred-to-write` options (see `python make_predictions.py --help`). To write the prediction scores of all nodes:

```
python make_predictions.py --config config-files/expc-core-bp.yaml --num-pred-to-write -1
```

If you do not need to specify the core and target species, meaning you have already built a network and wish to make predictions for all nodes in the network, you can use the standard `run_eval_algs.py` script in the [annotation_prediction pipeline](https://github.com/Murali-group/annotation_prediction).

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

Since AsyRW is run on only the last combination, we need to use a config file with only that combination. For example:
```
python -u run_experiments.py --config config-files/expc/expc-loso-ssnNbrs-bp.yaml
```

#### Run Using Screen or HPC Job Scheduler
Since we are evaluating many different algorithms on large datasets, we often need more compute power than a single machine. The script `start_jobs.py` allows you to run multiple evaluations/algorithms in parallel using either a screen sessions on a single machine, or by submitting jobs to a job scheduler in a high-performance computing cluster, and also automatically writes a job-specific config file and log file. This script ignores the `should_run` flag set in the config file, allowing you to specify which algorithms should be run directly.

If using an anaconda environment, include conda's `activate` command in the `--python` option to load the correct environment when running the jobs.

Example calls:

EXPC LOSO: 
```
python start_jobs.py \
  --config config-files/expc/expc-loso-net-combinations-bp.yaml \
  --alg fastsinksource --alg genemania \
  --alg birgrank --alg sinksource --alg localplus \
  --job-per-dataset --job-per-param \
  --qsub --cores 3 \
  --script-to-run run_experiments.py   \
  --python="source /data/jeff-law/tools/anaconda3/bin/activate mult-sp-pred; python"
```

EXPC eval COMP:
```
python start_jobs.py \
  --config config-files/expc-eval-comp/expc-eval-comp-net-comb-bp.yaml  \
  --alg fastsinksource --alg genemania \
  --alg birgrank --alg sinksource --alg localplus   \
  --job-per-dataset --job-per-param \
  --qsub --cores 3 \
  --script-to-run run_experiments.py \
  --python="source /data/jeff-law/tools/anaconda3/bin/activate mult-sp-pred; python" \
  --pass-to-script --skip-core
```

Note that this script is setup to submit jobs to the small, 6 node baobab cluster in the bioinformatics group at Virginia Tech, but it should be easy enough to modify for different compute clusters.

Once the jobs are submitted, you'll notice that 

### Plot the results
To recreate Figure 3a-c of the paper, use this command to make the fmax summary plot for EXPC LOSO:
```
python src/annotation_prediction/plot.py \
  --config config-files/expc/expc-loso-net-combinations-bp.yaml \
  --line --ci 95 --horiz --for-paper --local-line \
  --alg fastsinksource --alg genemania --alg birgrank \
  --alg sinksource --alg async_rw --alg localplus
```

Then run the same command with the different config files for EXPC eval COMP and EXPC eval ELEC. I used inkscape to combine each of the figures after that. Note that to get the y-axis the same for b anc c, I used `--ylim 0.58 0.85`.


## Cite
If you use FastSinkSource or other methods in this package, please cite:

Jeffrey N. Law, Shiv D. Kale, and T. M. Murali. [Accurate and Efficient Gene Function Prediction using a Multi-Bacterial Network](http://dx.doi.org/10.1093/bioinformatics/btaa885), _Bioinformatics_ (2020). doi.org/10.1093/bioinformatics/btaa885

