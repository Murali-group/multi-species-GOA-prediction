
# This script sets up the core/target species and networks and performs the coresponding evaluations

import argparse
import os, sys
import numpy as np
from scipy import sparse as sp
from tqdm import tqdm
import yaml
from collections import defaultdict
#import pandas as pd

# for some reason, the current directory isn't always in the path. this fixes that
if "" not in sys.path:
    sys.path.insert(0,"")
import run_experiments as run
from src.annotation_prediction.src import main as run_eval_algs
from src.annotation_prediction.src.evaluate import eval_leave_one_species_out as eval_loso
from src.annotation_prediction.src import setup_sparse_networks as setup
#import src.annotation_prediction.src.go_term_prediction_examples.go_term_prediction_examples as go_examples
#import src.annotation_prediction.src.algorithms.aptrank_birgrank.run_birgrank as run_birgrank
import src.annotation_prediction.src.utils.config_utils as config_utils
import src.annotation_prediction.src.utils.file_utils as utils
import src.annotation_prediction.src.algorithms.alg_utils as alg_utils


def main(config_map, **kwargs):
    input_settings, input_dir, output_dir, alg_settings, kwargs \
        = config_utils.setup_config_variables(config_map, **kwargs)

    for dataset in input_settings['datasets']:
        kwargs = run.apply_dataset_settings_to_kwargs(dataset, **kwargs)
        # get/set the minimum number of annotations to make predictions for a given term
        kwargs['num_ann_cutoff'] = kwargs.get('num_ann_cutoff', 10) 
        print({key: val for key,val in kwargs.items() \
                if key not in ['species_to_uniprot_idx', 'alg_taxon_terms_to_skip']})
        #print(kwargs)
        net_obj, ann_obj, eval_ann_obj = run_eval_algs.setup_dataset(dataset, input_dir, **kwargs) 
        # if there are no annotations, then skip this dataset
        if len(ann_obj.terms) == 0:
            print("No terms found. Skipping this dataset")
            continue

        # add the taxon file paths for this dataset to kwargs
        for arg in ['taxon_prot_file', 'core_taxons_file', 'target_taxons_file']:
            kwargs[arg] = "%s/%s" % (input_dir, dataset[arg]) if arg in dataset else None
        species_to_uniprot_idx = eval_loso.get_uniprot_species(kwargs['taxon_prot_file'], ann_obj)
        # set this in kwargs to use it in all functions
        kwargs['species_to_uniprot_idx'] = species_to_uniprot_idx

        # the outputs will follow this structure:
        # outputs/<net_version>/<exp_name>/<alg_name>/output_files
        out_dir = "%s/%s/%s/" % (output_dir, dataset['net_version'], dataset['exp_name'])
        alg_runners = run_eval_algs.setup_runners(alg_settings, net_obj, ann_obj, out_dir, **kwargs)

        # Generate predictions using the specified settings
        setup_networks_and_run(alg_runners, net_obj, ann_obj, eval_ann_obj, **kwargs)

        # if specified, write the SWSN combined network to a file
        save_net = dataset['net_settings'].get('save_net', None) if 'net_settings' in dataset else None
        if net_obj.weight_swsn is True and save_net is not None:
            for run_obj in alg_runners:
                out_file = "%s/%s/%s" % (input_dir, dataset['net_version'], save_net)
                # the SWSN network is part of the runner object. Need to organize that better
                run_obj.net_obj.save_net(out_file)


def setup_networks_and_run(alg_runners, net_obj, ann_obj, eval_ann_obj, **kwargs):
    # species_names is a dictionary of taxonomy ID to name, from the target_taxons_file
    # target_taxons is the list of taxons for which we want to make predictions
    species_names, target_taxons = eval_loso.get_selected_species(
        kwargs['species_to_uniprot_idx'], kwargs['target_taxons_file'], kwargs['taxons'])
    # keep track of the original to use it later
    orig_net_obj = net_obj
    # limit the networks to the given taxons 
    # TODO this is now the only way to run this script
    #if kwargs.get('core_taxons_file'):
    # read in the specified taxons from the file
    _, core_taxons = eval_loso.get_selected_species(
        kwargs['species_to_uniprot_idx'], kwargs['core_taxons_file'])
    core_taxon_prots = run.get_taxon_prots(
        len(net_obj.nodes), core_taxons, kwargs['species_to_uniprot_idx'])
    core_net_obj, core_ann_obj = run.limit_to_taxons(core_taxon_prots, net_obj=net_obj, ann_obj=ann_obj, **kwargs)

    # limit the original net to the core and target taxons
    target_taxon_prots = run.get_taxon_prots(
        len(orig_net_obj.nodes), target_taxons, kwargs['species_to_uniprot_idx'])
    # limit the SSN to the prots in the core+target
    core_target_prots = core_taxon_prots + target_taxon_prots
    orig_net_obj = run.limit_to_taxons(core_target_prots, net_obj=orig_net_obj, **kwargs)
    run_alg_target_sp(
        target_taxons, target_taxon_prots, core_taxon_prots,
        alg_runners, orig_net_obj, core_net_obj,
        core_ann_obj, **kwargs)


def run_alg_target_sp(
        target_taxons, target_taxon_prots, core_taxon_prots,
        alg_runners, orig_net_obj, core_net_obj,
        core_ann_obj, **kwargs):
    """
    Run the algorithms using the annotations of the core species only,
    and make predictions for all of the target species. 
    """
    # the annotations were alread limited to the core species
    core_ann_mat, terms = core_ann_obj.ann_matrix, core_ann_obj.terms
    num_train_pos_per_term = (core_ann_mat > 0).sum(axis=1)
    num_ann_cutoff = kwargs['num_ann_cutoff']
    terms_to_run = []
    for i, term in enumerate(terms):
        if num_train_pos_per_term[i] < num_ann_cutoff:
            continue
        terms_to_run.append(term)
    print("\t%d with >= %d annotations among the %d core terms " % (
        len(terms_to_run), num_ann_cutoff, len(terms)))

    if kwargs.get('ssnNbrs'):
        # For all other target species, precompute the scores and then transfer (ssnLocal),
        compute_core_scores_then_transfer(
            target_taxons, target_taxon_prots, core_taxon_prots,
            alg_runners, orig_net_obj, core_net_obj, core_ann_obj,
            terms_to_run, **kwargs)
    else:
        # connect all the target networks to the core, then compute the scores
        connect_target_species_run_eval(
            target_taxons, target_taxon_prots, core_taxon_prots,
            alg_runners, orig_net_obj, core_net_obj, core_ann_obj,
            terms_to_run, **kwargs)

    # now write the prediction scores to file
    for run_obj in alg_runners:
        num_pred_to_write = kwargs.get('num_pred_to_write',10)
        num_pred_to_write = 10 if num_pred_to_write is None else num_pred_to_write
        if kwargs.get('factor_pred_to_write') is not None:
            # make a dictionary with the # ann*factor for each term
            num_pred_to_write = {} 
            for i in range(run_obj.ann_matrix.shape[0]):
                y = run_obj.ann_matrix[i,:]
                positives = (y > 0).nonzero()[1]
                num_pred_to_write[run_obj.terms[i]] = len(positives) * kwargs['factor_pred_to_write']
        if num_pred_to_write != 0:
            out_file = "%s.txt" % (run_obj.out_pref)
            utils.checkDir(os.path.dirname(out_file)) 
            alg_utils.write_output(
                run_obj.term_scores, run_obj.ann_obj.terms, run_obj.ann_obj.prots,
                out_file, num_pred_to_write=num_pred_to_write)

    #eval_loso.write_stats_file(alg_runners, params_results)
    #print(params_results)
    

def print_net_stats(W):
    # count the nodes with at least one edge
    num_nodes = np.count_nonzero(W.sum(axis=0))
    # since the network is symmetric, the number of undirected edges is the number of entries divided by 2
    num_edges = (len(W.data) / 2)
    print("\t%s nodes, %s edges" % (num_nodes, num_edges))
    return num_nodes, num_edges


def compute_core_scores_then_transfer(
        target_taxons, target_taxon_prots, core_taxon_prots,
        alg_runners, orig_net_obj, core_net_obj, core_ann_obj,
        terms_to_run, **kwargs):
    params_results = defaultdict(int)
    print("\nComputing scores for the core species, then transferring to %d target species using ssnNbrs" % (len(target_taxons)))

    ssn = orig_net_obj.sparse_networks[0]

    # Setup/combine (e.g., SWSN) the core networks 
    train_mat = core_ann_obj.ann_matrix
    core_net_obj, _ = run.apply_net_comb_filters(
        train_mat, target_taxon_prots, core_net_obj, core_ann_obj, **kwargs)
    num_nodes, num_edges = print_net_stats(core_net_obj.W)
    params_results['core_num_nodes|edges'] = "%d|%d" % (num_nodes, num_edges)
    if orig_net_obj.weight_swsn:
        params_results['swsn_time'] += core_net_obj.swsn_time
        params_results['core_swsn_time'] = core_net_obj.swsn_time
        params_results['core_swsn_weights'] = core_net_obj.swsn_weights

    # Compute the scores for the core network.
    for i, run_obj in tqdm(enumerate(alg_runners)):
        alg = run_obj.name
        tqdm.write("Alg: %s" % (alg))

        run_obj.net_obj = core_net_obj
        run_obj.ann_obj = core_ann_obj

        # limit the run_obj to run on the terms for which there are annotations
        run_obj.terms_to_run = terms_to_run
        if kwargs.get('term') is not None:
            # if a subset of terms are specified, only run those
            run_obj.terms_to_run = kwargs['term']
        # limit the scores stored to the core 
        core_taxon_prot_idx = np.where(core_taxon_prots)[0]
        run_obj.target_prots = core_taxon_prot_idx
        # birgrank (and all 'gene-based' methods) uses the "nodes_to_run" option, so need to set that here as well
        if run_obj.get_alg_type() == 'gene-based':
            # set it so that it runs on all nodes in the core taxon. 
            run_obj.kwargs['nodes_to_run'] = core_taxon_prot_idx
            run_obj.kwargs['run_all_nodes'] = True

        # run each method
        run_obj.setupInputs()
        run_obj.run()
        run_obj.setupOutputs()
        if kwargs.get('verbose'):
            utils.print_memory_usage()

        # now that we have the scores, transfer them using ssnT and Local.
        # then limit it to just the edges from non-target prots to target prots
        core_to_target_ssn, _ = run.get_core_to_target_ssn(ssn, target_taxon_prots)
        # and run Local with that SSN
        run_obj = run.run_ssn_local(run_obj, core_to_target_ssn)

        if alg == "genemania":
            # BUG: for GM, there are unreachable target nodes. Those should get a score of the value of k
            # but they got a score of 0. Since most of the nodes in the target have a negative score,
            # this put unreachable left-out negative examples at a much higher score.
            # This function fixes that bug by setting the unknown examples to k in the target species
            run_obj.term_scores = run.fix_gm_default(run_obj, target_taxon_prots)

    return params_results


def connect_target_species_run_eval(
        target_taxons, target_taxon_prots, core_taxon_prots,
        alg_runners, orig_net_obj, core_net_obj, core_ann_obj, 
        terms_to_run, **kwargs):
    params_results = defaultdict(int)
    print("\nConnecting %d target taxons to the core species network" % (len(target_taxons)))

    # Build the SSN from the core to target species
    # need to remove the inter-target-species edges to get only core -> target edges
    # the first net is the SSN
    ssn = orig_net_obj.sparse_networks[0] if orig_net_obj.multi_net else orig_net_obj.W
    ssnC = core_net_obj.sparse_networks[0] if core_net_obj.multi_net else core_net_obj.W
    # Get the edges from the core species to all target species
    ssnT, _ = run.get_core_to_target_ssn(ssn, target_taxon_prots)
    # also get the species-specific SSN edges
    # but not edges between target species
    for t in target_taxons:
        curr_target_taxon_prots = run.get_taxon_prots(ssn.shape[0], [t], kwargs['species_to_uniprot_idx'])
        # This function call isn't right because its getting all edges to/from the target prots
        # in the SSN, meaning it's keeping edges between target species
        #curr_ssnT = get_edges_of_nodes(ssn, curr_target_taxon_prots)
        t_diag = sp.diags(curr_target_taxon_prots)
        intra_target_edges = t_diag.dot(ssn).dot(t_diag)
        ssnT += intra_target_edges
    # now rebuild the net obj to update the "normalized_nets"
    if orig_net_obj.multi_net: 
        # add ssnC here at least so it can be used when weighting the networks.
        sparse_nets = [ssnC+ssnT] + orig_net_obj.sparse_networks[1:]
        net_obj = setup.Sparse_Networks(
            sparse_nets, orig_net_obj.nodes, net_names=orig_net_obj.net_names,
            weight_method=orig_net_obj.weight_method)
    else:
        net_obj = setup.Sparse_Networks(ssnC+ssnT, orig_net_obj.nodes)

    train_mat = core_ann_obj.ann_matrix
    # now apply the network combination filters to get the specified combinations of networks
    net_obj = run.apply_net_comb_filters(
        train_mat, target_taxon_prots,
        net_obj, core_ann_obj, **kwargs)
    num_nodes, num_edges = print_net_stats(net_obj.W)
    params_results['core_target_num_nodes|edges'] = "%d|%d" % (num_nodes, num_edges)
    if orig_net_obj.weight_swsn:
        params_results['swsn_time'] += net_obj.swsn_time
        params_results['core_target_swsn_time'] = net_obj.swsn_time
        params_results['core_target_swsn_weights'] = net_obj.swsn_weights

    # Compute the scores for the core network.
    for i, run_obj in tqdm(enumerate(alg_runners)):
        alg = run_obj.name
        tqdm.write("Alg: %s" % (alg))
        run_obj.net_obj = net_obj
        run_obj.ann_obj = core_ann_obj
        print("%d pos, %d neg" % (len((train_mat > 0).data), len((train_mat < 0).data)))

        # limit the run_obj to run on the terms for which there are annotations
        run_obj.terms_to_run = terms_to_run
        if kwargs.get('term') is not None:
            # if a subset of terms are specified, only run those
            run_obj.terms_to_run = kwargs['term']

        # run each method
        run_obj.setupInputs()
        run_obj.run()
        run_obj.setupOutputs()
        if kwargs.get('verbose'):
            utils.print_memory_usage()
        print("\n\t%d total gene-term pairs with scores" % len(run_obj.term_scores.data))
    return params_results


def store_standalone_network(net_obj, out_file, forced=False):
    """
    Store the SWSN combined network for the core (pre-normalization) that is used to make predictions 
    """
    pass


def setup_parser(parser=None):
    if parser is None:
        parser = argparse.ArgumentParser(
            description='Script for setting up various string experiments')

        parser.add_argument('--config', required=True,
            help='Configuration file')
        parser.add_argument('--taxon', '-T', dest="taxons", type=str, action='append',
                help="Specify the species taxonomy ID for which to evaluate. Multiple may be specified. Otherwise, all species will be used")
        #parser.add_argument('--alg', action="append", 
        #    help="Name of algorithm to run. May specify multiple. Default is whatever is set to true in the config file")
        #parser.add_argument('--string-

        parser.add_argument('--forced', action="store_true", default=False,
            help='Overwrite the ExpressionData.csv file if it already exists.')
    # these options are now in the config file
    #parser.add_argument('--string-target', action='store_true', default=False,
    #        help='Keep STRING edges only for the target prots')
    #parser.add_argument('--ssn-target', action='store_true', default=False,
    #        help='Keep SSN edges only between target and non-target prots')
    parser.add_argument('--stats-only', action='store_true', default=False,
            help="Print statistics about network sizes. Skip prediction methods")
    parser.add_argument('--skip-core', action='store_true', default=False,
            help="Skip running and evaluating the species in the core. " +
            "Used for the COMP and ELEC evaluations to run/eval only the species without EXPC annotations")

    return parser


def parse_args(parser):
    opts = parser.parse_args()
    kwargs = vars(opts)
    # store the option as 'compute_smin' since that seems more intuitive
    kwargs['compute_smin'] = not kwargs['skip_smin']

    config_file = opts.config
    with open(config_file, 'r') as conf:
        #config_map = yaml.load(conf, Loader=yaml.FullLoader)
        config_map = yaml.load(conf)

    return config_map, kwargs


if __name__ == "__main__":
    # first load the run_eval_algs parser
    parser = run_eval_algs.setup_opts()
    parser = setup_parser(parser)
    config_map, kwargs = parse_args(parser)

    main(config_map, **kwargs)
