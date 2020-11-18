
# This script sets up the core/target species and networks and performs the coresponding evaluations

import argparse
import os, sys
import numpy as np
from scipy import sparse as sp
from tqdm import tqdm
import yaml
import time
from collections import defaultdict
#import pandas as pd

# for some reason, the current directory isn't always in the path. this fixes that
if "" not in sys.path:
    sys.path.insert(0,"")
from src.annotation_prediction.src import main as run_eval_algs
from src.annotation_prediction.src.evaluate import eval_leave_one_species_out as eval_loso
from src.annotation_prediction.src.evaluate import eval_utils as eval_utils
from src.annotation_prediction.src import setup_sparse_networks as setup
from src.annotation_prediction.src.algorithms import runner as runner
#import src.annotation_prediction.src.go_term_prediction_examples.go_term_prediction_examples as go_examples
#import src.annotation_prediction.src.algorithms.aptrank_birgrank.run_birgrank as run_birgrank
import src.annotation_prediction.src.utils.config_utils as config_utils
import src.annotation_prediction.src.utils.file_utils as utils
import src.annotation_prediction.src.utils.ontology_utils as go_utils
import src.annotation_prediction.src.algorithms.alg_utils as alg_utils


def main(config_map, **kwargs):
    input_settings, input_dir, output_dir, alg_settings, kwargs \
        = config_utils.setup_config_variables(config_map, **kwargs)

    for dataset in input_settings['datasets']:
        kwargs = apply_dataset_settings_to_kwargs(dataset, **kwargs)
        print({key: val for key,val in kwargs.items() \
                if key not in ['species_to_uniprot_idx', 'alg_taxon_terms_to_skip']})
        #print(kwargs)
        net_obj, ann_obj, eval_ann_obj = run_eval_algs.setup_dataset(dataset, input_dir, **kwargs) 
        # if there are no annotations, then skip this dataset
        if len(ann_obj.terms) == 0:
            print("No terms found. Skipping this dataset")
            continue

        if kwargs['ssn_only']:
            print("Keeping only the SSN (first network)")
            net_obj = setup.Sparse_Networks(net_obj.sparse_networks[0], net_obj.nodes)

        # add the taxon file paths for this dataset to kwargs
        for arg in ['taxon_prot_file', 'only_taxon_file']:
            kwargs[arg] = "%s/%s" % (input_dir, dataset[arg]) if arg in dataset else None
        species_to_uniprot_idx = eval_loso.get_uniprot_species(kwargs['taxon_prot_file'], ann_obj)
        # set this in kwargs to use it in all functions
        kwargs['species_to_uniprot_idx'] = species_to_uniprot_idx

        if eval_ann_obj is not None:
            # for ELEC, I'm getting a memory error. 
            # I should be able to limit the memory by limiting the eval_ann_matrix to only the taxons that are currently being evaluated
            _, eval_taxons = eval_loso.get_selected_species(
                    species_to_uniprot_idx, kwargs['only_taxon_file'], kwargs['taxons'])
            eval_taxon_prots = get_taxon_prots(len(net_obj.nodes), eval_taxons, species_to_uniprot_idx)
            eval_ann_obj = limit_to_taxons(eval_taxon_prots, ann_obj=eval_ann_obj, **kwargs)

        # if specified, remove negative examples that are neighbors of positive examples
        if kwargs.get('rem_neg_neighbors'):
            ann_obj = rem_neg_neighbors(net_obj, ann_obj)
            if eval_ann_obj is not None:
                eval_ann_obj = rem_neg_neighbors(net_obj, eval_ann_obj)
        # the outputs will follow this structure:
        # outputs/<net_version>/<exp_name>/<alg_name>/output_files
        out_dir = "%s/%s/%s/" % (output_dir, dataset['net_version'], dataset['exp_name'])
        alg_runners = run_eval_algs.setup_runners(alg_settings, net_obj, ann_obj, out_dir, **kwargs)

        # If computing the smin, compute the information content here
        #if kwargs.get('compute_smin'):
        # store all of the prediction scores for the pos/neg examples to compute the Smin
        for run_obj in alg_runners:
            run_obj.all_pos_neg_scores = sp.csr_matrix(ann_obj.ann_matrix.shape)
            # also store a version of the eval matrix with only the most specific terms per taxon
            run_obj.eval_mat = sp.csr_matrix(ann_obj.ann_matrix.shape)
        if eval_ann_obj is not None:
            term_ic_vec = eval_utils.compute_information_content(eval_ann_obj)
            # need to align the term_ic_vec to the ann_obj
            mapped_term_ic_vec = np.zeros(ann_obj.ann_matrix.shape[0])
            for i, t in enumerate(eval_ann_obj.terms):
                if t in ann_obj.term2idx:
                    mapped_term_ic_vec[ann_obj.term2idx[t]] = term_ic_vec[i]
            term_ic_vec = mapped_term_ic_vec
            # now make sure the shape matches the other term_ic_vec
            term_ic_vec = term_ic_vec.reshape(len(term_ic_vec), 1)
        else:
            term_ic_vec = eval_utils.compute_information_content(ann_obj)
        kwargs['term_ic_vec'] = term_ic_vec
        #print("ann_matrix shape: %s, term_ic_vec shape: %s" % (str(ann_obj.ann_matrix.shape), str(term_ic_vec.shape)))

        # Run LOSO validation for each species in the core
        loso_experiments(alg_runners, net_obj, ann_obj, eval_ann_obj, **kwargs)


def apply_dataset_settings_to_kwargs(dataset, **kwargs):
    # pull the specified options from the dataset in the config file to kwargs
    # *ssnC*: keep the SSN among nodes in the core
    # *stringC*: keep the STRING edges among nodes in the core species
    # *stringT*: keep the STRING edges among nodes in the target species
    # *ssnNbrs*: compute scores among the core species, and then transfer them to the target taxons
    # *ssn_only*: keep only the SSN (first network). Useful to keep the prots the same as SSN+STRING
    # *add_neighbor_edges*: integer. add the neighbors of non-target taxon nodes up to k steps away
    # *limit_to_taxons_file*: limit all of the networks to the subgraph of prots in the given species
    # oracle_weights: use the annotations of the target species when running SWSN
    # rem_neg_neighbors: if a negative example has a positive example as a neighbor in the SSN, relabel it as an unknown example
    # youngs_neg: for a term t, a gene g cannot be a negative for t if g shares an annotation with any gene annotated to t 
    # sp_leaf_terms_only: for a given species, limit the terms to only those that are the most specific, meaning remove the ancestors of all terms
    for arg in ['ssnC', 'ssn_only', 'add_neighbor_edges',
                'stringC', 'stringT',
                'limit_to_taxons_file', 
                'oracle_weights', 'rem_neg_neighbors', 'youngs_neg',
                'sp_leaf_terms_only', 'ssnNbrs', 'exp_name']:
        kwargs[arg] = dataset.get(arg)
    return kwargs


def loso_experiments(alg_runners, net_obj, ann_obj, eval_ann_obj, **kwargs):
    # species_names is a dictionary of taxonomy ID to name, from the only_taxon_file
    # taxons_to_run is the list of taxons to evaluate, which comes from the only_taxon_file if none are passed in via kwargs
    species_names, taxons_to_run = eval_loso.get_selected_species(
        kwargs['species_to_uniprot_idx'], kwargs['only_taxon_file'], kwargs['taxons'])
    kwargs['alg_taxon_terms_to_skip'] = eval_loso.get_already_run_terms(alg_runners, **kwargs) 
    kwargs['num_test_cutoff'] = kwargs.get('num_test_cutoff', 10) 
    # keep track of the original to use it later
    orig_net_obj = net_obj
    # limit the networks to the given taxons 
    # TODO this is now the only way to run this script
    #if kwargs.get('limit_to_taxons_file'):
    # read in the specified taxons from the file
    _, core_taxons = eval_loso.get_selected_species(
        kwargs['species_to_uniprot_idx'], kwargs['limit_to_taxons_file'])
    core_taxon_prots = get_taxon_prots(
        len(net_obj.nodes), core_taxons, kwargs['species_to_uniprot_idx'])
    core_net_obj, core_ann_obj = limit_to_taxons(core_taxon_prots, net_obj=net_obj, ann_obj=ann_obj, **kwargs)
    # can't print the # nodes and edges here since they can change for SWSN

    # run LOSO on only the core taxons, since we can run non-core species all at the same time
    core_taxons_to_run = set(taxons_to_run) & set(core_taxons)
    core_taxons_to_run = {t: species_names[t] for t in core_taxons_to_run}
    #core_taxons_to_run = {t: species_names[t] for t in taxons_to_run}
    print("EVAL CORE: %d core taxons, %d for which to run LOSO" % (len(core_taxons), len(core_taxons_to_run)))

    # skip running LOSO on the core
    if kwargs.get('skip_core'):
        print("\tSkipping the core taxons because --skip-core was specified")
        params_results = {}
    else:
        params_results = run_loso_eval(
            core_taxons_to_run, alg_runners, core_net_obj, 
            core_ann_obj, eval_ann_obj, core_taxon_prots, **kwargs)

    # Now run and evaluate for all the target taxons (not in the core)
    target_taxons = set(taxons_to_run) - set(core_taxons)
    print("EVAL TARGET: %d core taxons, %d target taxons for which to evaluate" % (len(core_taxons), len(target_taxons)))
    # if there are non-core species to evaluate, do that here
    if eval_ann_obj is not None and len(target_taxons) > 0 \
       and kwargs.get('limit_to_taxons_file'): 
        # limit the original net to the core and target taxons
        target_taxon_prots = get_taxon_prots(
            len(orig_net_obj.nodes), target_taxons, kwargs['species_to_uniprot_idx'])
        # limit the SSN to the prots in the core+target
        core_target_prots = core_taxon_prots + target_taxon_prots
        orig_net_obj = limit_to_taxons(core_target_prots, net_obj=orig_net_obj, **kwargs)

        curr_params_results = run_eval_target_sp(
            target_taxons, target_taxon_prots, core_taxon_prots,
            alg_runners, orig_net_obj, core_net_obj,
            core_ann_obj, eval_ann_obj, **kwargs)

        # combine the params_results from both
        for key, val in curr_params_results.items():
            if key not in params_results:
                params_results[key] = val
            else:
                params_results[key] += val

#    if kwargs.get('stats_only'):
#        # write the connected component stats to a file(?)
#        out_file = "outputs/stats/%s-cc-stats%s.tsv" % (
#            kwargs.get('exp_name'), kwargs.get('postfix'))
#        os.makedirs(os.path.dirname(out_file), exist_ok=True)
#        while os.path.isfile(out_file):
#            out_file = out_file.replace('.tsv','t.tsv')
#        print("writing %s" % (out_file))
#        df_cc_stats.to_csv(out_file, sep='\t')
#        return params_results

    if kwargs.get('compute_smin'):
        for run_obj in alg_runners:
            out_file = run_obj.out_pref + '-smin.txt'
            eval_utils.compute_smin(run_obj.all_pos_neg_scores, run_obj.eval_mat, kwargs['term_ic_vec'], out_file=out_file, verbose=True)

    # I'm not resetting the params_results in run_obj,
    # so those params_results were already appended to.
    # just built this dict to make sure its not empty
    alg_params_results = defaultdict(int)
    for run_obj in alg_runners:
        for key in run_obj.params_results:
            alg_params_results[key] += run_obj.params_results[key]
        # this is already set in eval_loso.run_and_eval_algs()
        #run_obj.out_pref = "%s/loso%s%s" % (
        #    run_obj.out_dir, run_obj.params_str, kwargs.get("postfix", ""))
    if len(params_results) != 0 or len(alg_params_results) != 0:
        # write the running times to a file
        eval_loso.write_stats_file(alg_runners, params_results)
        alg_params_results.update(params_results)
        print("Final running times: " + ' '.join([
            "%s: %0.4f" % (key, val) for key, val in sorted(alg_params_results.items()) if 'time' in key]))
    print("")
    return params_results


def run_eval_target_sp(
        target_taxons, target_taxon_prots, core_taxon_prots,
        alg_runners, orig_net_obj, core_net_obj,
        core_ann_obj, eval_ann_obj, **kwargs):
    """
    Run the algorithms using the annotations of the core species only,
    and evaluate on all of the target species. Should be disjoint from the core species.
    For example, compute scores with the EXPC core, and evaluate the species with COMP and no EXPC
    """
    # the annotations were alread limited to the core species
    core_ann_mat, terms = core_ann_obj.ann_matrix, core_ann_obj.terms
    eval_terms = eval_ann_obj.terms
    # get the test ann mat to evaluate and figure out which terms to run
    diag = sp.diags(target_taxon_prots)
    eval_ann_obj.ann_matrix = eval_ann_obj.ann_matrix.dot(diag)
    # need to re-align the test ann mat (e.g., COMP) so it matches the core ann matrix (e.g., EXPC)
    eval_ann_obj.reshape_to_terms(terms, core_ann_obj.dag_matrix)
    test_ann_mat = eval_ann_obj.ann_matrix
    # now figure out how many terms pass the # ann cutoffs
    # it will be split per taxon later
    num_train_pos_per_term = (core_ann_mat > 0).sum(axis=1)
    num_test_pos_per_term = (test_ann_mat > 0).sum(axis=1)
    num_test_cutoff = kwargs['num_test_cutoff']
    terms_to_run = []
    for i, term in enumerate(terms):
        if num_train_pos_per_term[i] < num_test_cutoff or \
           num_test_pos_per_term[i] < num_test_cutoff:
            continue
        terms_to_run.append(term)
    print("\t%d with >= %d annotations among the %d core and %d eval terms " % (
        len(terms_to_run), num_test_cutoff, len(terms), len(eval_terms)))

    if kwargs.get('sp_leaf_terms_only'):
        # limit the terms to only the most specific per species
        sp_leaf_terms = set()
        for t in target_taxons:
            taxon_eval_ann_obj = get_taxon_eval_ann_obj(t, eval_ann_obj, **kwargs)
            eval_terms = set() 
            # since the terms to evaluate could be different per runner (because of terms already in the output file),
            # loop over each of the algs
            for run_obj in alg_runners:
                eval_terms |= set(get_terms_to_eval(
                    t, run_obj.name, taxon_eval_ann_obj, **kwargs))
            leaf_terms = go_utils.get_most_specific_terms(eval_terms, ann_obj=taxon_eval_ann_obj)
            sp_leaf_terms.update(leaf_terms)
        terms_to_run = set(terms_to_run) & sp_leaf_terms 
        print("\t'sp_leaf_terms_only': %d terms are most specific across the %d target species" % (
            len(sp_leaf_terms), len(target_taxons)))

    # split the running time for the core and target species
    for run_obj in alg_runners:
        for key, val in run_obj.params_results.copy().items():
            if ('process_time' in key or 'wall_time' in key) and '_core_loso' not in key:
                run_obj.params_results[key+'_core_loso'] = val

    if kwargs.get('ssnNbrs'):
        # For all other target species, precompute the scores and then transfer (ssnLocal),
        params_results = compute_core_scores_then_transfer(
            target_taxons, target_taxon_prots, core_taxon_prots,
            alg_runners, orig_net_obj, core_net_obj, core_ann_obj, eval_ann_obj,
            terms_to_run, **kwargs)
    else:
        # connect all the target networks to the core, then compute the scores
        params_results = connect_target_species_run_eval(
            target_taxons, target_taxon_prots, core_taxon_prots,
            alg_runners, orig_net_obj, core_net_obj, core_ann_obj, eval_ann_obj,
            terms_to_run, **kwargs)

    # split the running time for the core and target species
    for run_obj in alg_runners:
        for key, val in run_obj.params_results.copy().items():
            if ('process_time' in key or 'wall_time' in key) \
                    and '_core_loso' not in key and '_target_only' not in key:
                try:
                    # to get the time for the target only, subtract the LOSO time of the core species
                    target_only = val - run_obj.params_results[key+'_core_loso']
                    run_obj.params_results[key+'_target_only'] = target_only
                # if for some reason key+'_core_loso' isn't in the params_results, don't let the code fail
                except KeyError:
                    pass
    return params_results


def connect_target_species_run_eval(
        target_taxons, target_taxon_prots, core_taxon_prots,
        alg_runners, orig_net_obj, core_net_obj, core_ann_obj, eval_ann_obj,
        terms_to_run, **kwargs):
    params_results = defaultdict(int)
    print("\nConnecting %d target taxons to the core species network" % (len(target_taxons)))

    # Build the SSN from the core to target species
    # need to remove the inter-target-species edges to get only core -> target edges
    # the first net is the SSN
    ssn = orig_net_obj.sparse_networks[0] if orig_net_obj.multi_net else orig_net_obj.W
    ssnC = core_net_obj.sparse_networks[0] if core_net_obj.multi_net else core_net_obj.W
    # Get the edges from the core species to all target species
    ssnT, _ = get_core_to_target_ssn(ssn, target_taxon_prots)
    # also get the species-specific SSN edges
    # but not edges between target species
    for t in target_taxons:
        curr_target_taxon_prots = get_taxon_prots(ssn.shape[0], [t], kwargs['species_to_uniprot_idx'])
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
    test_mat = eval_ann_obj.ann_matrix 
    # now apply the network combination filters to get the specified combinations of networks
    net_obj = apply_net_comb_filters(
        train_mat, target_taxon_prots,
        net_obj, core_ann_obj, **kwargs)
    num_nodes, num_edges = print_net_stats(
        net_obj.W, train_mat, test_mat,
        term_idx=[core_ann_obj.term2idx[t] for t in terms_to_run])
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

        # sum the boolean of the columns, then use nonzero to get the columns with a nonzero value
        nodes_to_eval = (test_mat != 0).sum(axis=0).nonzero()[1]
        # limit the scores stored to the prots which will be evaluated
        print("\trunning %s storing scores for only the %d pos/neg nodes" % (run_obj.name, len(nodes_to_eval)))
        run_obj.target_prots = nodes_to_eval
        run_obj.kwargs['nodes_to_run'] = nodes_to_eval

        # run each method
        run_obj.setupInputs()
        run_obj.run()
        run_obj.setupOutputs()
        if kwargs.get('verbose'):
            utils.print_memory_usage()

        # finally, evaluate each sp-taxon pair
        out_file = "%s/loso%s%s.txt" % (
            run_obj.out_dir, run_obj.params_str, kwargs.get("postfix", ""))
        run_obj.out_pref = out_file.replace('.txt','')
        utils.checkDir(os.path.dirname(out_file))

        # # 2020-09-01 TEMP: skip the all taxons. Taking too long for ELEC
        # now evaluate 
        # first evaluate all taxons together
        eval_utils.evaluate_ground_truth(
            run_obj, eval_ann_obj, out_file.replace('loso','all-taxons'),
            taxon='all', append=False, **kwargs)

        # evaluate each of the sp-term pairs
        for t in tqdm(sorted(target_taxons)):
            print("Taxon: %s" % (t))
            terms_to_eval = eval_taxon(t, run_obj, eval_ann_obj, terms_to_run, out_file, **kwargs) 

        if kwargs['compute_smin']:
            run_obj.all_pos_neg_scores += eval_utils.store_pos_neg_scores(run_obj.term_scores, test_mat)
            eval_utils.store_terms_eval_mat(run_obj, eval_ann_obj, test_mat, specific_terms=terms_to_run) 

    return params_results


def get_terms_to_eval(
        t, alg, eval_ann_obj, num_test_cutoff=10, alg_taxon_terms_to_skip=None, **kwargs):
    test_mat, terms = eval_ann_obj.ann_matrix, eval_ann_obj.terms
    # and only evaluate terms that have at least 10 ann
    num_test_pos_per_term = (test_mat > 0).sum(axis=1)
    terms_to_eval = [term for i, term in enumerate(terms) \
                        if num_test_pos_per_term[i] >= num_test_cutoff]
    print("\t%d/%d terms with >= %d annotations for taxon %s" % (
        len(terms_to_eval), len(terms), num_test_cutoff, t))
    if t in alg_taxon_terms_to_skip[alg]:  
        terms_to_skip = alg_taxon_terms_to_skip[alg][t]
        terms_to_eval = [term for term in terms_to_eval if term not in terms_to_skip]
        print("\t%d to eval that aren't in the output file yet." % (len(terms_to_eval)))
    return terms_to_eval


def compute_core_scores_then_transfer(
        target_taxons, target_taxon_prots, core_taxon_prots,
        alg_runners, orig_net_obj, core_net_obj, core_ann_obj, eval_ann_obj,
        terms_to_run, **kwargs):
    params_results = defaultdict(int)
    print("\nComputing scores for the core species, then transferring to %d target species using ssnLocal" % (len(target_taxons)))

    ssn = orig_net_obj.sparse_networks[0]

    # Setup/combine (e.g., SWSN) the core networks 
    train_mat = core_ann_obj.ann_matrix
    test_mat = eval_ann_obj.ann_matrix
    core_net_obj, _ = apply_net_comb_filters(
        train_mat, target_taxon_prots, core_net_obj, core_ann_obj, **kwargs)
    num_nodes, num_edges = print_net_stats(
        core_net_obj.W, train_mat, test_mat,
        term_idx=[core_ann_obj.term2idx[t] for t in terms_to_run])
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
        print("%d pos, %d neg" % (len((train_mat > 0).data), len((train_mat < 0).data)))

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
        core_to_target_ssn, _ = get_core_to_target_ssn(ssn, target_taxon_prots)
        # and run Local with that SSN
        run_obj = run_ssn_local(run_obj, core_to_target_ssn)

        if alg == "genemania":
            # BUG: for GM, there are unreachable target nodes. Those should get a score of the value of k
            # but they got a score of 0. Since most of the nodes in the target have a negative score,
            # this put unreachable left-out negative examples at a much higher score.
            # This function fixes that bug by setting the unknown examples to k in the target species
            run_obj.term_scores = fix_gm_default(run_obj, target_taxon_prots)

        # finally, evaluate each sp-taxon pair
        out_file = "%s/loso%s%s.txt" % (
            run_obj.out_dir, run_obj.params_str, kwargs.get("postfix", ""))
        run_obj.out_pref = out_file.replace('.txt','')
        utils.checkDir(os.path.dirname(out_file))

        # # 2020-09-01 TEMP: skip the all taxons. Taking too long for ELEC
        # now evaluate 
        # first evaluate all taxons together
        eval_utils.evaluate_ground_truth(
            run_obj, eval_ann_obj, out_file.replace('loso','all-taxons'),
            taxon='all', append=False, **kwargs)

        # evaluate each of the sp-term pairs
        for t in tqdm(sorted(target_taxons)):
            print("Taxon: %s" % (t))
            terms_to_eval = eval_taxon(t, run_obj, eval_ann_obj, terms_to_run, out_file, **kwargs) 

        if kwargs['compute_smin']:
            run_obj.all_pos_neg_scores += eval_utils.store_pos_neg_scores(run_obj.term_scores, test_mat)
            eval_utils.store_terms_eval_mat(run_obj, eval_ann_obj, test_mat, specific_terms=terms_to_run) 

    return params_results


def get_taxon_eval_ann_obj(t, eval_ann_obj, **kwargs):
    # get just this taxon from the eval ann mat 
    taxon_prots = get_taxon_prots(len(eval_ann_obj.prots), [t], kwargs['species_to_uniprot_idx']) 
    diag = sp.diags(taxon_prots)
    test_mat = eval_ann_obj.ann_matrix.dot(diag)
    taxon_eval_obj = setup.Sparse_Annotations(
        eval_ann_obj.dag_matrix, test_mat, eval_ann_obj.terms, eval_ann_obj.prots)
    return taxon_eval_obj


def eval_taxon(t, run_obj, eval_ann_obj, terms_to_run, out_file, **kwargs):
    taxon_eval_obj = get_taxon_eval_ann_obj(t, eval_ann_obj, **kwargs)
    # also limit the scores to the current taxon to make processing faster
    all_scores = run_obj.term_scores
    taxon_prots = get_taxon_prots(len(eval_ann_obj.prots), [t], kwargs['species_to_uniprot_idx']) 
    diag = sp.diags(taxon_prots)
    taxon_scores = all_scores.dot(diag)
    run_obj.term_scores = taxon_scores 

    # and only evaluate terms that have at least 10 ann in the target species
    terms_to_eval = get_terms_to_eval(
        t, run_obj.name, taxon_eval_obj, **kwargs)
    terms_to_eval = set(terms_to_run) & set(terms_to_eval)
    if len(terms_to_eval) == 0:
        return
    run_obj.terms_to_run = terms_to_eval

    eval_utils.evaluate_ground_truth(
        run_obj, taxon_eval_obj, out_file,
        taxon=t, append=True, **kwargs)

    run_obj.term_scores = all_scores
    return terms_to_eval


def run_ssn_local(run_obj, core_to_target_ssn):
    print("Running local to transfer the computed scores to the target species using the SSN edges")
    P = alg_utils.normalizeGraphEdgeWeights(core_to_target_ssn)
    run_obj.term_scores = run_obj.term_scores.tocsr() 
    print("\t%d total gene-term pairs with scores before" % len(run_obj.term_scores.data))

    start_wall_time = time.time()
    start_process_time = time.process_time()
    # multiply the entire term_scores matrix with P, the ssn, to get the new scores
    run_obj.term_scores = P.dot(run_obj.term_scores.T).T
    print("\t%d total gene-term pairs with scores after" % len(run_obj.term_scores.data))

    # and store the times
    wall_time = time.time() - start_wall_time
    process_time = time.process_time() - start_process_time
    alg_name = "%s%s" % (run_obj.name, run_obj.params_str)
    print("\t%0.4f process time before" % (run_obj.params_results["%s_process_time"%alg_name]))
    run_obj.params_results["%s_wall_time"%alg_name] += wall_time
    run_obj.params_results["%s_process_time"%alg_name] += process_time
    print("\t%0.4f process time after" % (run_obj.params_results["%s_process_time"%alg_name]))
    return run_obj


def run_loso_eval(
        taxons_to_run, alg_runners, net_obj, 
        ann_obj, eval_ann_obj,
        core_taxon_prots, **kwargs):

    if kwargs.get('stats_only'):
        # keep track of some additional stats
        df_cc_stats = pd.DataFrame()
    for run_obj in alg_runners:
        if run_obj.name == 'async_rw':
            kwargs['async_rw'] = True
            if len(alg_runners) > 1:
                print("ERROR: 'async_rw' has a special network setup and should not be run at the same time as other methods.\nQuitting")
                sys.exit()

    print("\nRunning LOSO on %d core species" % (len(taxons_to_run)))
    # now loop through each of the taxons, subset the networks according to the options, and run LOSO validation
    params_results = defaultdict(int)
    for t, species_name in tqdm(sorted(taxons_to_run.items())):
        tqdm.write("\n" + "-"*30)
        tqdm.write("Taxon: %s - %s" % (t, species_name))

        if kwargs.get('verbose'):
            utils.print_memory_usage()
        # split the annotation matrix into testing (the annotations of this taxon) and training (all other annotations)
        # sp_terms are the terms for this species that pass the annotation cutoffs
        taxon_prot_idx = list(kwargs['species_to_uniprot_idx'][t])
        taxon_prots = get_taxon_prots(len(net_obj.nodes), [t], kwargs['species_to_uniprot_idx'])
        if kwargs.get('keep_ann'):
            # this function checks every term for overlap between the train and test annotations (super slow)
            # which is really only needed if there could be overlap (e.g., keeping all EXPC annotations and evaluating all COMP)
            train_ann_mat, test_ann_mat, sp_terms = eval_loso.leave_out_taxon(
                t, ann_obj, kwargs['species_to_uniprot_idx'],
                eval_ann_obj=eval_ann_obj, **kwargs)
        else:
            # split the annotation matrix (and eval ann mat) into train and test.
            # Don't need to worry about overlap since the train prots and test prots are disjoint
            train_ann_mat, test_ann_mat, sp_terms = eval_loso.split_ann_mat_train_test(
                taxon_prots, ann_obj, eval_ann_obj=eval_ann_obj, **kwargs)
        if kwargs.get('verbose'):
            utils.print_memory_usage()

        tqdm.write("\t%d/%d terms with >= %d annotations" % (len(sp_terms), len(ann_obj.terms), kwargs['num_test_cutoff']))
        if len(sp_terms) == 0:
            tqdm.write("\tskipping")
            continue
        # subset the terms if we want only the most specific species-term pairs
        if kwargs.get('sp_leaf_terms_only'):
            leaf_terms = go_utils.get_most_specific_terms(sp_terms, ann_obj=ann_obj)
            print("\t'sp_leaf_terms_only': %d/%d terms are most specific, or leaf terms" % (len(leaf_terms), len(sp_terms)))
            sp_terms = leaf_terms 

        # before weighting the networks, make sure at least one algorithms haven't already been run on these species term pairs
        if not kwargs.get('stats_only'): 
            terms_to_run = set()
            for run_obj in alg_runners:
                if t in kwargs['alg_taxon_terms_to_skip'][run_obj.name]:  
                    terms_to_run.update(set(sp_terms) - set(kwargs['alg_taxon_terms_to_skip'][run_obj.name][t]))
                else:
                    # this alg hasn't been run yet
                    terms_to_run.update(set(sp_terms))
            if len(terms_to_run) == 0:
                tqdm.write("\t0 terms to run. Either no algs were specified, or the terms are in the output file already. skipping")
                continue
            # if we're computing the smin, we need to make sure to run all terms
            elif kwargs.get('compute_smin'):
                terms_to_run = sp_terms

        # now subset the network(s) and weight them accordingly
        pos_train_ann = (train_ann_mat > 0).astype(int)
        # sum over the columns to get the # ann per gene
        ann_prots = np.ravel(pos_train_ann.sum(axis=0))
        tqdm.write("\t%d taxon prots, %d core prots, %d non-taxon prots with an annotation" % (
            taxon_prots.sum(), (core_taxon_prots - taxon_prots).sum(),
            np.count_nonzero(ann_prots)))
        new_net_obj = net_obj
        # UPDATE 2019-11-22: running for non-core taxons was moved outside of this loop
        ## add the extra SSN Target edges, if the network was already limitted, and this taxon isn't in the original set
        #if kwargs.get('limit_to_taxons_file'):
        #    #if t not in core_taxons:
        #        new_net_obj = limit_to_taxons(core_taxon_prots + taxon_prots, net_obj=orig_net_obj, **kwargs)
        #        #new_sparse_nets = add_ssn_target_taxon(orig_sparse_nets, sparse_nets, core_taxon_prots, taxon_prots, **kwargs)
        # oracle_weights is a special option to give the SWSN weighting the test annotations instead of the training annotations
        curr_ann_mat = train_ann_mat
        if kwargs.get('oracle_weights'):
            tqdm.write("Using the target taxon's annotations to weight the network")
            curr_ann_mat = test_ann_mat
        ann_obj.terms_to_run = terms_to_run if not kwargs.get('stats_only') else ann_obj.terms
        if kwargs.get('ssnNbrs'):
            new_net_obj, ssn_net_obj = apply_net_comb_filters(
                    curr_ann_mat, taxon_prots, new_net_obj, ann_obj, **kwargs)
        else:
            new_net_obj = apply_net_comb_filters(
                    curr_ann_mat, taxon_prots, new_net_obj, ann_obj, **kwargs)

        if not net_obj.weight_gmw:
            num_nodes, num_edges = print_net_stats(
                new_net_obj.W, train_ann_mat, test_ann_mat,
                term_idx=[ann_obj.term2idx[t] for t in sp_terms])
            params_results['%s_num_nodes|edges'%t] = "%d|%d" % (num_nodes, num_edges)
            if net_obj.weight_swsn:
                params_results['swsn_time'] += new_net_obj.swsn_time
                params_results['%s_swsn_weights'%t] = new_net_obj.swsn_weights
        if kwargs.get('stats_only'):
            df['taxon'] = t
            df_cc_stats = pd.concat([df_cc_stats, df])
            print("Skipping prediction methods")
            continue

        # now run each method
        for i, run_obj in enumerate(alg_runners):
            alg = run_obj.name
            tqdm.write("Alg: %s" % (alg))

            # and set the network to the SWSN weighted network in the runners
            run_obj.net_obj = new_net_obj
            # limit the run_obj to run on the terms for which there are annotations
            if t in kwargs['alg_taxon_terms_to_skip'][alg] \
                    and not run_obj.kwargs.get('debug_scores'):  
                terms_to_skip = kwargs['alg_taxon_terms_to_skip'][alg][t]
                sp_terms = [term for term in sp_terms if term not in terms_to_skip]
                tqdm.write("\t%d to run that aren't in the output file yet." % (len(sp_terms)))
                if len(sp_terms) == 0:
                    continue
            if kwargs.get('term') is not None:
                # if a subset of terms are specified, only run those
                print("\tspecified terms: %s" % (', '.join(kwargs['term'])))
                print("\t%d/%d specified terms overlap with the %d terms for this taxon" % (
                        len(set(sp_terms) & set(kwargs['term'])),
                        len(kwargs['term']), len(sp_terms)))
                sp_terms = list(set(sp_terms) & set(kwargs['term']))
                if len(sp_terms) == 0:
                    print("\tskipping")
                    continue
            # limit the run_obj to run on the terms for which there are annotations
            run_obj.terms_to_run = sp_terms
            # limit the scores stored to the current taxon's prots
            run_obj.target_prots = taxon_prot_idx

            new_run_obj = run_obj
            if kwargs.get('ssnNbrs'):
                # TODO make this less hacky
                # had to hack the runner object so that Local will be run after the alg is run on the core. Makes the rest of the pipeline work
                if run_obj.name in ['sinksource_bounds', 'sinksourceplus_bounds']:
                    print("skipping %s because it doesn't make sense in this context" % (run_obj.name))
                    continue
                # each of the run objects has setup_inputs called first, then run, then setup outputs. 
                # the setup outputs function is empty for FSS, GM and birgrank, so I can replace it with a call to transfer the scores using the SSN
                # after which the run_obj will be evaluated like any other
                new_run_obj = Runner_then_local(
                    run_obj.name, run_obj.net_obj, run_obj.ann_obj,
                    run_obj.out_dir, run_obj.params, **run_obj.kwargs)
                new_run_obj.ssn_net_obj = ssn_net_obj
                new_run_obj.terms_to_run = sp_terms
                new_run_obj.out_dir = run_obj.out_dir
                new_run_obj.params_results = run_obj.params_results
                # keep the scores of only the nodes in the core species when running so that they can be transferred
                curr_core_taxon_prots = core_taxon_prots - taxon_prots
                curr_core_taxon_prot_idx = np.where(curr_core_taxon_prots)[0]
                new_run_obj.target_prots = curr_core_taxon_prot_idx
                # also keep track of the target taxon prots
                new_run_obj.target_taxon_prots = taxon_prots
                # birgrank (and all 'gene-based' methods) uses the "nodes_to_run" option, so need to set that here as well
                if run_obj.get_alg_type() == 'gene-based':
                    # set it so that it runs on all nodes in the core taxon. 
                    new_run_obj.kwargs['nodes_to_run'] = curr_core_taxon_prot_idx
                    new_run_obj.kwargs['run_all_nodes'] = True
                if kwargs.get('compute_smin'):
                    new_run_obj.all_pos_neg_scores = run_obj.all_pos_neg_scores
                    new_run_obj.eval_mat = run_obj.eval_mat

            # now run the rest of the loso pipeline
            eval_loso.run_and_eval_algs(
                new_run_obj, ann_obj, 
                train_ann_mat, test_ann_mat,
                taxon=t, **kwargs)
            if kwargs.get('compute_smin') and kwargs.get('ssnNbrs'):
                run_obj.all_pos_neg_scores = new_run_obj.all_pos_neg_scores
                run_obj.eval_mat = new_run_obj.eval_mat
                run_obj.out_pref = new_run_obj.out_pref

            # to make sure the run_obj stats can be written, 
            # replace the params results of the orig run_obj with this new one
            alg_runners[i].params_results = new_run_obj.params_results
    return params_results


# redefine the setupOutputs function as a hack
class Runner_then_local(runner.Runner):

    #def __init__(self):
    #    self.origSetupOutputs = runner.Runner.setupOutputs

    def setupOutputs(self, **kwargs):
        # TODO move this back to the top of the function
        # call the original function
        #runner.Runner.setupOutputs(**kwargs) 
## this function is a hack to get the run_obj to transfer the scores it computed to a new set of nodes using local.
## The only issue is it will take more RAM since scores are stored for all nodes.
## Maybe I could 
#def hack_run_obj_to_run_local(run_obj, **kwargs):
        print("Running local to transfer the computed scores to the target species using the SSN edges")
        P = alg_utils.normalizeGraphEdgeWeights(self.ssn_net_obj.W)
        # make sure the target taxon doesn't already have scores
        target_taxon_prots = self.target_taxon_prots.astype(int)
        diag = sp.diags(target_taxon_prots)
        target_taxon_scores = self.term_scores.dot(diag)
        assert len(target_taxon_scores.data) == 0, \
            "Target taxon already has non-zero scores"

        self.term_scores = self.term_scores.tocsr()

        start_wall_time = time.time()
        start_process_time = time.process_time()
        print("\t%d total gene-term pairs with scores before" % len(self.term_scores.data))
        # multiply the entire term_scores matrix with P, the ssn, to get the new scores
        self.term_scores = P.dot(self.term_scores.T).T
        print("\t%d total gene-term pairs with scores after" % len(self.term_scores.data))

        if self.name == "genemania":
            # BUG: for GM, there are unreachable target nodes. Those should get a score of the value of k
            # but they got a score of 0. Since most of the nodes in the target have a negative score,
            # this put unreachable left-out negative examples at a much higher score.
            # This function fixes that bug by setting the unknown examples to k in the target species
            self.term_scores = fix_gm_default(self, target_taxon_prots)

        # and store the times
        wall_time = time.time() - start_wall_time
        process_time = time.process_time() - start_process_time
        alg_name = "%s%s" % (self.name, self.params_str)
        print("\t%0.4f process time before" % (self.params_results["%s_process_time"%alg_name]))
        self.params_results["%s_wall_time"%alg_name] += wall_time
        self.params_results["%s_process_time"%alg_name] += process_time
        print("\t%0.4f process time after" % (self.params_results["%s_process_time"%alg_name]))
        # call the original function for the debug_scores
        runner.Runner.setupOutputs(self, **kwargs) 


def fix_gm_default(run_obj, target_taxon_prots):
    # BUG: for GM, there are unreachable target nodes. Those should get a score of the value of k
    # but they got a score of 0. Since most of the nodes in the target have a negative score,
    # this put unreachable left-out negative examples at a much higher score.
    # I will fix this by setting the unknown examples to k in the target species
    term_num_pos = np.ravel((run_obj.ann_obj.ann_matrix > 0).astype(int).sum(axis=1))
    term_num_neg = np.ravel((run_obj.ann_obj.ann_matrix < 0).astype(int).sum(axis=1))
    #k = (num_pos - num_neg) / float(num_pos + num_neg)
    term_k = (term_num_pos - term_num_neg) / (term_num_pos + term_num_neg)
    # remove the terms for which we don't need scores
    idx2run = [run_obj.ann_obj.term2idx[t] for t in run_obj.terms_to_run]
    idx2run2 = np.ones(len(term_k), dtype=bool)
    idx2run2[idx2run] = False
    term_k[idx2run2] = 0
    # expand the term_k vector to have the same value for each prot 
    num_terms, num_prots = run_obj.term_scores.shape
    term_k_mat = np.ones((int(num_terms), int(target_taxon_prots.astype(int).sum())))
    term_k_mat = (term_k_mat.T*term_k).T
    # get the unreachable target nodes  (i.e., have a score = 0)
    unreachable_nodes = run_obj.term_scores[:,target_taxon_prots.astype(bool)] == 0
    print(unreachable_nodes.shape)
    print(term_k_mat.shape)
    # and get the default value of k for those nodes only
    target_default_scores = unreachable_nodes.astype(int).multiply(term_k_mat)
    # now set the default scores of the target prots to the term k value
    # transpose to get the prots on the rows
    gs_target_default = sp.lil_matrix((num_prots, num_terms))
    gs_target_default[target_taxon_prots.astype(bool)] = target_default_scores.T
    # then add these values to the term_scores matrix to set them as the default values
    run_obj.term_scores += gs_target_default.T
    print("\t%d total gene-term pairs with scores after setting unreachable nodes to GM default" % len(run_obj.term_scores.data))
    return run_obj.term_scores


def print_net_stats(W, train_ann_mat, test_ann_mat, term_idx=None):
    """
    Get statistics about the # of connected components, and the number of left-out positive examples that are in a component without a positive example
    *term_idx*: indices of terms for which to limit the train and test annotations 
        e.g., current terms for a species
    """
    # count the nodes with at least one edge
    num_nodes = np.count_nonzero(W.sum(axis=0))
    # since the network is symmetric, the number of undirected edges is the number of entries divided by 2
    num_edges = (len(W.data) / 2)
    print("\t%s nodes, %s edges" % (num_nodes, num_edges))
    # get the number of connected components (ccs), as well as the nodes in each cc
    num_ccs, cc_labels = sp.csgraph.connected_components(W, directed=False, return_labels=True)
    # now split the nodes into their respective connected components
    ccs = defaultdict(set)
    for i, label in enumerate(cc_labels):
        ccs[label].add(i)
    cc_sizes = [len(nodes) for cc, nodes in ccs.items()]
    print("%s connected_components; sizes: max: %d, 75%%: %d, median: %d, 25%%: %d" % (
        num_ccs, np.max(cc_sizes), np.percentile(cc_sizes, 75),
        np.median(cc_sizes), np.percentile(cc_sizes, 25)))
    print("\t%0.2f%% nodes in the largest cc" % ((np.max(cc_sizes) / float(num_nodes))*100))
    print("top 20 connected_components sizes:")
    print(' '.join(str(cc_size) for cc_size in sorted(cc_sizes, reverse=True)[:20]))

    # to get stats about the train and test pos and neg examples,
    # first limit the annotations to the terms with annotations for this species
    curr_train_mat, curr_test_mat = train_ann_mat, test_ann_mat
    if term_idx is not None:
        # get just the rows of the terms specified
        curr_train_mat = train_ann_mat[term_idx]
        curr_test_mat = train_ann_mat[term_idx]
    # sum over the columns to get the prots with at least 1 positive example
    train_pos_prots = np.ravel((curr_train_mat > 0).astype(int).sum(axis=0))
    test_pos_prots = np.ravel((curr_test_mat > 0).astype(int).sum(axis=0))
    check_frac_ccs(ccs, train_pos_prots, test_pos_prots)
    # also check for test neg prots
    test_neg_prots = np.ravel((curr_test_mat < 0).astype(int).sum(axis=0))
    check_frac_ccs(ccs, train_pos_prots, test_neg_prots, test_string='neg')
    train_neg_prots = np.ravel((curr_train_mat < 0).astype(int).sum(axis=0))
    check_frac_ccs(ccs, train_pos_prots, test_neg_prots, train_neg_prots=train_neg_prots, train_string='pos or neg', test_string='neg')

    # TODO move this somewhere else
#    # also check how many ccs have only positive or only negative examples
#    # and the proportion of positive to negative examples in ccs
#    train_neg_prots = np.ravel((curr_train_mat < 0).astype(int).sum(axis=0))
#    train_neg_prot_idx = set(list(train_neg_prots.nonzero()[0]))
#    cc_stats = {}
#    for cc, nodes in ccs.items():
#        num_pos_in_cc = len(nodes & train_pos_prot_idx)
#        num_neg_in_cc = len(nodes & train_neg_prot_idx)
#        if num_pos_in_cc > 1 or num_neg_in_cc > 1:
#            cc_stats[cc] = {
#                'num_pos': num_pos_in_cc,
#                'num_neg': num_neg_in_cc,
#                'pos/(pos+neg)': num_pos_in_cc / (num_pos_in_cc+num_neg_in_cc)}
#    df = pd.DataFrame(cc_stats).T
#    df[['num_pos', 'num_neg']] = df[['num_pos', 'num_neg']].astype(int)
#    df['pos/(pos+neg)'] = df['pos/(pos+neg)'].round(3)
#    df.sort_values('pos/(pos+neg)', inplace=True, ascending=False)
#    print("head and tail of ratio of pos to neg per cc:")
#    print(df.head())
#    print(df.tail())

    #print("%d ccs have no train ann, %d have test ann, %d have test ann and no train ann" % (
    #    len(ccs_no_train_ann), len(ccs_test_ann), len(ccs_no_train_ann & ccs_test_ann)))
    #return num_nodes, num_edges, df
    return num_nodes, num_edges


def check_frac_ccs(ccs, train_pos_prots, test_pos_prots, train_neg_prots=None, train_string='pos', test_string='pos'):
    # sum over the columns to get the prots with at least 1 positive example
    train_pos_prot_idx = set(list(train_pos_prots.nonzero()[0]))
    test_pos_prot_idx = set(list(test_pos_prots.nonzero()[0]))
    ccs_no_train_pos = set(cc for cc, nodes in ccs.items() \
                           if len(nodes & train_pos_prot_idx) == 0)
    ccs_test_pos     = set(cc for cc, nodes in ccs.items() \
                           if len(nodes & test_pos_prot_idx) != 0)
    if train_neg_prots is not None:
        train_neg_prot_idx = set(list(train_neg_prots.nonzero()[0]))
        ccs_no_train_neg = set(cc for cc, nodes in ccs.items() \
                            if len(nodes & train_neg_prot_idx) == 0)
        # limit the ccs considered to those that have no training positive or negative examples
        ccs_no_train_pos = ccs_no_train_pos & ccs_no_train_neg
    # for the target species, get the prots which have at least one annotation 
    # and find the connected_components with target species annotations, yet no train positive annotations
    # in other words, get the percentage of nodes with a test ann that are in a cc with no train ann
    nodes_ccs_test_pos = set()
    for cc in ccs_no_train_pos & ccs_test_pos:
        nodes_ccs_test_pos.update(ccs[cc] & test_pos_prot_idx)
    print("%d/%d (%0.2f%%) test %s prots are in a cc with no train %s" % (
        len(nodes_ccs_test_pos), len(test_pos_prot_idx),
        (len(nodes_ccs_test_pos) / float(len(test_pos_prot_idx)))*100,
        test_string, train_string))


def get_gmw_weights(train_ann_mat, net_obj, ann_obj):
    # each term gets its own weight
    term_weights = {}
    for term in tqdm(ann_obj.terms_to_run):
        idx = ann_obj.term2idx[term]
        y = train_ann_mat[idx,:].toarray()[0]
        _, _, weights = net_obj.weight_GMW(y, term)
        term_weights[term] = weights
    return term_weights


def limit_to_taxons(taxon_prots, net_obj=None, ann_obj=None, **kwargs):

    # get 1s at the indices of the taxon prots along the diagonal
    diag = sp.diags(taxon_prots)
    if net_obj is not None:
        sparse_nets = net_obj.sparse_networks if net_obj.multi_net else [net_obj.W]
        print("\tlimitting %d nets to %d taxon_prots" % (len(sparse_nets), np.count_nonzero(taxon_prots)))
        new_sparse_nets = []
        for sparse_net in sparse_nets:
            # get only the edges between the specified nodes
            new_sparse_net = diag.dot(sparse_net).dot(diag)
            new_sparse_nets.append(new_sparse_net)
        new_net_obj = setup.Sparse_Networks(new_sparse_nets, net_obj.nodes, 
                net_names=net_obj.net_names, weight_method=net_obj.weight_method, 
                verbose=kwargs.get('verbose',False))

    if ann_obj is not None:
        # also limit the annotations to these prots
        new_ann_mat = ann_obj.ann_matrix.dot(diag)
        new_ann_obj = setup.Sparse_Annotations(
                ann_obj.dag_matrix, new_ann_mat, ann_obj.terms, ann_obj.prots)
        print("\t%d pos annotations reduced to %d" % (
            len((ann_obj.ann_matrix > 0).astype(int).data),
            len((new_ann_mat > 0).astype(int).data)))

    if net_obj is not None and ann_obj is not None:
        return new_net_obj, new_ann_obj
    elif ann_obj is not None:
        return new_ann_obj
    else:
        return new_net_obj


def apply_net_comb_filters(
        train_ann_mat, taxon_prots, net_obj, ann_obj,
        stringC=True, stringT=True, 
        ssnC=True, 
        ssnNbrs=False, 
        **kwargs):
    """
    Apply the network combination filters
    *taxon_prots*: array with 1's at the indices of proteins/nodes in the target taxon, 0's at the rest
    *stringC*: keep STRING edges for core taxons
    *stringT*: keep STRING edges for target taxon
    *ssnC*: keep SSN edges between core taxon prots
    *ssnNbrs: 
    """
    sparse_nets = net_obj.sparse_networks if net_obj.multi_net is True else [net_obj.W]
    # the first net should be the SSN. it will be modified later
    ssn = sparse_nets[0]
    # only run SWSN and modify the string networks if they are passed in
    if net_obj.multi_net is True:
        # the rest are the string nets
        string_nets = sparse_nets[1:]
        string_net_names = net_obj.net_names[1:]
        if (ssnNbrs and not ssnC) or kwargs.get('async_rw'):
            # create two separate net objects. One for the core network edges, and one for the SSN edges.
            # if SSN should not be included in the core, then remove it
            print("\tkeeping only the STRING edges")
            net_obj = setup.Sparse_Networks(
                    string_nets, net_obj.nodes, 
                    net_names=string_net_names, 
                    weight_method=net_obj.weight_method)
        if net_obj.weight_swsn:
            # get the weights to use from SWSN
            W, process_time = net_obj.weight_SWSN(train_ann_mat)
            weights = net_obj.swsn_weights
        elif net_obj.weight_gmw:
            term_weights = get_gmw_weights(train_ann_mat, net_obj, ann_obj)
        # setup the STRING networks.
        # specifying neither of these options will give stringT+stringC
        if not stringC or not stringT:
            if not stringC:
                print("\tremoving STRING edges which are among the core")
                prots_to_keep = taxon_prots 
            if not stringT:
                print("\tremoving STRING edges of the target taxon")
                # flip the 0s and 1s
                prots_to_keep = 1 - taxon_prots
            # get a matrix with 1s at the diagonal positions of nodes to keep
            diag = sp.diags(prots_to_keep)
            new_string_nets = []
            for string_net in string_nets:
                # use matrix multiplication to get only the edges between nodes to keep
                new_string_net = diag.dot(string_net).dot(diag)
                new_string_nets.append(new_string_net) 
            string_nets = new_string_nets

    # now setup the SSN
    if ssnNbrs and net_obj.multi_net:
        core_to_target_ssn, target_ssn = get_core_to_target_ssn(ssn, taxon_prots)
        ssn_net_obj = setup.Sparse_Networks(core_to_target_ssn, net_obj.nodes)
        if not ssnC:
            # keep only the STRING edges for the core
            sparse_nets = string_nets
        else:
            # include the ssn core
            # remove the ssnT edges (SSN edges to/from target nodes) to get just the core ssn 
            core_ssn = ssn - target_ssn
            if kwargs.get('async_rw'):
                sparse_nets = string_nets
                ssn = core_ssn
            else:
                sparse_nets = [core_ssn] + string_nets
    elif not ssnC:
        new_ssn = remove_core_ssn(
            ssn, train_ann_mat, taxon_prots, **kwargs) 
        if kwargs.get('async_rw'):
            sparse_nets = string_nets if net_obj.multi_net else [new_ssn]
            ssn = new_ssn
        else:
            sparse_nets = [new_ssn] + string_nets if net_obj.multi_net else [new_ssn]
    elif net_obj.multi_net:
        # otherwise, just combine the original ssn with the string nets 
        if kwargs.get('async_rw'):
            sparse_nets = string_nets
        else:
            sparse_nets = [ssn] + string_nets

    # now integrate the new nets using the previously computed weights
    # and make a new network object with that weighted network
    if net_obj.multi_net is True:
        if net_obj.weight_swsn:
            temp_net_obj = setup.Sparse_Networks(sparse_nets, net_obj.nodes, weight_method='swsn')
            W = temp_net_obj.combine_using_weights(weights)
            new_net_obj = setup.Sparse_Networks(W, net_obj.nodes)
            # keep the weight string as SWSN for writing the output file
            new_net_obj.weight_str = temp_net_obj.weight_str
            new_net_obj.swsn_time = process_time
            new_net_obj.swsn_weights = weights
        elif net_obj.weight_gmw:
            new_net_obj = setup.Sparse_Networks(
                sparse_nets, net_obj.nodes, weight_method='gmw', term_weights=term_weights)
        if kwargs.get('async_rw'):
            new_net_obj.SSN = ssn
        if ssnNbrs:
            return new_net_obj, ssn_net_obj
    else:
        W = sparse_nets[0]
        new_net_obj = setup.Sparse_Networks(W, net_obj.nodes)
        if kwargs.get('async_rw'):
            new_net_obj.SSN = ssn
    return new_net_obj


def get_core_to_target_ssn(ssn, target_prots):
    """
    Get the SSN edges from the core (non-target nodes) to the target nodes

    *returns*: the core_to_target_ssn, 
        as well as the network of all edges connected to the target prots
    """
    # get a matrix with 1s at the diagonal positions of nodes of the target taxon
    diag = sp.diags(target_prots)
    # get edges between the target taxon only
    target_only_ssn = diag.dot(ssn).dot(diag)
    # this next multiplication will get the edges to/from the target nodes.
    # it also doubles the edges between the given nodes,
    # so subtract edges between those nodes to get them back to their original weight
    target_ssn = diag.dot(ssn) + ssn.dot(diag) - target_only_ssn
    # subtract the taxon only edges once again to get only the edges between the target taxon and the core taxons
    core_to_target_ssn = target_ssn - target_only_ssn
    return core_to_target_ssn, target_ssn


def remove_core_ssn(ssn, train_ann_mat, taxon_prots, 
        add_neighbor_edges=0, **kwargs):
    """
    *add_neighbor_edges*: Add edges of non-target neighbors that are the given # of steps away 
    """
    # now remove the SSN edges which are between non-target species 
    #tqdm.write("Removing SSN edges between non-target taxon prots")
    # replacing with print temporarily
    print("Removing SSN edges between core taxon prots")
    new_ssn = get_edges_of_nodes(ssn, taxon_prots)
    # *ssn_target_ann_only*: also remove nodes that have no annotations in train_ann_mat
    # if ssn_target_ann_only:
    #     # UPDATE: now also remove edges to core prots which don't have any annotations
    #     # a node can only be a negative example if it also has a positive example
    #     pos_train_ann = (train_ann_mat > 0).astype(int)
    #     # sum over the columns to get the # ann per gene
    #     ann_prots = np.ravel(pos_train_ann.sum(axis=0))
    #     # limit to 1 and add the taxon prots
    #     ann_prots = (ann_prots > 0).astype(int)
    #     ann_prots += taxon_prots.astype(int)
    #     print("\tkeeping %d prots that have an annotation, or are from the target taxon" % (np.sum(ann_prots)))
    #     # flip the 0s and 1s
    #     non_ann_prots = 1 - ann_prots
    #     print("\tremoving edges to/from non-target taxon prots with no annotations")
    #     # now set all of the non-annotated prot rows and columns to 0
    #     diag = sp.diags(non_ann_prots)
    #     # there's no edges between the non_ann_prots, so no need to handle those
    #     ssn_to_remove = diag.dot(new_ssn) + new_ssn.dot(diag)
    #     new_ssn = new_ssn - ssn_to_remove
#    # not going to use this
#    if add_neighbor_edges is not None and add_neighbor_edges > 0:
#        for i in range(add_neighbor_edges):
#            # add the edges to the neighbors of the non-target nodes
#            # nodes are the indices that have at least one incident edge
#            nodes = (np.ravel(new_ssn.sum(axis=0)) > 0).astype(int)
#            # subtract the taxon prots to get only the non-taxon-prots, 
#            # and make sure there are no -1s (taxon prots which have no edges)
#            non_taxon_prots = ((nodes - taxon_prots) > 0).astype(int)
#            print("\tadding neighbors of %d non-taxon prots (iteration %d)" % (
#                    np.count_nonzero(non_taxon_prots), i+1))
#            # now get the edges incident on those nodes in the original ssn 
#            non_taxon_prot_edges = get_edges_of_nodes(ssn, non_taxon_prots)
#            # make sure none of the edges are duplicated
#            # do that by subtracting the original edges, and getting all > 0 indices
#            edges_to_add = non_taxon_prot_edges - new_ssn
#            edges_to_add[edges_to_add < 0] = 0
#            edges_to_add.eliminate_zeros()
#            new_ssn = new_ssn + edges_to_add
#            print("\tssn %s nodes, %s edges" % (np.count_nonzero(new_ssn.sum(axis=0)), len(new_ssn.data) / 2))

    new_ssn = new_ssn.tocsr()
    new_ssn.eliminate_zeros()
    return new_ssn


#def add_ssn_target_taxon(orig_sparse_nets, sparse_nets, net_taxon_prots, taxon_prots, **kwargs):
#    #orig_ssn = orig_sparse_nets[0]
#    # limit the 
#    # get the SSN edges to the target taxon
#    #target_taxon_ssn = get_edges_of_nodes(orig_ssn, taxon_prots)
#    # limit these edges to the nodes of teh target taxon and the net_taxon prots
#    tqdm.write("\tadding %d nodes and %d edges for the target taxon SSN" % (
#        np.count_nonzero(target_taxon_ssn.sum(axis=0)), len(target_taxon_ssn.data) / 2))
#    # this function should only be called if the taxon is already in the ssn
#    # TODO if these edges are already in the net_obj, we don't want to add them again
#    new_sparse_nets = []
#    ssn = net_obj.sparse_networks[0]
#    new_ssn = ssn + target_taxon_ssn
#    new_sparse_nets[0] = ssn 
#
#    return new_net_obj 


def rem_neg_neighbors(net_obj, ann_obj, cutoff=0):
    print("Relabelling negative examples as unknown if they are neighbors of positive examples in the SSN")
    neg_mat = (ann_obj.ann_matrix < 0).astype(int)
    num_neg = len(neg_mat.data) 
    # relabel negative examples as unknown examples if they are neighbors of a positive node
    # get the ssn from the list of networks
    if net_obj.multi_net is True:
        ssn = net_obj.sparse_networks[0]
    else:
        ssn = net_obj.W
    # get an unweighted version of the network
    W = (ssn > cutoff).astype(int) 
    # get the neighbors of the positives first
    pos_mat = (ann_obj.ann_matrix > 0).astype(int)
    pos_neighbor_mat = W.dot(pos_mat.T).T
    # then limit them to the non-positive neighbors, and set them to 0
    # or I could add 1 to them, and then mod 1 to get everything back to -1, 0 and 1
    new_ann_mat = ann_obj.ann_matrix + pos_neighbor_mat
    # get the matrix of only negative examples
    new_neg_mat = (new_ann_mat < 0).astype(int)
    new_num_neg = len(new_neg_mat.data)
    # and add them back together
    new_ann_mat = pos_mat - new_neg_mat
    new_ann_obj = setup.Sparse_Annotations(
            ann_obj.dag_matrix, new_ann_mat, ann_obj.terms, ann_obj.prots)
    print("\t%d negative examples neighboring positive examples relabeled to unknown examples (%d negative examples before, %d after)." % (
        num_neg - new_num_neg, num_neg, new_num_neg))
    return new_ann_obj


def get_most_specific_ann(pos_mat, dag_matrix):
    # full_prop_ann_mat will have the annotations of the ancestor terms of all prot-term pairs
    # thus, we can remove the ancestor annotations to get only the most specific annotations 
    full_prop_ann_mat = sp.csr_matrix(pos_mat.shape).T
    #last_prop_ann_mat = sp.csr_matrix(pos_mat.shape)
    last_prop_ann_mat = pos_mat.T
    prop_ann_mat = pos_mat.copy().T
    # propagate all of the annotations up the DAG
    i = 0
    while True:
        i += 1
        #print("\tpropagating annotations (%s)..." % (i))
        prop_ann_mat = prop_ann_mat.dot(dag_matrix)
        diff = prop_ann_mat != last_prop_ann_mat
        if diff is True or diff.nnz != 0:
            # full_prop_ann_mat doesn't get the initial values of prop_ann_mat. Only the next level up
            full_prop_ann_mat += prop_ann_mat
            last_prop_ann_mat = prop_ann_mat
        else:
            break
    #print("\tdone!")
    # now change values > 1 in the full_prop_ann_mat to 1s 
    full_prop_ann_mat = (full_prop_ann_mat > 0).astype(int).T
    # and subtract them from the pos_mat to get the most specific ones
    spec_ann_mat = pos_mat - full_prop_ann_mat
    spec_ann_mat = (spec_ann_mat > 0).astype(int) 
    spec_ann_mat.eliminate_zeros() 
    return spec_ann_mat


def get_taxon_prots(num_nodes, taxons, species_to_uniprot_idx):
    taxon_prots = np.zeros(num_nodes)
    for t in taxons:
        taxon_prot_idx = list(species_to_uniprot_idx[t])
        taxon_prots[taxon_prot_idx] = 1
    return taxon_prots


def get_edges_of_nodes(W, nodes):
    # get a matrix with 1s at the diagonal positions of nodes of the target taxon
    diag = sp.diags(nodes)
    # this next step will keep the edges to/from the given nodes.
    # it also doubles the edges between the given nodes,
    # so subtract edges between those nodes to get them back to their original weight
    new_W = diag.dot(W) + W.dot(diag)
    node_node_edges = diag.dot(W).dot(diag)
    new_W = new_W - node_node_edges
    return new_W


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
