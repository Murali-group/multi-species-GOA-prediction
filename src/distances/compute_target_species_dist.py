#! /usr/bin/python

# script to analyze the results of the leave-one-species-out evaluation
# specifically, compute the distance of left-out positives to positives,
#     left-out negatives to positives, left-out positives to negatives,
#     and left-out negatives to negatives


import argparse
import os, sys
import numpy as np
from scipy import sparse as sp
from tqdm import tqdm
import yaml
#from collections import defaultdict
import fcntl

# for some reason, the current directory isn't always in the path. this fixes that
if "" not in sys.path:
    sys.path.insert(0,"")
from scripts import experiments
from FastSinkSource import run_eval_algs
from FastSinkSource.src.evaluate import eval_leave_one_species_out as eval_loso
#from FastSinkSource.src import setup_sparse_networks as setup
#import FastSinkSource.src.go_term_prediction_examples.go_term_prediction_examples as go_examples
import FastSinkSource.src.algorithms.alg_utils as alg_utils
#import FastSinkSource.src.algorithms.aptrank_birgrank.run_birgrank as run_birgrank


def main(config_map, **kwargs):
    input_settings = config_map['input_settings']
    input_dir = input_settings['input_dir']
    alg_settings = config_map['algs']
    output_settings = config_map['output_settings']
    # combine the evaluation settings in the config file and the kwargs
    kwargs.update(config_map['eval_settings'])

    for dataset in input_settings['datasets']:
        # pull the ssn_target_only and string_target_only options from the config file
        kwargs['ssn_target_only'] = dataset.get('ssn_target_only')
        kwargs['ssn_target_ann_only'] = dataset.get('ssn_target_ann_only')
        kwargs['string_target_only'] = dataset.get('string_target_only')
        # oracle_weights: use the annotations of the target species when running SWSN
        kwargs['oracle_weights'] = dataset.get('oracle_weights')
        # rem_neg_neighbors: if a negative example has a positive example as a neighbor in the SSN, relabel it as an unknown example
        kwargs['rem_neg_neighbors'] = dataset.get('rem_neg_neighbors')
        # youngs_neg: for a term t, a gene g cannot be a negative for t if g shares an annotation with any gene annotated to t 
        kwargs['youngs_neg'] = dataset.get('youngs_neg')
        net_obj, ann_obj, eval_ann_obj = run_eval_algs.setup_dataset(dataset, input_dir, alg_settings, **kwargs) 
        # if there are no annotations, then skip this dataset
        if len(ann_obj.goids) == 0:
            print("No terms found. Skipping this dataset")
            continue

        # make sure the birgrank matrices are carried over
        #if 'birgrank' in alg_settings or 'aptrank' in alg_settings:
        #    dag_mat, pos_mat, dag_goids = ann_obj.dag_matrix, ann_obj.pos_matrix, ann_obj.dag_goids
        # if specified, remove negative examples that are neighbors of positive examples
        if kwargs.get('rem_neg_neighbors'):
            ann_obj = experiments.rem_neg_neighbors(net_obj, ann_obj)
            if eval_ann_obj is not None:
                eval_ann_obj = experiments.rem_neg_neighbors(net_obj, eval_ann_obj)
        if kwargs.get('youngs_neg'):
            obo_file = kwargs['youngs_neg']
            ann_obj = experiments.youngs_neg(ann_obj, obo_file, "%s/%s" % (input_dir,dataset['pos_neg_file']))
            if eval_ann_obj is not None:
                eval_ann_obj = experiments.youngs_neg(eval_ann_obj, obo_file, "%s/%s" % (input_dir,dataset['pos_neg_file']))
        #if 'birgrank' in alg_settings or 'aptrank' in alg_settings:
        #    ann_obj.dag_matrix, ann_obj.pos_matrix, ann_obj.dag_goids = dag_mat, pos_mat, dag_goids
        # the outputs will follow this structure:
        # outputs/<net_version>/<exp_name>/<alg_name>/output_files
        out_dir = "%s/%s/%s/" % (output_settings['output_dir'], dataset['net_version'], dataset['exp_name'])
        alg_runners = run_eval_algs.setup_runners(alg_settings, net_obj, ann_obj, out_dir, **kwargs)

        # add the taxon file paths for this dataset to kwargs
        for arg in ['taxon_file', 'only_taxon_file']:
            kwargs[arg] = "%s/%s" % (input_dir, dataset[arg]) 
        taxon_net_ann = get_loso_net_ann(alg_runners, net_obj, ann_obj, eval_ann_obj, **kwargs)

        exp_name = dataset['exp_name'].split('/')[-1]
        net_version = dataset['net_version'].split('/')[-1]
        out_dir = "outputs/viz/distances/%s/%s" % (exp_name, net_version)
        os.makedirs(out_dir, exist_ok=True)
        # now compute the shortest paths
        for new_net_obj, train_ann_mat, test_ann_mat, taxon in taxon_net_ann:
            compute_target_nontarget_distances(
                new_net_obj, ann_obj, train_ann_mat, test_ann_mat, out_dir, taxon) 


def compute_target_nontarget_distances(
        net_obj, ann_obj, train_ann_mat, test_ann_mat, out_dir, taxon, unweighted=False):
    P = alg_utils.normalizeGraphEdgeWeights(net_obj.W)
    nodes = net_obj.nodes

    pos_train_mat = (train_ann_mat > 0).astype(int)
    neg_train_mat = (train_ann_mat < 0).astype(int)
    pos_test_mat = (test_ann_mat > 0).astype(int)
    #neg_test_mat = (test_ann_mat < 0).astype(int)
    goid_pos_target_lengths = {}
    goid_neg_target_lengths = {}
    # for each term, compute the shortest path from the left-out positives to the pos of other species
    for goid in tqdm(ann_obj.goids):
        goid_idx = ann_obj.goid2idx[goid]
        pos_train_prots = pos_train_mat[goid_idx].nonzero()[1]
        neg_train_prots = neg_train_mat[goid_idx].nonzero()[1]
        pos_test_prots = pos_test_mat[goid_idx].nonzero()[1]
        #neg_test_prots = neg_test_mat[goid_idx].nonzero()[1]

        start_end_sets = [
            # find the shortest paths from the left-out positives to the positives of the other species
            (pos_test_prots, pos_train_prots, goid_pos_target_lengths),
            # from the left-out positives to the negatives
            (pos_test_prots, neg_train_prots, goid_neg_target_lengths),]
        start_idx = list(pos_test_prots)
        # returns the distance of each source (row) to every other node in the graph
        dist_matrix = sp.csgraph.shortest_path(P, directed=True, indices=start_idx)
        for _, end_idx, goid_any_target_lengths in start_end_sets:
            # only store the length of the shortest path to any of the targets
            any_target_lengths = {}
            for i in range(len(start_idx)):
                distances = dist_matrix[i][end_idx]
                if len(distances) > 0:
                    min_idx = np.argmin(distances)
                    # convert the distance back to the normalized weights?
                    #min_dist = 10**(-distances[min_idx])
                    min_dist = distances[min_idx]
                    min_dist = float("%0.6e"%min_dist)
                    end_node = nodes[end_idx[min_idx]]
                else:
                    min_dist = float("inf")
                    end_node = '-'
                any_target_lengths[(nodes[start_idx[i]],end_node)] = min_dist
            goid_any_target_lengths[goid] = any_target_lengths

    # distance from left-out positives to positives
    pos_pos_file = "%s/pos-pos%s.tsv" % (out_dir, "-unw" if unweighted else "")
    # distance from left-out positives to negatives
    pos_neg_file = "%s/pos-neg%s.tsv" % (out_dir, "-unw" if unweighted else "")
    for out_file, goid_any_target_lengths in [
            (pos_pos_file, goid_pos_target_lengths),
            (pos_neg_file, goid_neg_target_lengths)]:
        if len(ann_obj.goids) <= 2:
            out_file = out_file.replace(".tsv","-%s-%s.tsv" % (taxon, '-'.join(ann_obj.goids)))
        # now write to file
        write_path_lengths(goid_any_target_lengths, out_file, taxon=taxon, write_type='a')


def get_loso_net_ann(
        alg_runners, net_obj, ann_obj, eval_ann_obj, **kwargs):

    # first load the species per uniprot ID
    species_to_uniprot_idx = eval_loso.get_uniprot_species(kwargs['taxon_file'], ann_obj)
    selected_species, taxons = eval_loso.get_selected_species(species_to_uniprot_idx, kwargs['only_taxon_file'], kwargs['taxons'])

    for t in tqdm(sorted(taxons)):
        tqdm.write("\n" + "-"*30)
        tqdm.write("Taxon: %s - %s" % (
            t, selected_species[t]))

        # leave out the annotations for this taxon ID
        train_ann_mat, test_ann_mat, sp_goterms = eval_loso.leave_out_taxon(
            t, ann_obj, species_to_uniprot_idx,
            eval_ann_obj=eval_ann_obj, **kwargs)

        tqdm.write("\t%d/%d goterms with >= %d annotations" % (len(sp_goterms), len(ann_obj.goids), kwargs['num_test_cutoff']))
        if len(sp_goterms) == 0:
            tqdm.write("\tskipping")
            continue

        # remove edges not part of this target taxon
        taxon_prot_idx = np.asarray(list(species_to_uniprot_idx[t]))
        taxon_prots = np.zeros(len(net_obj.nodes))
        taxon_prots[taxon_prot_idx] = 1 
        pos_train_ann = (train_ann_mat > 0).astype(int)
        # sum over the columns to get the # ann per gene
        ann_prots = np.ravel(pos_train_ann.sum(axis=0))
        tqdm.write("\t%d taxon prots, %d non-taxon prots with an annotation" % (taxon_prots.sum(), np.count_nonzero(ann_prots)))
        # oracle_weights is a special option to give the SWSN weighting the test annotations instead of the training annotations
        curr_ann_mat = train_ann_mat
        if kwargs.get('oracle_weights'):
            tqdm.write("Using the target taxon's annotations to weight the network")
            curr_ann_mat = test_ann_mat
        new_net_obj = experiments.limit_net_to_target_taxon(
                curr_ann_mat, taxon_prots, net_obj, **kwargs)

        yield new_net_obj, train_ann_mat, test_ann_mat, t

        #print_net_stats(new_net_obj.W, taxon_prot_idx, train_ann_mat, test_ann_mat)
        #if kwargs.get('stats_only'):
        #    print("Skipping prediction methods")
        #    continue


def write_path_lengths(goid_path_lengths, out_file, taxon='-', goid='-', write_type='w'):
    """
    *write_type*: either 'w' for write, or 'a' for append
    """
    if write_type == 'w':
        print("\tWriting to %s" % (out_file))
    else:
        print("\tAppending to %s" % (out_file))
    os.makedirs(os.path.dirname(out_file), exist_ok=True)
    with open(out_file, write_type) as out:
        # lock it to avoid scripts trying to write at the same time
        fcntl.flock(out, fcntl.LOCK_EX)
        for goid, path_lengths in goid_path_lengths.items():
            out.write(''.join("%s\t%s\t%s\t%s\t%s\n" % (taxon, goid,
                p1, p2, str(length)) for (p1,p2), length in sorted(path_lengths.items())))
        fcntl.flock(out, fcntl.LOCK_UN)


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
    # add extra options here

    return parser


if __name__ == "__main__":
    # first load the run_eval_algs parser
    parser = run_eval_algs.setup_opts()
    parser = setup_parser(parser)
    opts = parser.parse_args()
    kwargs = vars(opts)
    config_file = opts.config
    with open(config_file, 'r') as conf:
        config_map = yaml.load(conf, Loader=yaml.FullLoader)

    main(config_map, **kwargs)
