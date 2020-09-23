#!/usr/bin/python

print("Importing libraries")

from optparse import OptionParser
import os
import sys
import gzip
from tqdm import tqdm
import utils.file_utils as utils
#sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
#import fungcat_settings as f_settings
#import matplotlib
#matplotlib.use("Agg")
#import matplotlib.pyplot as plt
import pandas as pd
from collections import defaultdict
## From SVN: src/python/CellCycle 
#from annotations import Annotations
#from annotations import GOdag


def main(args):
    global species_to_name

    ## Parse command line args.
    usage = '%s [options]\n' % (sys.argv[0])
    parser = OptionParser(usage=usage)
    parser.add_option('-G', '--gaf-file', 
            help="Directory containing the annotations for each species (in GAIN format)")
    parser.add_option('-o', '--out-file', type='string', metavar='STR', default="outputs/viz/stats/goa-stats.txt",
            help="Directory to write the annotations stats to. Default: outputs/viz/stats/goa-stats.txt")
    #parser.add_option('-o'
    # default is non-iea
    #parser.add_option('', '--ev-codes',
    #                  help="Comma-separated evidence codes to consider when parsing the GAF file")
    # the taxon ID is in the GAF file
    parser.add_option('-P', '--uniprot-proteomes',
            help="File containing the proteome of each uniprot ID.")
    parser.add_option('', '--uniprot-taxon',
            help="File containing the taxon ID of each uniprot ID.")
    parser.add_option('-T', '--taxon-file', type='string', 
            help="taxon IDs and the species name for which to get stats")

    (opts, args) = parser.parse_args()

    #ev_codes = opts.ev_codes.split(',')
    #print("Using the following evidence codes: %s" % (','.join(ev_codes)))

    print("Reading %s" % (opts.uniprot_proteomes))
    uniprot_to_proteomes = utils.readDict(opts.uniprot_proteomes, 1,2)
    # a protein can have multiple ',' spearated proteomes
    # get rid of the extra stuff after the ':' here
    #uniprot_to_proteomes = {u: ','.join(set(p.split(':')[0] for p in proteomes.split(', '))) \
    #        for u, proteomes in uniprot_to_proteomes.items()}
    uniprot_to_taxon = utils.readDict(opts.uniprot_taxon, 1,2)
    taxon_proteomes = defaultdict(set)
    for u, p in uniprot_to_proteomes.items():
        t = uniprot_to_taxon[u]
        taxon_proteomes[t].add(p)
    print("%d uniprot_ids, %d proteomes" % (len(uniprot_to_proteomes), len(set(uniprot_to_proteomes.values()))))
    taxon_to_name = utils.readDict(opts.taxon_file, 1, 2)
    #print '\n'.join(['%s\t%s' % (s, name) for s, name in species_to_name])
    print("getting the stats for %d taxon ids" % (len(taxon_to_name)))
    df = goa_stats(opts.gaf_file, taxon_to_name, uniprot_to_proteomes)

    s = pd.Series({p: t for t, proteomes in taxon_proteomes.items() for p in proteomes})
    # set all of the proteomes which have no data to 0
    #df = df.reindex(s.index)
    df['#taxon'] = s
    print(df.head())
    df.reset_index(inplace=True, col_fill="proteome")
    df.set_index('#taxon', inplace=True)
    df['species'] = pd.Series(taxon_to_name)
    df.rename(columns={'index': "proteome"}, inplace=True)
    out_file = opts.out_file
    os.makedirs(os.path.dirname(out_file), exist_ok=True)  # only works with python3
    print("writing %s" % (out_file))
    df.to_csv(out_file, sep='\t', index_label="#taxon")


def goa_stats(gaf_file, selected_strains, uniprot_to_proteomes):
    """ For each strain/taxonomy ID in the given selected_strains list, 
        print some overview statistics about the number of annotations.
        Also write a table for the # of proteins annotated to each function
    """
    prots_skipped = set()
    lines_skipped = 0
    # get the MF stats of the selected_strains
    num_prots_per_term = defaultdict(int)
    num_noniea_prots_per_term = defaultdict(int)
    num_ann_per_prot = defaultdict(int)
    num_f_ann_per_proteome = defaultdict(int)
    num_p_ann_per_proteome = defaultdict(int)
    exp_num_f_ann_per_proteome = defaultdict(int)
    exp_num_p_ann_per_proteome = defaultdict(int)
    hierarchies = {} 
    gonames = {} 
    #df_per_strain = {}
    # keep track of which species each proteome belongs to
    taxon_proteomes = defaultdict(set)
    exp_evidence_codes = set(['EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'IEP'])

    print("Reading %s" % (gaf_file))
    with gzip.open(gaf_file, 'r') as f:
        # calculated previously
        total_lines = 697908163
        for line in tqdm(f, total=total_lines):
            line = line.decode()
            if line[0] == "!" or line[0] == "#" or line.rstrip() == '':
                continue

            # columns of this file are explained in the README: ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/README and http://geneontology.org/page/go-annotation-file-gaf-format-21
            # 
            #	Column  Contents
            #	1       DB
            #	2       DB_Object_ID
            #	3       DB_Object_Symbol
            #	4       Qualifier
            #	5       GO_ID
            #	6       DB:Reference
            #	7       Evidence Code
            #	8       With (or) From
            #	9       Aspect
            #	10      DB_Object_Name
            #	11      DB_Object_Synonym
            #	12      DB_Object_Type
            #	13      Taxon and Interacting taxon
            #	14      Date
            #	15      Assigned_By
            #	16      Annotation_Extension
            #	17      Gene_Product_Form_ID
            cols = line.rstrip().split('\t')
            aspect = cols[8]
            if aspect == 'C':
                continue
            # if multiple taxon are present, the first one is of this entry 
            curr_taxon = cols[12].split('|')[0].split(':')[-1]

            if selected_strains is not None and curr_taxon not in selected_strains:
                #print("Skipping taxon %s" % (curr_taxon))
                #print line
                #sys.exit()
                continue

            prot = cols[1]
            proteome = uniprot_to_proteomes.get(prot)
            goid = cols[4]
            hierarchies[goid] = aspect
            ev_code = cols[6]
            num_prots_per_term[goid] += 1
            if ev_code == "IEA":
                continue
            if proteome is None:
                prots_skipped.add(prot)
                lines_skipped += 1 
                continue
            num_noniea_prots_per_term[goid] += 1
            num_ann_per_prot[prot] += 1
            if aspect == 'P':
                num_p_ann_per_proteome[proteome] += 1
            else:
                num_f_ann_per_proteome[proteome] += 1
            if ev_code in exp_evidence_codes:
                if aspect == 'P':
                    exp_num_p_ann_per_proteome[proteome] += 1
                else:
                    exp_num_f_ann_per_proteome[proteome] += 1
            taxon_proteomes[curr_taxon].add(proteome) 
            if curr_taxon == "316": 
                print("316 Here %s" % (ev_code))
            #print(prot, curr_taxon, proteome, goid, aspect, ev_code)
            #print(exp_num_p_ann_per_proteome, exp_num_f_ann_per_proteome)
            #print(num_p_ann_per_proteome, num_f_ann_per_proteome)
            #sys.exit()

    print("%d non-iea proteins (%d total lines) skipped that weren't in a proteome" % (len(prots_skipped), lines_skipped))
    #out_str = "#taxon\tspecies\tproteome(s)\tnum_p_ann\tnum_f_ann\texp_num_p_ann\texp_num_f_ann\n"
    df = pd.DataFrame({
        'num_p_ann': num_p_ann_per_proteome,
        'num_f_ann': num_f_ann_per_proteome,
        'exp_num_p_ann': exp_num_p_ann_per_proteome,
        'exp_num_f_ann': exp_num_f_ann_per_proteome,
        }, dtype=int)
    df.fillna(0, inplace=True)
    return df
    #for taxon in selected_strains:
    #    out_str += "%s\t%s\t%s\t%d\t%d\t%d\t%d"



#    df = pd.concat(df_per_strain)
#    # move the species index back as a column
#    df = df.reset_index()
#    # the species index is set to level_0 by default
#
#    # TODO figure out how to get all of these counts into one big table. Should be doable with pandas
#    # count the number of annotations per species
#    print("Number of annotations per species:")
#    #print("\tMF:")
#    # not sure how to do this with pandas...
#    #print df[df['hierarchy'] == 'f'][['level_0', 'goid', 'orf']].groupby('level_0').nunique()
#    #print("\tBP:")
#    #print df[df['hierarchy'] == 'p'][['level_0', 'goid', 'orf']].groupby('level_0').nunique()
#    # count the number of proteins with an annotation per species 
#    print('\t'.join(['Taxon', 'Name', '# annotations', '# f annotations', '# p annotations', 'MF > 0 proteins annotated', 'BP > 0 proteins annotated',]))
#    for species in sorted(num_ann_per_species):
#        print('\t'.join([species, species_to_name[species], str(num_ann_per_species[species]),
#                         str(num_f_ann_per_species[species]), str(num_p_ann_per_species[species])]))
#
#    print("Number of proteins with an annotation per species:")
#    print("\tMF:")
#    print df[df['hierarchy'] == 'f'][['level_0', 'orf']].groupby('level_0').nunique()
#    print("\tBP:")
#    print df[df['hierarchy'] == 'p'][['level_0', 'orf']].groupby('level_0').nunique()
#    print("\tMF non-IEA:")
#    print df[(df['hierarchy'] == 'f') & (df['evidencecode'] != 'IEA')][['level_0', 'orf']].groupby('level_0').nunique()
#    print("\tBP non-IEA:")
#    print df[(df['hierarchy'] == 'p') & (df['evidencecode'] != 'IEA')][['level_0', 'orf']].groupby('level_0').nunique()
#
#    # repeat for only the GO experimental evidence codes (from here: http://geneontology.org/page/experimental-evidence-codes)
#    exp_evidence_codes = ['EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'IEP']
#    print("Number of proteins with an annotation from an experimental evidence code per species:")
#    print("\tMF:")
#    print df[(df['hierarchy'] == 'f') & (df['evidencecode'].isin(exp_evidence_codes))][['level_0', 'orf']].groupby('level_0').nunique()
#    print("\tBP:")
#    print df[(df['hierarchy'] == 'p') & (df['evidencecode'].isin(exp_evidence_codes))][['level_0', 'orf']].groupby('level_0').nunique()
#
#    df = pd.DataFrame({"num_prots": num_prots_per_term, "num_noniea_prots": num_noniea_prots_per_term, "hierarchy": hierarchies, "goname": gonames})
#
#    print("Stats for proteins:")
#    prot_df = pd.Series(num_ann_per_prot)
#    for cutoff in [0, 50, 100, 200]:
#        #print("\t%d terms have > %d proteins annotated" % (len(df[df['num_prots'] > cutoff]), cutoff))
#        print("\t%d proteins have > %d annotations" % (len(prot_df[prot_df > cutoff]), cutoff))
#
#    print("Stats for all GO terms:")
#    for cutoff in [0, 100, 500, 1000]:
#        print("\t%d terms have > %d proteins annotated" % (len(df[df['num_prots'] > cutoff]), cutoff))
#
#    print("Limiting to BP GO terms:")
#    bp_df = df[df["hierarchy"] == "p"]
#    for cutoff in [0, 100, 500, 1000]:
#        print("\t%d terms have > %d proteins annotated" % (len(bp_df[bp_df['num_prots'] > cutoff]), cutoff))
#
#    print("Limiting to MF GO terms:")
#    mf_df = df[df["hierarchy"] == "f"]
#    for cutoff in [0, 100, 500, 1000]:
#        print("\t%d terms have > %d proteins annotated" % (len(mf_df[mf_df['num_prots'] > cutoff]), cutoff))
#
#    print("Limiting to CC GO terms:")
#    cc_df = df[df["hierarchy"] == "c"]
#    for cutoff in [0, 100, 500, 1000]:
#        print("\t%d terms have > %d proteins annotated" % (len(cc_df[cc_df['num_prots'] > cutoff]), cutoff))
#
#    # now write these to an output file
#    utils.checkDir(out_dir)
#    out_file = "%s/ann-counts-%d-strains.tsv" % (out_dir, len(selected_strains))
#    print("Writing all annotation counts to: %s" % (out_file))
#    df.to_csv(out_file, sep='\t')
#    out_file = "%s/ann-counts-f-%d-strains.tsv" % (out_dir, len(selected_strains))
#    print("Writing MF annotation counts to: %s" % (out_file))
#    mf_df.to_csv(out_file, sep='\t')
#    out_file = "%s/ann-counts-p-%d-strains.tsv" % (out_dir, len(selected_strains))
#    print("Writing BP annotation counts to: %s" % (out_file))
#    bp_df.to_csv(out_file, sep='\t')
#
##print("terms with more than 1000 proteins annotated: ")
##print s[s > 1000]
#
##    # also plot the numbers
##    out_file = "test.png"
##    print("plotting %s" % (out_file))
##    plt.plot()
##    s[s < 2000].plot.hist(bins=20)
##    plt.savefig(out_file)
##    plt.close()


if __name__ == '__main__':
    main(sys.argv)
