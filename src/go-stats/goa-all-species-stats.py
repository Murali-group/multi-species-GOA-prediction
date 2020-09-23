#!/usr/bin/python

print("Importing libraries")

from optparse import OptionParser
from collections import defaultdict
import os
import sys
import gzip
from tqdm import tqdm
import utils.file_utils as utils


expc_ev_codes = {'EXP','IDA','IPI','IMP','IGI','IEP','TAS','IC'}
# 2019-01-31 UPDATE: High Throughput Evidence Codes
htpe_ev_codes = {'HTP', 'HDA', 'HMP', 'HGI', 'HEP'}
comp_ev_codes = {'ISS','ISO','ISA','ISM','IGC','IBA','IBD','IKR','IRD','RCA'}
elec_ev_codes = {'IEA'}
# don't include ND and NAS
all_ev_codes = expc_ev_codes | htpe_ev_codes | comp_ev_codes | elec_ev_codes

def main(args):
    ## Parse command line args.
    usage = '%s [options]\nCounts the number of BP and MF annotations in the GAF file.' % (sys.argv[0])
    parser = OptionParser(usage=usage)
    parser.add_option('-i', '--goa-annotations', 
                      help="File containing goa annotations downloaded from UniProt in GAF format. Should be compressed with gzip (i.e., ends with .gz)")
    parser.add_option('-o', '--out-file', type='string', metavar='STR',
                      help="File for which to write the # annotations per species.")
    parser.add_option('-T', '--only-taxons', type='string', 
                      help="File with taxonomy IDs in the first column for which to limit the stats. Either this, or --taxon-categories is required.")
    parser.add_option('', '--only-uniprots',
            help="Limit the uniprot IDs to those contained in this file. Can only be used in conjunction with the --only-taxons file.")
    # the categories.dmp file can be found here: ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxcat.tar.gz
    parser.add_option('-t', '--taxon-categories', type='string', #default='inputs/goa/2017_09/categories.dmp',
                      help="File containing the category (e.g., B for Bacteria, E for Eukaryote) for each taxon. Will limit to B. Either this or --only-taxons is required")
    parser.add_option('-H', '--hierarchy', action='append', type='string', #default='inputs/goa/2017_09/categories.dmp',
            help="GO Hierarchy to use for counting annotations. Can specify multiple (e.g., -H P -H F) Default: P")

    (opts, args) = parser.parse_args(args)

    if opts.goa_annotations is None or opts.out_file is None:
        sys.exit("--goa-annotations (-i) and --out-dir (-o) required")
    if opts.only_taxons is None and opts.taxon_categories is None:
        sys.exit("--only-taxons or --taxon-categories required")
    # check the hierarchy options
    if opts.hierarchy is None:
        opts.hierarchy = ['P']
    else:
        if len(set(opts.hierarchy) - set(['P','F','C'])) != 0:
            sys.stderr.write("ERROR: illegal value(s) for --hierarchy: %s\n" % (opts.hierarchy))
            sys.stderr.write("\tAllowed values: 'P', 'F', 'C'\n")
            sys.exit(1)
    print("Limitting annotations to the following GO hierarchies: %s" % (opts.hierarchy))
    opts.hierarchy = set(opts.hierarchy)

    # useful to find the top 200 species of bacteria with EXPC and COMP annotations
    if opts.taxon_categories is not None:
        print("Reading taxon categories from %s" % (opts.taxon_categories))
        taxon_categories = {}
        with open(opts.taxon_categories, 'r') as f:
            for line in f:
                line = line.rstrip().split('\t')
                taxon_categories[line[2]] = line[0]
                taxon_categories[line[1]] = line[0]

    # useful to limit the stats to a specific subset of species (e.g., core)
    if opts.only_taxons is not None:
        print("Reading taxons from %s" % (opts.only_taxons))
        only_taxons = utils.readDict(opts.only_taxons, 1, 2)
        print("\tlimitting to %d taxon IDs" % (len(only_taxons)))
    if opts.only_uniprots is not None:
        print("Limitting uniprot IDs for each taxon using %s" % (opts.only_uniprots))
        only_uniprots = utils.readItemSet(opts.only_uniprots, 1)
        uniprot_to_species = utils.readDict(opts.only_uniprots, 1,2)
        species_to_uniprot = defaultdict(set)
        for p in uniprot_to_species:
            species_to_uniprot[uniprot_to_species[p]].add(p)
        for t in species_to_uniprot:
            species_to_uniprot[t].discard(None) 

    prot_num_expc_ann = defaultdict(int)
    taxon_num_expc_ann = defaultdict(int)
    taxon_num_htpe_ann = defaultdict(int)
    taxon_num_comp_ann = defaultdict(int)
    taxon_num_elec_ann = defaultdict(int)
    taxon_num_all_ann = defaultdict(int)
    taxon_prot_expc_ann = defaultdict(set)
    taxon_prot_htpe_ann = defaultdict(set)
    taxon_prot_comp_ann = defaultdict(set)
    taxon_prot_elec_ann = defaultdict(set)
    taxon_prot_all_ann = defaultdict(set)
    uniprots_skipped = set() 
    taxons_skipped = set()
    num_not_ann = 0

    os.makedirs(os.path.dirname(opts.out_file), exist_ok=True)
    #out_file = "%s/all-bacteria-stats.txt" % (opts.out_dir)
    out_file = opts.out_file

    print("Reading %s" % (opts.goa_annotations))
    with gzip.open(opts.goa_annotations, 'rb') as f:
        total_lines = None
        # it takes a while to count the number of lines, so I stored the line count from a previous run here
        if 'all' in opts.goa_annotations:
            total_lines = 697908163
        elif 'gcrp' in opts.goa_annotations:
            total_lines = 170779190
        for line in tqdm(f, total=total_lines):
            line = line.decode('UTF-8')
            if line.rstrip() == '' or line[0] == '!':
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
            #aspect = line[8]
            #if aspect != 'P':
            #    continue
            # if multiple taxon are present, the first one is of this entry
            curr_taxon = cols[12].split('|')[0].split(':')[-1]
            # skip all non-bacteria taxon
            if opts.taxon_categories is not None \
                    and (curr_taxon not in taxon_categories \
                    or taxon_categories[curr_taxon] != "B"):
                taxons_skipped.add(curr_taxon)
                continue

            if opts.only_taxons and curr_taxon not in only_taxons:
                taxons_skipped.add(curr_taxon)
                continue

            hierarchy = cols[8]
            # limit to P and F
            #if hierarchy not in ['P', 'F']:
            if hierarchy not in opts.hierarchy:
                continue

            uniprotID = cols[1]
            if opts.only_uniprots and uniprotID not in only_uniprots:
                uniprots_skipped.add(uniprotID)
                continue

            if cols[3] == "NOT":
                num_not_ann += 1
                continue
            #goid = cols[4]
            evidencecode = cols[6]
            if evidencecode in all_ev_codes:
                taxon_num_all_ann[curr_taxon] += 1 
                taxon_prot_all_ann[curr_taxon].add(uniprotID)
            if evidencecode in expc_ev_codes:
                taxon_num_expc_ann[curr_taxon] += 1 
                taxon_prot_expc_ann[curr_taxon].add(uniprotID)
                prot_num_expc_ann[uniprotID] += 1 
            elif evidencecode in htpe_ev_codes:
                taxon_num_htpe_ann[curr_taxon] += 1 
                taxon_prot_htpe_ann[curr_taxon].add(uniprotID)
            elif evidencecode in comp_ev_codes:
                taxon_num_comp_ann[curr_taxon] += 1 
                taxon_prot_comp_ann[curr_taxon].add(uniprotID)
            elif evidencecode in elec_ev_codes:
                taxon_num_elec_ann[curr_taxon] += 1 
                taxon_prot_elec_ann[curr_taxon].add(uniprotID)

    print("Finished")
    print("")
    print("\t%d taxons skipped" % (len(taxons_skipped)))
    print("\t%d NOT annotations skipped" % (num_not_ann))
    if opts.only_uniprots:
        print("\t%d uniprot IDs not in the %d uniprot IDs of the %d specified species" % (len(uniprots_skipped), len(only_uniprots), len(only_taxons)))
    print("writing to file %s" % (out_file))
    with open(out_file, 'w') as out:
        out.write("taxon%s\t# all ann\t# ELEC ann\t# COMP ann\t# HTP ann\t# EXPC ann\t# total prots\t# all prots\t# ELEC prots\t# COMP prots\t# HTP prots\t# EXPC prots\n"%(
            "\tname" if opts.only_taxons else ""))
        for t in sorted(taxon_num_all_ann, key=taxon_num_all_ann.get, reverse=True):
            out.write("%s%s\t%s\n" % (
                t, 
                "\t%s. %s"%(only_taxons[t][0], only_taxons[t].split(' ')[1]) if opts.only_taxons else "",
                #"\t%s"%only_taxons[t] if opts.only_taxons else "",
                '\t'.join(str(x) for x in [
                    taxon_num_all_ann[t], taxon_num_elec_ann[t], 
                    taxon_num_comp_ann[t], taxon_num_htpe_ann[t], taxon_num_expc_ann[t],
                    len(species_to_uniprot[t]),
                    len(taxon_prot_all_ann[t]), len(taxon_prot_elec_ann[t]), 
                    len(taxon_prot_comp_ann[t]), len(taxon_prot_htpe_ann[t]), 
                    len(taxon_prot_expc_ann[t])])))

    # Make a histogram of the number of EXPC annotations per prot 
    num_anns = sorted(set(prot_num_expc_ann.values()))
    num_expc_ann_prot = {i: set(p for p, num in prot_num_expc_ann.items() if num == i) for i in num_anns}
    print("\n# prots\tfrac prots\t# annotations")
    num_prots = len(prot_num_expc_ann)
    for i, prots in num_expc_ann_prot.items():
        print("%d\t%0.3f\t%d" % (len(prots), len(prots)/float(num_prots), i))


if __name__ == '__main__':
    main(sys.argv)
