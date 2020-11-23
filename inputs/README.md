## Download
https://zenodo.org/record/4280990/files/2020-fastsinksource-data.zip?download=1


## 200 Bacterial Species
The taxonomy IDs and names of the 200 bacterial species, as well as the number of proteins, and EXPC, COMP and ELEC annotations are available below. The sequences, as well as namespace mappings from UniProt ID to taxonomy ID, gene name, and STRING ID, are also available.

- 200 species BP counts: `annotations/sp200-BP-ann-counts.tsv`
- 200 species MF counts: `annotations/sp200-MF-ann-counts.tsv`
- Sequences: `sequences/2019-11-04-uniprot-sp200-sequences-shortid.fasta.gz`
- Namespace Mapping: `sequences/sp200-uniprot-mapping.tab.gz`

### EXPC Core, COMP and ELEC Target Groups
The BP and MF core and target species groups are in the taxons directory:

- 200 species taxonomy IDs and names: taxons/sp200-taxons.txt
- Core species (EXPC):
  - BP (25): `taxons/expc/expc-bp-taxons.txt`
  - MF (40): `taxons/expc/expc-mf-taxons.txt`

## Sequence Similarity Network (SSN)
The 200-species SSN is available as a tab-delimited edge list, where nodes are in the UniProt namespace and the weight of an edge is the −log of the BLAST E-value. Note that we used an E-value cutoff of 0.1, and that in the case of BLAST reciprocal comparisons (i.e., A vs B and B vs A), we kept the weaker score of the two.

- 200 species SSN: `networks/sp200-eval-e0_1-net.txt.gz`

## STRING Networks
In STRING v11, data was available for 173/200 species. They are all concatenated into a single, tab-delimited file with a STRING channel per column, and the nodes are mapped to the UniProt namespace. We used a cutoff of 700 on the "combined" STRING edge score.

- STRING networks: `networks/sp200-stringv11-700-net.txt.gz`

## EXPC Core Network
The SSN-C+STRING-C network (combined using SWSN) for the core EXPC species are available here:

- BP: `networks/expc-core/expc-bp-core-swsn-net.txt.gz`
- MF: `networks/expc-core/expc-mf-core-swsn-net.txt.gz`

## GO Annotations
Also available are the following annotations files.

- A GAF file containing the GO annotations of each gene for the 200 species.
  - `annotations/2019_10-sp200-goa.gaf.gz`
- The positive and negative examples for each BP and MF GO term, generated from the GAF file, are available in the folders listed below. In a "pos-neg" file, the rows are the genes (UniProt IDs), the columns are the terms, and the values are 1/0/-1 corresponding to positive, unknown, and negative examples, respectively.
  - `annotations/expc`
  - `annotations/comp`
  - `annotations/elec`

> Note that the negative examples in these files do not yet match what was used in the paper. 
> Specifically, these files use the definition that a gene _g_ is a negative example for a term _t_ 
>     "if _g_ was not directly annotated to _t_ or to an ancestor or descendant of _t_ in the GO DAG, and had at least one other annotation"
> 
> In order to match what was used in the paper, the "Youngs negatives" method is applied in the scripts, 
> which relabels _g_ as an unknown example if it shares an annotation with any of the genes annotated to _t_.

