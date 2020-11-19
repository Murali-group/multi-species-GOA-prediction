# The zip file has the following directores: annotations, networks, sequences, taxons
# This script will put these directories inside the inputs folder
# and setup the network files the way the scripts are expecting.

#mkdir inputs; cd inputs
# TODO download the data from zenodo
#wget
#unzip

# make the input directory file structure the config files are expecting
mkdir net-versions; cd net-versions
# rather than copy the networks, make symbolic links
mkdir 2019_10-sp200-eval-e0_1; cd 2019_10-sp200-eval-e0_1
ln -s ../../networks/sp200-eval-e0_1-net.txt.gz
cd ..

mkdir 2019_10-sp200-eval-e0_1-stringv11-700; cd 2019_10-sp200-eval-e0_1-stringv11-700
ln -s ../../networks/sp200-eval-e0_1-net.txt.gz
ln -s ../../networks/sp200-stringv11-700-net.txt.gz

cd ../../
gunzip sequences/sp200-uniprot-mapping.tab.gz
