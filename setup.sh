#!/bin/bash

# The zip file has the following directores: annotations, networks, sequences, taxons
# This script will put these directories inside the inputs folder
# and setup the network files the way the scripts are expecting.

cd inputs
echo "Downloading the data from zenodo"
wget https://zenodo.org/record/4280990/files/2020-fastsinksource-data.zip?download=1
echo "unzipping"
unzip 2020-fastsinksource-data.zip
mv 2020-fastsinksource-data/* .
rmdir 2020-fastsinksource-data
rm 2020-fastsinksource-data.zip

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

