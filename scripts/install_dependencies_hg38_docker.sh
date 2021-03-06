#!/bin/bash

#Usage: bash install_dependencies.sh $path_to_setup_dir $path_to_DNASCAN_dir $path_to_ANNOVAR $path_to_gatk_download
#Example: bash install_dependencies.sh /home/local/ /home/DNAscan /home/annovar /home/gatk_download_dir

INSTALL_DIR=$1

DNASCAN_DIR=$2

apt-get install -y update

apt-get install -y vim

apt-get install -y python3

apt-get install -y perl

apt-get install -y wget bzip2

apt-get install -y ttf-dejavu

mkdir $INSTALL_DIR

mkdir $INSTALL_DIR/humandb

cd $INSTALL_DIR

wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh

chmod +x Miniconda3-latest-Linux-x86_64.sh

./Miniconda3-latest-Linux-x86_64.sh -b -p $INSTALL_DIR/Miniconda3/

export PATH=$INSTALL_DIR/Miniconda3/bin:$PATH

echo export PATH=$INSTALL_DIR/Miniconda3/bin:$PATH >> ~/.bashrc

conda config --add channels conda-forge

conda config --add channels defaults

conda config --add channels r

conda config --add channels bioconda

conda install -y python=3.8.0

conda install -y pysam=0.16.0.1

conda install -y samtools=1.10

conda install -y freebayes>=1.0.2

conda install -y bedtools>=2.25.0

conda install -y vcftools>=0.1.13

conda install -y bcftools=1.10.2

conda install -y hisat2>=2.1.0

conda install -y bwa=0.7.17

conda install -y rtg-tools=3.12

conda install -y multiqc>=1.2

conda install -y fastqc=0.11.9

conda install -y expansionhunter=4.0.2

conda install -y wham=1.8.0

conda install -y sambamba>=0.6.6

conda install -y samblaster>=0.1.24

cd $DNASCAN_DIR

mkdir hg38

cd hg38

wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz

gzip -d hg38.fa.gz

samtools faidx hg38.fa

nohup bwa index hg38.fa &

nohup hisat2-build hg38.fa hg38 &

apt-get update -qq

apt-get install -y -qq bzip2 gcc g++ make python zlib1g-dev

cd $INSTALL_DIR

mkdir manta

cd manta

wget https://github.com/Illumina/manta/releases/download/v1.6.0/manta-1.6.0.centos6_x86_64.tar.bz2

tar -xjf manta-1.6.0.centos6_x86_64.tar.bz2

export PATH=$INSTALL_DIR/manta-1.6.0.centos6_x86_64/bin:$PATH

echo export PATH=$INSTALL_DIR/manta-1.6.0.centos6_x86_64/bin:$PATH >> ~/.bashrc

cd $DNASCAN_DIR

sed "s|path_reference = \"\"|path_reference = \"$DNASCAN_DIR\/hg38\/hg38.fa\"|" scripts/paths_configs.py > scripts/paths_configs.py_temp

sed "s|path_hisat_index = \"\"|path_hisat_index = \"$DNASCAN_DIR\/hg38\/hg38\"|" scripts/paths_configs.py_temp > scripts/paths_configs.py

sed "s|path_bwa_index = \"\"|path_bwa_index = \"$DNASCAN_DIR\/hg38\/hg38.fa\"|" scripts/paths_configs.py > scripts/paths_configs.py_temp

mv scripts/paths_configs.py_temp scripts/paths_configs.py

chmod +x scripts/*

export PATH=$DNASCAN_DIR/scripts/:$PATH

echo export PATH=$DNASCAN_DIR/scripts/:$PATH >> ~/.bashrc

echo "###########################################IMPORTANT######################################################"
echo "Hisat2-build and bwa-index are still creating their indexes. Please wait untill they complete their task."
echo "You can check whether or not they are still running using the 'top' command"
echo "##########################################################################################################"
