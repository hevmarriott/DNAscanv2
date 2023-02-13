#!/bin/bash

#Usage: bash install_dependencies.sh $path_to_setup_dir $path_to_DNASCAN_dir $path_to_ANNOVAR $path_to_gatk_download
#Example: bash install_dependencies.sh /home/local/ /home/DNAscan /home/annovar /home/gatk_download_dir

INSTALL_DIR=$1

DNASCAN_DIR=$2

apt-get install -y update

apt-get install -y vim

apt-get install -y python2

apt-get install -y perl

apt-get install -y wget bzip2

apt-get install -y ttf-dejavu

apt-get install -y wget gzip make git tcl tcllib tar gcc g++

mkdir $INSTALL_DIR

mkdir $INSTALL_DIR/humandb

cd $DNASCAN_DIR

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

conda install -y samtools=1.9

conda install -y bedtools=2.25.0

conda install -y sambamba=0.7.1

conda install -y samblaster=0.1.26

conda install -y vcftools=0.1.16

conda install -y bcftools=1.9

conda install -y hisat2=2.2.1

conda install -y bwa=0.7.17

conda install -y rtg-tools=3.12

conda install -y multiqc=1.10.1

conda install -y fastqc=0.11.9

conda install -y expansionhunter=3.2.2

conda install -y delly=0.8.3

conda install -y pysimplegui=4.60.4

conda install -y survivor=1.0.7

conda install -y perl=5.26.2=h470a237_0

conda install -y biopython=1.78

conda install -y perl-app-cpanminus

cpan YAML::XS

cpan Sort::Key::Natural

cd $DNASCAN_DIR

mkdir hg19

cd hg19

wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz

tar -zxvf chromFa.tar.gz

for i in chr1.fa chr2.fa chr3.fa chr4.fa chr5.fa chr6.fa chr7.fa chr8.fa chr9.fa chr10.fa chr11.fa chr12.fa chr13.fa chr14.fa chr15.fa chr16.fa chr17.fa chr18.fa chr19.fa chr20.fa chr21.fa chr22.fa chrY.fa chrX.fa chrM.fa; do cat $i >> hg19.fa ; rm $i ; done

rm chr*

samtools faidx hg19.fa

nohup bwa index hg19.fa &

nohup hisat2-build hg19.fa hg19 &

apt-get update -qq

apt-get install -y -qq bzip2 gcc g++ make python zlib1g-dev

cd $INSTALL_DIR

wget https://github.com/Illumina/strelka/releases/download/v2.9.10/strelka-2.9.10.centos6_x86_64.tar.bz2

tar -xjf strelka-2.9.10.centos6_x86_64.tar.bz2

export PATH=$INSTALL_DIR/strelka-2.9.10.centos6_x86_64/bin:$PATH

echo export PATH=$INSTALL_DIR/strelka-2.9.10.centos6_x86_64/bin:$PATH >> ~/.bashrc

wget https://github.com/Illumina/manta/releases/download/v1.6.0/manta-1.6.0.centos6_x86_64.tar.bz2

tar -xjf manta-1.6.0.centos6_x86_64.tar.bz2

export PATH=$INSTALL_DIR/manta-1.6.0.centos6_x86_64/bin:$PATH

echo export PATH=$INSTALL_DIR/manta-1.6.0.centos6_x86_64/bin:$PATH >> ~/.bashrc

wget https://github.com/Illumina/ExpansionHunterDenovo/releases/download/v0.9.0/ExpansionHunterDenovo-v0.9.0-linux_x86_64.tar.gz

tar -xzf ExpansionHunterDenovo-v0.9.0-linux_x86_64.tar.gz

export PATH=$INSTALL_DIR/ExpansionHunterDenovo-v0.9.0-linux_x86_64/:$PATH

echo export PATH=$INSTALL_DIR/ExpansionHunterDenovo-v0.9.0-linux_x86_64/:$PATH >> ~/.bashrc

git clone https://github.com/lgmgeo/AnnotSV

cd AnnotSV

make PREFIX=. install

make PREFIX=. install-human-annotation

export ANNOTSV=$INSTALL_DIR/AnnotSV/

cd $ANNOTSV/share/AnnotSV/Annotations_Human/Genes/GRCh37

wget https://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz

cd $INSTALL_DIR

git clone https://github.com/mobidic/knotAnnotSV

export PATH=$INSTALL_DIR/knotAnnotSV/:$PATH

echo export PATH=$INSTALL_DIR/knotAnnotSV/:$PATH >> ~/.bashrc

cd $DNASCAN_DIR

sed "s|path_reference = \"\"|path_reference = \"$DNASCAN_DIR\/hg19\/hg19.fa\"|" scripts/paths_configs.py > scripts/paths_configs.py_temp

sed "s|path_hisat_index = \"\"|path_hisat_index = \"$DNASCAN_DIR\/hg19\/hg19\"|" scripts/paths_configs.py_temp > scripts/paths_configs.py

sed "s|path_bwa_index = \"\"|path_bwa_index = \"$DNASCAN_DIR\/hg19\/hg19.fa\"|" scripts/paths_configs.py > scripts/paths_configs.py_temp

sed "s|path_annovar_db = \"\"|path_annovar_db = \"$INSTALL_DIR\/humandb\/\"|" scripts/paths_configs.py_temp > scripts/paths_configs.py

sed "s|dnascan_dir = \"\"|dnascan_dir = \"$DNASCAN_DIR\/\"|" scripts/paths_configs.py > scripts/paths_configs.py_temp

sed "s|path_manta = \"\"|path_manta = \"$INSTALL_DIR\/manta-1.6.0.centos6_x86_64\/bin\/\"|" scripts/paths_configs.py_temp > scripts/paths_configs.py

sed "s|path_expansionHunter_jsons = \"\"|path_expansionHunter_jsons = \"$DNASCAN_DIR\/repeats\/hg19\/\"|" scripts/paths_configs.py > scripts/paths_configs.py_temp

sed "s|path_to_db = \"\"|path_to_db = \"$DNASCAN_DIR\/db\/\"|" scripts/paths_configs.py_temp > scripts/paths_configs.py

sed "s|path_annotsv = \"\"|path_annotsv = \"$INSTALL_DIR\/AnnotSV\/\"|" scripts/paths_configs.py > scripts/paths_configs.py_temp

sed "s|path_knotannotsv = \"\"|path_knotannotsv = \"$INSTALL_DIR\/knotAnnotSV\/\"|" scripts/paths_configs.py_temp > scripts/paths_configs.py

sed "s|path_strelka = \"\"|path_strelka = \"$INSTALL_DIR\/strelka-2.9.10.centos6_x86_64\/bin\/\"|" scripts/paths_configs.py > scripts/paths_configs.py_temp

sed "s|path_delly_exclude_regions = \"\"|path_delly_exclude_regions = \"$DNASCAN_DIR\/db\/delly_hg19.excl.tsv\"|"  scripts/paths_configs.py_temp > scripts/paths_configs.py

sed "s|path_expansionHunterDenovo_dir = \"\"|path_expansionHunterDenovo_dir = \"$INSTALL_DIR\/ExpansionHunterDenovo-v0.9.0-linux_x86_64\/\"|" scripts/paths_configs.py > scripts/paths_configs.py_temp

mv scripts/paths_configs.py_temp scripts/paths_configs.py

chmod +x scripts/*

export PATH=$DNASCAN_DIR/scripts/:$PATH

echo export PATH=$DNASCAN_DIR/scripts/:$PATH >> ~/.bashrc

echo "###########################################IMPORTANT######################################################"
echo "Hisat2-build and bwa-index are still creating their indexes. Please wait untill they complete their task."
echo "You can check whether or not they are still running using the 'top' command"
echo "##########################################################################################################"
