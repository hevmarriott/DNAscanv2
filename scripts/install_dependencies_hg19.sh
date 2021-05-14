#!/bin/bash

#Usage: bash install_dependencies.sh $path_to_setup_dir $path_to_DNASCAN_dir $path_to_ANNOVAR $path_to_gatk_download
#Example: bash install_dependencies.sh /home/local/ /home/DNAscan /home/annovar /home/gatk_download_dir

INSTALL_DIR=$1

DNASCAN_DIR=$2

ANNOVAR_DIR=$3

GATK_DOWNLOAD_DIR=$4

NUM_CPUS=$5

sudo apt-get update

sudo apt-get install -y vim

sudo apt-get install -y python2

sudo apt-get install -y perl

sudo apt-get install -y ttf-dejavu

sudo apt-get install -y wget bzip2 gzip

mkdir $INSTALL_DIR

mkdir $INSTALL_DIR/humandb

cd $DNASCAN_DIR

chmod +x $ANNOVAR_DIR/*

#nohup $ANNOVAR_DIR/annotate_variation.pl -buildver hg19 -downdb -webfrom annovar cadd $INSTALL_DIR/humandb/ &

$ANNOVAR_DIR/annotate_variation.pl -buildver hg19 -downdb -webfrom annovar refGene $INSTALL_DIR/humandb/

$ANNOVAR_DIR/annotate_variation.pl -buildver hg19 -downdb -webfrom annovar exac03 $INSTALL_DIR/humandb/

$ANNOVAR_DIR/annotate_variation.pl -buildver hg19 -downdb -webfrom annovar dbnsfp30a $INSTALL_DIR/humandb/

$ANNOVAR_DIR/annotate_variation.pl -buildver hg19 -downdb -webfrom annovar clinvar_20210123 $INSTALL_DIR/humandb/

$ANNOVAR_DIR/annotate_variation.pl -buildver hg19 -downdb -webfrom annovar avsnp147 $INSTALL_DIR/humandb/

cd $INSTALL_DIR

wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

chmod +x Miniconda3-latest-Linux-x86_64.sh

bash Miniconda3-latest-Linux-x86_64.sh -b -p $INSTALL_DIR/Miniconda3/

export PATH=$INSTALL_DIR/Miniconda3/bin:$PATH

echo export PATH=$INSTALL_DIR/Miniconda3/bin:$PATH >> ~/.bashrc

conda config --add channels conda-forge

conda config --add channels defaults

conda config --add channels r

conda config --add channels bioconda

conda install -y python=3.8.0

conda install -y pysam=0.16.0.1

conda install -y samtools=1.9

conda install -y freebayes=1.0.2

conda install -y bedtools=2.25.0

conda install -y vcftools=0.1.16

conda install -y bcftools=1.9

conda install -y hisat2=2.2.1

conda install -y bwa=0.7.17

conda install -y rtg-tools=3.12

conda install -y multiqc=1.10.1

conda install -y fastqc=0.11.9

conda install -y expansionhunter=3.2.2

conda install -y wham=1.8.0

conda install -y sambamba=0.7.1

conda install -y samblaster=0.1.26

cd $DNASCAN_DIR

mkdir hg19

cd hg19

wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz

tar -zxvf chromFa.tar.gz

for i in chr1.fa chr2.fa chr3.fa chr4.fa chr5.fa chr6.fa chr7.fa chr8.fa chr9.fa chr10.fa chr11.fa chr12.fa chr13.fa chr14.fa chr15.fa chr16.fa chr17.fa chr18.fa chr19.fa chr20.fa chr21.fa chr22.fa chrY.fa chrX.fa chrM.fa; do cat $i >> hg19.fa ; rm $i ; done

rm chr*

samtools faidx hg19.fa

nohup bwa index hg19.fa &

nohup hisat2-build -p $NUM_CPUS hg19.fa hg19 &

sudo apt-get update -qq

sudo apt-get install -y -qq bzip2 gcc g++ make python zlib1g-dev

cd $INSTALL_DIR

wget https://github.com/Illumina/manta/releases/download/v1.6.0/manta-1.6.0.centos6_x86_64.tar.bz2

tar -xjf manta-1.6.0.centos6_x86_64.tar.bz2

export PATH=$INSTALL_DIR/manta-1.6.0.centos6_x86_64/bin:$PATH

echo export PATH=$INSTALL_DIR/manta-1.6.0.centos6_x86_64/bin:$PATH >> ~/.bashrc

cd $DNASCAN_DIR

sed "s|path_reference = \"\"|path_reference = \"$DNASCAN_DIR\/hg19\/hg19.fa\"|" scripts/paths_configs.py > scripts/paths_configs.py_temp

sed "s|path_hisat_index = \"\"|path_hisat_index = \"$DNASCAN_DIR\/hg19\/hg19\"|" scripts/paths_configs.py_temp > scripts/paths_configs.py

sed "s|path_bwa_index = \"\"|path_bwa_index = \"$DNASCAN_DIR\/hg19\/hg19.fa\"|" scripts/paths_configs.py > scripts/paths_configs.py_temp

sed "s|path_annovar = \"\"|path_annovar = \"$ANNOVAR_DIR\/\"|" scripts/paths_configs.py_temp > scripts/paths_configs.py

sed "s|path_annovar_db = \"\"|path_annovar_db = \"$INSTALL_DIR\/humandb\/\"|" scripts/paths_configs.py > scripts/paths_configs.py_temp

sed "s|path_gatk = \"\"|path_gatk = \"$GATK_DOWNLOAD_DIR\/\"|" scripts/paths_configs.py_temp >  scripts/paths_configs.py

sed "s|dnascan_dir = \"\"|dnascan_dir = \"$DNASCAN_DIR\/\"|" scripts/paths_configs.py > scripts/paths_configs.py_temp

sed "s|path_manta = \"\"|path_manta = \"$INSTALL_DIR\/manta-1.6.0.centos6_x86_64\/bin\/\"|" scripts/paths_configs.py_temp > scripts/paths_configs.py

sed "s|path_samblaster = \"\"|path_samblaster = \"$INSTALL_DIR\/Miniconda3\/bin\/\"|" scripts/paths_configs.py > scripts/paths_configs.py_temp

sed "s|path_sambamba = \"\"|path_sambamba = \"$INSTALL_DIR\/Miniconda3\/bin\/\"|" scripts/paths_configs.py_temp > scripts/paths_configs.py

sed "s|path_whamg = \"\"|path_whamg = \"$INSTALL_DIR\/Miniconda3\/bin\/\"|" scripts/paths_configs.py > scripts/paths_configs.py_temp

sed "s|path_expansionHunter = \"\"|path_expansionHunter = \"$INSTALL_DIR\/Miniconda3\/bin\/\"|" scripts/paths_configs.py_temp > scripts/paths_configs.py

sed "s|path_fastqc = \"\"|path_fastqc = \"$INSTALL_DIR\/Miniconda3\/bin\/\"|" scripts/paths_configs.py > scripts/paths_configs.py_temp

sed "s|path_multiqc = \"\"|path_multiqc = \"$INSTALL_DIR\/Miniconda3\/bin\/\"|" scripts/paths_configs.py_temp > scripts/paths_configs.py

sed "s|path_rtg = \"\"|path_rtg = \"$INSTALL_DIR\/Miniconda3\/bin\/\"|" scripts/paths_configs.py > scripts/paths_configs.py_temp

sed "s|path_bwa = \"\"|path_bwa = \"$INSTALL_DIR\/Miniconda3\/bin\/\"|" scripts/paths_configs.py_temp > scripts/paths_configs.py

sed "s|path_hisat = \"\"|path_hisat = \"$INSTALL_DIR\/Miniconda3\/bin\/\"|" scripts/paths_configs.py > scripts/paths_configs.py_temp

sed "s|path_bcftools = \"\"|path_bcftools = \"$INSTALL_DIR\/Miniconda3\/bin\/\"|" scripts/paths_configs.py_temp > scripts/paths_configs.py

sed "s|path_samtools = \"\"|path_samtools = \"$INSTALL_DIR\/Miniconda3\/bin\/\"|" scripts/paths_configs.py > scripts/paths_configs.py_temp

sed "s|path_freebayes = \"\"|path_freebayes = \"$INSTALL_DIR\/Miniconda3\/bin\/\"|" scripts/paths_configs.py_temp > scripts/paths_configs.py

sed "s|path_java = \"\"|path_java = \"$INSTALL_DIR\/Miniconda3\/bin\/\"|" scripts/paths_configs.py > scripts/paths_configs.py_temp

sed "s|path_vcftools = \"\"|path_vcftools = \"$INSTALL_DIR\/Miniconda3\/bin\/\"|" scripts/paths_configs.py_temp > scripts/paths_configs.py

sed "s|path_tabix = \"\"|path_tabix = \"$INSTALL_DIR\/Miniconda3\/bin\/\"|" scripts/paths_configs.py > scripts/paths_configs.py_temp

sed "s|path_bedtools = \"\"|path_bedtools = \"$INSTALL_DIR\/Miniconda3\/bin\/\"|" scripts/paths_configs.py_temp > scripts/paths_configs.py

mv scripts/paths_configs.py_temp scripts/paths_configs.py

chmod +x scripts/*

export PATH=$DNASCAN_DIR/scripts/:$PATH

echo export PATH=$DNASCAN_DIR/scripts/:$PATH >> ~/.bashrc

echo "###########################################IMPORTANT######################################################"
echo "Hisat2-build and bwa-index are still creating their indexes. Please wait untill they complete their task."
echo "You can check whether or not they are still running using the 'top' command"
echo "##########################################################################################################"
