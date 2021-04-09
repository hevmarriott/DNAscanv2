#!/bin/bash

#Usage: bash install_annovar_and_gatk.sh $path_to_setup_dir $path_to_DNASCAN_dir $path_to_ANNOVAR $path_to_gatk_download
#Example: bash install_annovar_and_gatk.sh /home/local/ /home/DNA-NGS_scan /home/annovar /home/gatk_download_dirINSTALL_DIR=$1

INSTALL_DIR=$1

DNASCAN_DIR=$2

ANNOVAR_DIR=$3

GATK_DOWNLOAD_DIR=$4

mkdir $INSTALL_DIR/humandb

cd $DNASCAN_DIR

chmod +x $ANNOVAR_DIR/*

nohup $ANNOVAR_DIR/annotate_variation.pl -buildver hg19 -downdb -webfrom annovar cadd $INSTALL_DIR/humandb/ &

$ANNOVAR_DIR/annotate_variation.pl -buildver hg19 -downdb -webfrom annovar refGene $INSTALL_DIR/humandb/

$ANNOVAR_DIR/annotate_variation.pl -buildver hg19 -downdb -webfrom annovar exac03 $INSTALL_DIR/humandb/

$ANNOVAR_DIR/annotate_variation.pl -buildver hg19 -downdb -webfrom annovar dbnsfp30a $INSTALL_DIR/humandb/

$ANNOVAR_DIR/annotate_variation.pl -buildver hg19 -downdb -webfrom annovar clinvar_20210123 $INSTALL_DIR/humandb/

$ANNOVAR_DIR/annotate_variation.pl -buildver hg19 -downdb -webfrom annovar avsnp147 $INSTALL_DIR/humandb/

gatk-register $GATK_DOWNLOAD_DIR/GenomeAnalysisTK-3.8.tar.bz2 

sed "s|path_annovar = \"\"|path_annovar = \"$ANNOVAR_DIR\/\"|" scripts/paths_configs.py_temp > scripts/paths_configs.py

sed "s|path_gatk = \"\"|path_gatk = \"$INSTALL_DIR\/Miniconda3\/opt\/gatk-4.1.9.0\/\"|" scripts/paths_configs.py_temp >  scripts/paths_configs.py




echo "###########################################IMPORTANT######################################################"
echo "Annovar is still downloading the CADD database (wget). Please wait until it completes its tasks."
echo "You can check whether or not they are still running using the 'top' command"
echo "##########################################################################################################"
