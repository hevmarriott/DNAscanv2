#!/bin/bash
#Usage: bash GUI_advanced_options.sh  $path_to_bed_file $path_to_gene_list $hisat_options $bwa_options $annotsv_options $melt_options
path_gene_list=$2
path_bed=$1
hisat_custom_options=$3
bwa_custom_options=$4
annotsv_custom_options=$5
melt_custom_options=$6
sed "s|path_gene_list = \"\"|path_gene_list = \"$path_gene_list\"|" scripts/paths_configs.py > scripts/paths_configs.py_temp
sed "s|path_bed = \"\"|path_bed = \"$path_bed\"|" scripts/paths_configs.py_temp > scripts/paths_configs.py
sed "s|hisat_custom_options = \"\"|hisat_custom_options = \"$hisat_custom_options\"|" scripts/paths_configs.py > scripts/paths_configs.py_temp
sed "s|bwa_custom_options = \"\"|bwa_custom_options = \"$bwa_custom_options\"|" scripts/paths_configs.py_temp > scripts/paths_configs.py
sed "s|annotsv_custom_options = \"\"|annotsv_custom_options = \"$annotsv_custom_options\"|" scripts/paths_configs.py > scripts/paths_configs.py_temp
sed "s|melt_custom_options = \"\"|melt_custom_options = \"$melt_custom_options\"|" scripts/paths_configs.py_temp > scripts/paths_configs.py
chmod +x scripts/*
