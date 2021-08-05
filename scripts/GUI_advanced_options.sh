#!/bin/bash

#Usage: bash add_dependencies.sh $path_to_gene_list $path_to_bed_file $hisat_options $bwa_options $annotsv_options $melt_options $RG_ID $RG_LB $RG_SM $RG_PU $RG_PL

path_gene_list=$1

path_bed=$2

hisat_custom_options=$3

bwa_custom_options=$4

annotsv_custom_options=$5

melt_custom_options=$6

RG_ID=$7

RG_LB=$8

RG_SM=$9

RG_PU=$10

RG_PL=$11

sed "s|path_gene_list = \"\"|path_gene_list = \"$path_gene_list\"|" scripts/paths_configs.py > scripts/paths_configs.py_temp

sed "s|path_bed = \"\"|path_bed = \"$path_bed\"|" scripts/paths_configs.py_temp > scripts/paths_configs.py

sed "s|hisat_custom_options = \"\"|hisat_custom_options = \"$hisat_custom_options\"|" scripts/paths_configs.py > scripts/paths_configs.py_temp

sed "s|bwa_custom_options = \"\"|bwa_custom_options = \"$bwa_custom_options\"|" scripts/paths_configs.py_temp > scripts/paths_configs.py

sed "s|annotsv_custom_options = \"\"|annotsv_custom_options = \"$annotsv_custom_options\"|" scripts/paths_configs.py > scripts/paths_configs.py_temp

sed "s|melt_custom_options = \"\"|melt_custom_options = \"$melt_custom_options\"|" scripts/paths_configs.py_temp > scripts/paths_configs.py

sed "s|RG_ID = \"\"|RG_ID = \"$RG_ID\"|" scripts/paths_configs.py > scripts/paths_configs.py_temp

sed "s|RG_LB = \"\"|RG_LB = \"$RG_LB\"|" scripts/paths_configs.py_temp > scripts/paths_configs.py

sed "s|RG_SM = \"\"|RG_SM = \"$RG_SM\"|" scripts/paths_configs.py > scripts/paths_configs.py_temp

sed "s|RG_PU = \"\"|RG_PU = \"$RG_PU\"|" scripts/paths_configs.py_temp > scripts/paths_configs.py

sed "s|RG_PL = \"\"|RG_PL = \"$RG_PL\"|" scripts/paths_configs.py > scripts/paths_configs.py_temp

mv scripts/paths_configs.py_temp scripts/paths_configs.py

chmod +x scripts/*
