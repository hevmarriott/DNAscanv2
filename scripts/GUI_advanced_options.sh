#!/bin/bash

#Usage: bash GUI_advanced_options.sh $path_to_gene_list $path_to_bed_file $hisat_options $bwa_options $annotsv_options $melt_options

args=()
[  "x$1" = xFALSE ] ||  args+=( --optflag1 "$1" )
[  "x$2" = xFALSE ] ||  args+=( --optflag2 "$2" )
[  "x$3" = xFALSE ] ||  args+=( --optflag3 "$3" )
[  "x$4" = xFALSE ] ||  args+=( --optflag4 "$4" )
[  "x$5" = xFALSE ] ||  args+=( --optflag5 "$5" )
[  "x$6" = xFALSE ] ||  args+=( --optflag6 "$6" )

GUI_advanced_options "${args[@]}"

path_gene_list=$1

path_bed=$2

hisat_custom_options=$3

bwa_custom_options=$4

annotsv_custom_options=$5

melt_custom_options=$6

sed -i "s|path_gene_list = \"\"|path_gene_list = \"$path_gene_list\"|" scripts/paths_configs.py

sed -i "s|path_bed = \"\"|path_bed = \"$path_bed\"|"  scripts/paths_configs.py

sed -i "s|hisat_custom_options = \"\"|hisat_custom_options = \"$hisat_custom_options\"|" scripts/paths_configs.py 

sed -i "s|bwa_custom_options = \"\"|bwa_custom_options = \"$bwa_custom_options\"|" scripts/paths_configs.py

sed -i "s|annotsv_custom_options = \"\"|annotsv_custom_options = \"$annotsv_custom_options\"|" scripts/paths_configs.py 

sed -i "s|melt_custom_options = \"\"|melt_custom_options = \"$melt_custom_options\"|" scripts/paths_configs.py 

chmod +x scripts/*
