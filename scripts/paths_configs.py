#!/usr/bin/env python3
#DNAscan uses several tools and files. DNAscan paths to software binaries always endwith "/" while paths to files always end with the filename.
#Instructions about how to get the needed files are commented below near the coresponding path line.

import os.path

#Configs

RG_ID = ""

RG_LB = ""

RG_PL = ""

RG_PU = ""

RG_SM = ""

tmp_dir = ""

num_cpu = "4"

RAM_GB = ""


#Paths

path_to_db = ""

dnascan_dir = ""

path_python2 = ""

path_hisat = ""

path_samtools = ""

path_rtg = ""

path_multiqc = ""

path_fastqc = ""

path_sambamba = ""

path_strelka = ""

path_annovar = ""

path_annovar_db = ""

path_annotsv = ""

path_knotannotsv = ""

#annovar_protocols and annovar_operations can be modified according to the wanted annovar databases. For instructions about
#how to download the annovar databases and to modify annovar_protocols and annovar_operations please visit
#http://annovar.openbioinformatics.org/en/latest/ .

annovar_protocols = "refGene,dbnsfp33a,clinvar_20210501,intervar_20180118,avsnp147,exac03,ALL.sites.2015_08,gnomad211_genome"

annovar_operations = "g,f,f,f,f,f,f,f"

path_expansionHunter = ""

path_expansionHunter_jsons = ""

path_expansionHunterDenovo_dir = ""

path_samblaster = ""

path_java = ""

path_SURVIVOR = ""

path_vcftools = ""

path_scripts = "scripts/"

#path_bed with test data ="data/test_data.bed"
path_bed = ""

#path_gene_list with test data="data/list_of_genes.txt"
path_gene_list = ""

#e.g. hg19 can be downloaded from http://hgdownload.cse.ucsc.edu/goldenPath/hg19/
path_reference = ""

#e.g. hg19 index can be downloaded from ftp://ftp.ccb.jhu.edu/pub/data/hisat_indexes/hg19_hisat.tar.gz
#e.g. the index can be created running "./$path_hisat/hisat-build $path_reference index_base"
path_hisat_index = ""

path_manta = ""

path_delly = ""

path_delly_exclude_regions = ""

path_melt = ""

path_bedtools = ""

path_tabix = ""

#A viral index for the whole NCBI database of viral complete genomes can be downloaded from XXX
#Otherwise the index can be created running "./$path_hisat/hisat-build $path_reference index_base"
path_virus_index = ""

#A bacterial index for the whole NCBI database of bacterial complete genomes can be downloaded from XXX
#Otherwise the index can be created running "./$path_hisat/hisat-build $path_reference index_base"
path_virus_bacteria_index = ""

#The index can be created running "./$path_hisat/hisat-build $path_custom_reference index_base"
path_custom_microbes_index = ""

path_bwa = ""

path_bcftools = ""

#hg19 index can be downloaded from XXX
#the index can be created running "./$path_bwa/bwa index $path_reference"
path_bwa_index = ""

#custom tool options


hisat_custom_options = ""

bwa_custom_options = ""

melt_custom_options = ""

annotsv_custom_options = ""
