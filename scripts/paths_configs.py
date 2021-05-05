#!/usr/bin/env python3
#DNAscan uses several tools and files. DNAscan paths to software binaries always endwith "/" while paths to files always end with the filename.
#Instructions about how to get the needed files are commented below near the coresponding path line.

import os.path

#Configs 

RG_ID = "4"

RG_LB = "lib1"

RG_PL = "illumina"

RG_PU = "unit1"

RG_SM = "20"

tmp_dir = ""

num_cpu = "1"


#Paths

path_to_db = ""

dnascan_dir = ""

path_hisat = ""

path_samtools = ""

path_rtg = ""

path_multiqc = ""

path_fastqc = ""

path_sambamba = ""

path_freebayes = ""

path_annovar = ""

path_annovar_db = ""

#annovar_protocols and annovar_operations can be modified according to the wanted annovar databases. For instructions about
#how to download the annovar databases and to modify annovar_protocols and annovar_operations please visit 
#http://annovar.openbioinformatics.org/en/latest/ . 

annovar_protocols = "refGene,dbnsfp30a,clinvar_20210123,avsnp147,exac03"

annovar_operations = "g,f,f,f,f"

path_expansionHunter = ""

path_expansionHunter_jsons = ""

path_samblaster = ""

path_java = ""

path_vcftools = ""

path_scripts = "scripts/"

path_gatk = ""

path_bed = "data/test_data.bed"

path_gene_list = "data/list_of_genes.txt"

#path_bed="/users/k1513213/brc_scratch/indels_project/gene_list_positions_sorted_no_overlap.bed"

#hg19 can be downloaded from http://hgdownload.cse.ucsc.edu/goldenPath/hg19/
path_reference = ""

#hg19 index can be downloaded from ftp://ftp.ccb.jhu.edu/pub/data/hisat_indexes/hg19_hisat.tar.gz
#the index can be created running "./$path_hisat/hisat-build $path_reference index_base"
path_hisat_index = ""

path_manta = ""

path_whamg = ""

whamg_exclude_regions = "GL000207.1,GL000226.1,GL000229.1,GL000231.1,GL000210.1,GL000239.1,GL000235.1,GL000201.1,GL000247.1,GL000245.1,GL000197.1,GL000203.1,GL000246.1,GL000249.1,GL000196.1,GL000248.1,GL000244.1,GL000238.1,GL000202.1,GL000234.1,GL000232.1,GL000206.1,GL000240.1,GL000236.1,GL000241.1,GL000243.1,GL000242.1,GL000230.1,GL000237.1,GL000233.1,GL000204.1,GL000198.1,GL000208.1,GL000191.1,GL000227.1,GL000228.1,GL000214.1,GL000221.1,GL000209.1,GL000218.1,GL000220.1,GL000213.1,GL000211.1,GL000199.1,GL000217.1,GL000216.1,GL000215.1,GL000205.1,GL000219.1,GL000224.1,GL000223.1,GL000195.1,GL000212.1,GL000222.1,GL000200.1,GL000193.1,GL000194.1,GL000225.1,GL000192.1,NC_007605"

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

freebayes_custom_options = ""

GATK_HC_custom_options = ""
