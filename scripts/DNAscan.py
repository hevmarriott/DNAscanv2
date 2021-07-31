#!/usr/bin/env python3

################################################################
# Program: DNAscan
# Version 0.1
# Author: Alfredo Iacoangeli (alfredo.iacoangeli@kcl.ac.uk)
#################################################################

################################################################
# Script structure:
# 1. Import python modules
# 2. Define paths_configs variables
# 3. Define options from command line
# 4. Parse options from command line
# 5. Create working dir tree
# 6. Bed splitting
# 7. Remove duplicates
# 8. Alignment
#   8.1 Aligns paired end reads
#       8.1.1 Fast mode alignment
#       8.1.2 Normal and intensive mode alignment
#   8.2 Aligns single end reads
#       8.2.1 Fast mode alignment
#       8.2.2 Normal and intensive mode alignment
# 9. If input file is a sam file, it converts it into bam
# 10. Variant (snv and indel) calling
#   10.1 Strelka snv and indel calling
#   10.2 Strelka SNV and indel calling using indel sites (only performed in intensive mode)
#       10.2.1 Identification of potential indel sites
#       10.2.2 Indel calling on selected positions (from previous step 10.1.1)
#       10.2.3 Merging of snv and indel files
# 11. Perform variant hard filtering
# 12. Perform known expansions search with ExpansionHunter
# 13. Structural Variant calling with Manta (normal and intensive mode)
#   13.1 Manta only (fast mode)
#   13.2 Delly additionally calls inversion and deletion variants (in normal mode) and general structural variant classes (intensive mode)
#   13.3 SV calls are merged to create union callset (normal and intensive modes)
# 14. Transposable element insertion detection with MELT
# 15. Annotation with Annovar -  with optional missense variant prioritisation according to ACMG guidelines (intervar_20180118 database)
#   15.1 AnnotSV for structural variant annotation and prioritisation
# 16. Microbes screening
#   16.1 Extract non human reads
#   16.2 Identifies present non-human microbes
#       16.2.1 Identifies present viruses
#       16.2.2 Generates virus report
#       16.2.3 Identifies present bacteria
#       16.2.4 Generates bacteria report
#       16.2.5 Identifies present user-selected microbes
#       16.2.6 Generates user-selected microbes report
# 17. Alignment report generation ( samtools flagstat and stats )
# 18. Sequencing data report generation ( fastqc )
# 19. Snv and indel calling report generation ( bcftools stats )
# 20. Html report generation ( Multiqc )
# 21. Annotated variants results report generation
#   21.1 knotAnnotSV annotated structural variants report generation
# 22. Starting iobio services
#################################################################

# 1. Import necessary python modules

import argparse
import sys
import string
import os
import os.path
import re
import csv
import paths_configs
import subprocess
import pysam
import gzip
import loadInNS
import parsing_file

from argparse import RawTextHelpFormatter

# 2. Define paths_configs variables from paths_configs.py

var = loadInNS.load('scripts/paths_configs.py')
locals().update(var)

# 3. Define options variables from command line

parser = parsing_file.create_parser()

requiredNamed = parser.add_argument_group('required named arguments')

requiredNamed.add_argument(
    '-out',
    required=True,
    action="store",
    dest="out",
    help=
    'path to the output folder. It has to end in /" e.g. /home/user/local/test_folder/'
)

requiredNamed.add_argument(
    '-in',
    required=True,
    action="store",
    dest="input_file",
    help='input file [string]')

# 4. Parse command line options

args = parser.parse_args()
mode = args.mode
alsgenescanner = args.alsgenescanner
exome = args.exome
format = args.format
paired = args.paired
reference = args.reference
input_file = args.input_file
input_file2 = args.input_file2
expansion = args.expansion
out = args.out + '/'
vcf = args.vcf
SV = args.SV
MEI = args.MEI
BED = args.BED
bacteria = args.bacteria
virus = args.virus
alignment = args.alignment
iobio = args.iobio
variantcalling = args.variantcalling
annotation = args.annotation
results_report = args.results_report
alignment_report = args.alignment_report
sequencing_report = args.sequencing_report
calls_report = args.calls_report
sample_name = args.sample_name
rm_dup = args.rm_dup
debug = args.debug
filter_string = args.filter_string
custom_microbes = args.custom_microbes
RG = args.RG
tmp = args.tmp
ref_file = args.ref_file
dnascan_main_dir = args.dnascan_main_dir

if ref_file:
    path_reference = ref_file

if dnascan_main_dir:
    dnascan_dir = dnascan_main_dir

variant_results_file = ""

if format == "vcf":
    variant_results_file = input_file

    if not vcf:
        vcf = variant_results_file


print("\n################################################")

# 5. Create working dir tree

print(
    "\nCreating working directory tree...\n"
)

os.system(
    "mkdir %s ; mkdir  %slogs ; mkdir  %sreports ; mkdir  %sresults ; mkdir %stmp"
    % (out, out, out, out, out))

print(
                "############DNAscan Options############ \n\n DNAscan is running an analysis with the following specifics:\n"
            )

options_log = open('%s/logs/options.log' %(out), 'w')

for arg in vars(args):
        print(arg,':    ', getattr(args,arg))
        options_log.write(arg + ':    ' + str(getattr(args,arg)) + "\n")
options_log.close()

print('\nOptions saved onto %s/logs/options.log \n' %(out))

# Function for checking output result files
def is_variant_file_OK(file, t, s):
    if os.path.isfile(file) == True:
        if os.path.getsize(file) != 0:
            if t == "bam":
                with pysam.AlignmentFile(file, 'rb') as f:
                    if f.nreferences == 0 and f.mapped == 0 and f.unmapped == 0:
                        sys.exit("\nWARNING: %s does not contain data, therefore DNAscan will now terminate.\n" % file)
                    else:
                        print("\n%s has sufficient data for DNAscan to continue...\n" % file)
                f.close()
            elif t == "Vcf":
                with gzip.open(file, 'r') as f:
                    if any(not line.startswith(b"#") for line in f):
                        print("\n%s has sufficient data for DNAscan to continue...\n" % file)
                    else:
                        if s == "variantcalling":
                            sys.exit("\nWARNING: %s only contains the header and no data, therefore DNAscan will now terminate.\n" % file)
                        else:
                            print("\nWARNING: %s only contains the header and no data.\n" % file)
                f.close()
        else:
            sys.exit("WARNING: %s is empty - DNAscan will now terminate.\n" % file)
    else:
        sys.exit("WARNING: %s does not exist - DNAscan will now terminate.\n" % file)

# 6. Bed splitting: splitting the analysis region into subsets of equal length to distribute the work across the available threads.
# To do this DNAscan uses a bed file.
# If bed file is not provided, it generates one starting from the
# reference genome index. The pipeline uses a bed file to split the
# analyses in several subprocesses

if alsgenescanner:

    alignment = True
    variantcalling = True
    annotation = True
    #SV = True
    #expansion = True
    BED = True
    path_expansionHunter_jsons = '%s/repeats/%s/' % (dnascan_dir, reference)
    path_gene_list = ""
    annovar_operations = "g,f,f,f"
    annovar_protocols = "refGene,dbnsfp33a,clinvar_20210501,intervar_20180118"

    if reference == "hg19":
        path_bed = '%s/als_gene_scanner_hg19.bed' % (path_to_db)
    else:
        path_bed = '%s/als_gene_scanner_hg38.bed' % (path_to_db)

# Y. adapt DB to reference

if reference == "grch37" or  reference == "grch38" :
    if reference == "grch37":
        ref_hg = "hg19"
    else:
        ref_hg = "hg38"
    print(
        "\nAdapting the database provided (exome, gene_db or gene_list) to the reference genome...\n"
    )

    os.system("zcat %s/exome_%s.bed.gz | sed 's/chr//g' | bgzip -c > %s/exome_%s.bed.gz" %(path_to_db,ref_hg,path_to_db,reference) )
    os.system("zcat %s/%s_gene_db.txt.gz | sed 's/chr//g' | bgzip -c > %s/%s_gene_db.txt.gz" %(path_to_db,ref_hg,path_to_db,reference) )
    os.system("cp %s/%s_gene_names.txt.gz  %s/%s_gene_names.txt.gz" %(path_to_db,ref_hg,path_to_db,reference) )

    if annotation == True:

        print(
                "\n\nWARNING: The annotation step can only be performed using hg19 and hg38. Unfortunately Annovar cannot be used with grch37 and grch38. DNAscan will not perform annotation.\n\n"
            )

        annotation = False

if BED or path_gene_list:
    if len(path_bed) != 0:
        # splitting the analysis region into subsets of equal length to
        # distribute the work across the available threads.
        if len(path_gene_list) != 0:
            print(
                "\n\nWARNING: Both a bed file and a list of genes were provided. DNAscan will ignore the list of genes.\n\n"
            )
        os.system(
            "awk \'{i=$2; while (i < $3) {print $1\"\t\"i\"\t\"i+1 ;  i++}}\' %s > %stmp/tmp.bed"
            % (path_bed, out))
        os.system(
            "split -d -l `wc -l %stmp/tmp.bed | awk '{if ($1/%s > int($1/%s)) print int($1/%s)+1; else print int($1/%s)}'` %stmp/tmp.bed %stmp/"
            % (out, num_cpu, num_cpu, num_cpu, num_cpu, out, out))
        os.system("rm %stmp/tmp.bed" % (out))

        i = 0
        zero = "0"

        while i < int(num_cpu):
            if i > 9:
                zero = ""

            os.system("%sbedtools merge -i %stmp/%s%s > %stemp%s.bed" %
                      (path_bedtools, out, zero, str(i), out, str(int(i) + 1)))
            #os.system("rm %stmp/%s%s" %(out,zero,str(i)))

            i += 1

    else:
        if len(path_gene_list) != 0:
            os.system(
                "zgrep -iwf %s %s%s_gene_names.txt.gz | awk '{print $2}' > %smatched_genes.txt"
                % (path_gene_list, path_to_db, reference, out))
            os.system(
                "zgrep -viwf %smatched_genes.txt %s  > %sunmatched_genes.txt" %
                (out, path_gene_list, out))

            if os.stat("%sunmatched_genes.txt" % (out)).st_size != 0:
                print(
                    "\n\nWARNING: some genes provided in the gene list were not found, please check which ones in %sunmatched_genes.txt "
                    % (out))
            os.system(
                "zgrep -iwf %s %s%s_gene_names.txt.gz | awk '{print $1}' > %smatched_genes_codes.txt"
                % (path_gene_list, path_to_db, reference, out))
            os.system(
                "zgrep -wf %smatched_genes_codes.txt %s%s_gene_db.txt | awk '{i=1; while (i<= int($8)) {n=split($9,a,/,/);n=split($10,b,/,/); print $2\"\t\"a[i]\"\t\"b[i]; i+=1}}' > %scustom_tmp.bed "
                % (out, path_to_db, reference, out))
            os.system(
                "%sbedtools sort -i  %scustom_tmp.bed> %scustom_sorted.bed" %
                (path_bedtools, out, out))
            os.system("%sbedtools merge -i  %scustom_sorted.bed> %scustom.bed"
                      % (path_bedtools, out, out))
            os.system("rm %scustom_sorted.bed %scustom_tmp.bed" % (out, out))
            os.system(
                "awk \'{i=$2; while (i < $3) {print $1\"\t\"i\"\t\"i+1 ;  i++}}\' %scustom.bed > %stmp/tmp.bed"
                % (out, out))
            os.system(
                "split -d -l `wc -l %stmp/tmp.bed | awk '{if ($1/%s > int($1/%s)) print int($1/%s)+1; else print int($1/%s)}'` %stmp/tmp.bed %stmp/"
                % (out, num_cpu, num_cpu, num_cpu, num_cpu, out, out))
            #os.system("rm %stmp/tmp.bed" %( out))

            i = 0
            zero = "0"

            while i < int(num_cpu):
                if i > 9:
                    zero = ""

                os.system(
                    "%sbedtools merge -i %stmp/%s%s > %stemp%s.bed" %
                    (path_bedtools, out, zero, str(i), out, str(int(i) + 1)))
                # os.system("rm %stmp/%s%s" %(out,zero,str(i)))

                i += 1

        else:
            sys.exit(
                '\n\nERROR: the BED flag was used but neither a bed file nor a gene list was provided. Either specify a bed file in paths_configs.py or remove the BED flag so DNAscan can create a BED file from the reference genome.\n\n'
            )

#    os.system(
#        "%sbedtools makewindows -n %s -i winnum -b %s | awk \'{print $1\"\t\"$2\"\t\"$3 >> \"%stemp\"$4\".bed\"}\'" %
#        (path_bedtools, num_cpu, path_bed, out))
else:
    # If bed file is not provided, it generates one starting from the
    # reference genome index or the default exome bed. The pipeline uses a bed file to split the
    # analyses in several subprocesses
    if exome == "True":
        os.system("zcat %sdb/exome_%s.bed.gz > %stmp/exome_%s.bed" %
                  (dnascan_dir, reference, out, reference))
        path_bed = "%stmp/exome_%s.bed" % (out, reference)

        if len(path_bed) != 0:
            # splitting the analysis region into subsets of equal length to
            # distribute the work across the available threads.
            os.system(
                "awk \'{i=$2; while (i < $3) {print $1\"\t\"i\"\t\"i+1 ;  i++}}\' %s > %stmp/tmp.bed"
                % (path_bed, out))
            os.system(
                "split -d -l `wc -l %stmp/tmp.bed | awk '{if ($1/%s > int($1/%s)) print int($1/%s)+1; else print int($1/%s)}'` %stmp/tmp.bed %stmp/"
                % (out, num_cpu, num_cpu, num_cpu, num_cpu, out, out))
            os.system("rm %stmp/tmp.bed" % (out))

            i = 0
            zero = "0"

            while i < int(num_cpu):
                if i > 9:
                    zero = ""

                os.system(
                    "%sbedtools merge -i %stmp/%s%s > %stemp%s.bed" %
                    (path_bedtools, out, zero, str(i), out, str(int(i) + 1)))
                os.system("rm %stmp/%s%s" % (out, zero, str(i)))

                i += 1

    else:
        os.system("cat %s.fai | awk '{print $1\"\t0\t\"$2}' > %sreference.bed"
                  % (path_reference, out))

        path_bed = "%sreference.bed" % (out)

        os.system(
            "%sbedtools makewindows -n %s -i winnum -b %s | awk \'{print $1\"\t\"$2\"\t\"$3 >> \"%stemp\"$4\".bed\"}\'"
            % (path_bedtools, num_cpu, path_bed, out))

    BED = True

# 7. Remove duplicates command line.
# The output from HISAT2 and BWA is piped into $samblaster_cmq during the alignment
if rm_dup == "True":
    print(
        "\nDuplicates will be removed after alignment.\n"
    )
    if exome == "True":
        samblaster_cmq = "%ssamblaster --ignoreUnmated |" % (path_samblaster)
    else:
        samblaster_cmq = "%ssamblaster |" % (path_samblaster)
        samblaster_bwa = "%ssamblaster --ignoreUnmated |" % (path_samblaster)

else:
    samblaster_cmq = ""

if tmp == True :
    tmp_dir = tmp_dir
else:
    tmp_dir = "/tmp"


print("\nTemporary files will be held at this directory: %s\n" % (tmp_dir))
# 8. Alignment
# Performs alignment if input sequencing data is in fastq format

if alignment:
    if format == "fastq" and "alignment.log" not in os.listdir(out + "logs"):

        # 8.1 Aligns paired end reads
        if paired == "1":

            # 8.1.1 Fast mode uses HISAT2 only to align all reads
            if mode == "fast":
                print(
                    "\nPerforming paired read alignment with HISAT2 in fast mode...\n"
                )

                os.system(
                    "%shisat2 %s --no-spliced-alignment -p %s -x %s -1 %s -2 %s | %s %ssamtools view -@ %s -Sb -  | %ssambamba sort -t %s  --tmpdir=%s -o %ssorted.bam  /dev/stdin"
                    % (path_hisat,hisat_custom_options, num_cpu, path_hisat_index, input_file,
                       input_file2, samblaster_cmq, path_samtools, num_cpu,
                       path_sambamba, num_cpu, tmp_dir, out))

                bam_file = "%ssorted.bam" % (out)

                is_variant_file_OK(bam_file, "bam", "alignment")

                os.system("touch  %slogs/alignment.log" % (out))

                print(
                    "\nCompleted paired read alignment with HISAT2.\n"
                )

            # 8.1.2 Normal and intensive modes use HISAT2 to align all reads,
            # then soft-clipped and unaligned reads are realigned with BWA mem
            if mode == "normal" or mode == "intensive":
                # Intensive mode uses GATK haplotype caller which does not
                # support read group missing anymore. Default values for RG id,
                # lb, pl, pu, and sm are defined in paths_configs.py. Please change
                # them as needed. Also, whamg structural variant calling requires RG tags in the input bam file/s.
                if mode == "intensive":
                    RG = True
                if RG:
                    rg_option_hisat2 = " --rg-id %s --rg LB:%s --rg PL:%s  --rg PU:%s --rg SM:%s " % (
                        RG_ID, RG_LB, RG_PL, RG_PU, RG_SM)
                    rg_option_bwa = " -R '@RG\\tID:%s\\tLB:%s\\tPL:%s\\tPU:%s\\tSM:%s' " % (
                        RG_ID, RG_LB, RG_PL, RG_PU, RG_SM)
                else:
                    rg_option_hisat2 = ""
                    rg_option_bwa = ""

                print(
                    "\nPerforming paired read alignment with HISAT2 in normal/intensive mode...\n"
                )
                print(
                    "%shisat2 %s %s  --no-softclip --no-spliced-alignment -p %s -x %s -1 %s -2 %s | %s %ssamtools view -@ %s -Sb -  | %ssambamba sort -t %s --tmpdir=%s -o %ssorted.bam /dev/stdin; %ssamtools index -@ %s %ssorted.bam"
                    % (path_hisat, hisat_custom_options, rg_option_hisat2, num_cpu, path_hisat_index,
                       input_file, input_file2, samblaster_cmq, path_samtools,
                       num_cpu, path_sambamba, num_cpu, tmp_dir, out, path_samtools,
                       num_cpu, out))
                os.system(
                    "%shisat2 %s %s  --no-softclip --no-spliced-alignment -p %s -x %s -1 %s -2 %s | %s %ssamtools view -@ %s -Sb -  | %ssambamba sort -t %s --tmpdir=%s -o %ssorted.bam /dev/stdin; %ssamtools index -@ %s %ssorted.bam"
                    % (path_hisat, hisat_custom_options, rg_option_hisat2, num_cpu, path_hisat_index,
                       input_file, input_file2, samblaster_cmq, path_samtools,
                       num_cpu, path_sambamba, num_cpu, tmp_dir,out, path_samtools,
                       num_cpu, out))
                print(
                    "%ssamtools view -@ %s -bhf 4 %ssorted.bam | samtools bam2fq - > %sunaligned_reads.fq"
                    % (path_samtools, num_cpu, out, out))
                os.system(
                    "%ssamtools view -@ %s -bhf 4 %ssorted.bam | samtools bam2fq - > %sunaligned_reads.fq"
                    % (path_samtools, num_cpu, out, out))
                print(
                    "\nPerforming paired read alignment of soft-clipped and unaligned HISAT2 reads with BWA-MEM in normal/intensive mode...\n"
                )
                print(
                    "%sbwa mem %s %s -t %s %s %sunaligned_reads.fq | %s %ssamtools view -@ %s -Sb -  | %ssambamba sort -t %s --tmpdir=%s -o %ssorted_bwa.bam  /dev/stdin ; %ssamtools index -@ %s %ssorted_bwa.bam "
                    % (path_bwa, bwa_custom_options, rg_option_bwa, num_cpu, path_bwa_index, out,
                       samblaster_bwa, path_samtools, num_cpu, path_sambamba,
                       num_cpu,tmp_dir, out, path_samtools, num_cpu, out))
                os.system(
                    "%sbwa mem %s %s -t %s %s %sunaligned_reads.fq | %s %ssamtools view -@ %s -Sb -  | %ssambamba sort -t %s --tmpdir=%s -o %ssorted_bwa.bam  /dev/stdin ; %ssamtools index -@ %s %ssorted_bwa.bam "
                    % (path_bwa, bwa_custom_options, rg_option_bwa, num_cpu, path_bwa_index, out,
                       samblaster_bwa, path_samtools, num_cpu, path_sambamba,
                       num_cpu, tmp_dir, out, path_samtools, num_cpu, out))
                os.system("%ssamtools view -H %ssorted.bam > %sheader.txt" %
                          (path_samtools, out, out))
                print("\nMerging HISAT2 and BWA-MEM aligned reads...\n")
                os.system(
                    "%ssamtools merge -c -@ %s -f -h %sheader.txt %ssorted_merged.bam %ssorted.bam  %ssorted_bwa.bam"
                    % (path_samtools, num_cpu, out, out, out, out))

                if not debug:
                    os.system(
                        "rm %ssorted.bam*  %ssorted_bwa.bam* %sunaligned_reads.fq %sheader.txt"
                        % (out, out, out, out))

                os.system("%ssamtools index -@ %s %ssorted_merged.bam" %
                          (path_samtools, num_cpu, out))

                bam_file = "%ssorted_merged.bam" % (out)

                is_variant_file_OK(bam_file, "bam", "alignment")

                os.system("touch  %slogs/alignment.log" % (out))

                print(
                    "\nCompleted paired read alignment with HISAT2 and BWA-MEM.\n"
                )

        # 8.2 Aligns single end reads
        if paired == "0":

            # 8.2.1 Fast mode uses HISAT2 only to align all reads
            if mode == "fast":
                print(
                    "\nPerforming single-end read alignment with HISAT2 in fast mode...\n"
                )
                os.system(
                    "%shisat2 %s --no-spliced-alignment -p %s -x %s -U %s | %s %ssamtools view -@ %s -Sb -  | %ssambamba sort -t %s --tmpdir=%s -o %ssorted.bam /dev/stdin"
                    % (path_hisat, hisat_custom_options, num_cpu, path_hisat_index, input_file,
                       samblaster_cmq, path_samtools, num_cpu, tmp_dir, path_sambamba,
                       num_cpu, out))

                bam_file = "%ssorted.bam" % (out)

                is_variant_file_OK(bam_file, "bam", "alignment")

                os.system("touch  %salignment.log" % (out))

                print(
                    "\nCompleted single-end read alignment with HISAT2.\n"
                )

            # 8.2.2 Normal and intensive modes use HISAT2 to align all reads,
            # then soft-clipped and unaligned reads are realigned with BWA mem
            if mode == "normal" or mode == "intensive":
                # Intensive mode uses GATK haplotype caller which does not
                # support read group missing anymore. Default values for RG id,
                # lb, pl, pu, and sm are defined in paths_configs.py. Please change
                # them as needed.
                if mode == "intensive":
                    RG = True

                if RG:
                    rg_option_hisat2 = " --rg-id %s --rg LB:%s --rg PL:%s  --rg PU:%s --rg SM:%s " % (
                        RG_ID, RG_LB, RG_PL, RG_PU, RG_SM)
                    rg_option_bwa = " -R '@RG\\tID:%s\\tLB:%s\\tPL:%s\\tPU:%s\\tSM:%s' " % (
                        RG_ID, RG_LB, RG_PL, RG_PU, RG_SM)
                else:
                    rg_option_hisat2 = ""
                    rg_option_bwa = ""

                print(
                    "\nPerforming single-end read alignment with HISAT2 in normal/intensive mode...\n"
                )

                os.system(
                    "%shisat2 %s --no-softclip --no-spliced-alignment -p %s -x %s -U %s | %s %ssamtools view -Sb -  | %ssambamba sort -t %s --tmpdir=%s -o %ssorted.bam /dev/stdin; %ssamtools index -@ %s %ssorted.bam"
                    % (path_hisat, hisat_custom_options, rg_option_hisat2, num_cpu, path_hisat_index,
                       input_file, samblaster_cmq, path_samtools,
                       path_sambamba, num_cpu,tmp_dir, out, path_samtools, num_cpu,
                       out))
                os.system(
                    "%ssamtools view -@ %s -bhf 4 %ssorted.bam | samtools bam2fq - > %sunaligned_reads.fq"
                    % (path_samtools, num_cpu, out, out))
                print(
                    "\nPerforming single-end read alignment of soft-clipped and unaligned HISAT2 reads with BWA-MEM in normal/intensive mode...\n"
                )
                os.system(
                    "%sbwa mem %s %s -t %s %s %sunaligned_reads.fq | %s %ssamtools view -Sb -  | %ssambamba sort -t %s --tmpdir=%s -o %ssorted_bwa.bam /dev/stdin; %ssamtools index -@ %s %ssorted_bwa.bam "
                    % (path_bwa, bwa_custom_options, rg_option_bwa, num_cpu, path_bwa_index, out,
                       samblaster_cmq, path_samtools, path_sambamba, num_cpu, tmp_dir,
                       out, path_samtools, num_cpu, out))
                os.system("%ssamtools view -H %ssorted.bam > %sheader.txt" %
                          (path_samtools, out, out))
                print("\nMerging HISAT2 and BWA-MEM aligned reads...\n")
                os.system(
                    "%ssamtools merge -c -@ %s -f -h %sheader.txt %ssorted_merged.bam %ssorted.bam  %ssorted_bwa.bam"
                    % (path_samtools, num_cpu, out, out, out, out))
                os.system("%ssamtools index -@ %s %ssorted_merged.bam " %
                          (path_samtools, num_cpu, out))

                if not debug:
                    os.system(
                        "rm %ssorted.bam*  %ssorted_bwa.bam* %sunaligned_reads.fq "
                        % (out, out, out))

                bam_file = "%ssorted_merged.bam" % (out)

                is_variant_file_OK(bam_file, "bam", "alignment")

                os.system("touch  %slogs/alignment.log" % (out))

                print(
                    "\nCompleted single-end read alignment with HISAT2 and BWA-MEM.\n"
                )

    else:
        if format != "fastq":
            print(
                "WARNING: Fastq format input data is required if you want to perform the alignment stage\n"
            )

        if "alignment.log" in os.listdir(out + "logs"):
            print(
                "WARNING: The presence of alignment.log in logs is telling you that the alignment was already peformed, please remove alignment.log if you wish to perform this stage anyway\n"
            )

            if mode == "normal" or mode == "intensive":
                bam_file = "%ssorted_merged.bam" % (out)

            if mode == "fast":
                bam_file = "%ssorted.bam" % (out)

# 9. Convert input sam file into bam

if format == "sam" and "sam2bam.log" not in os.listdir(out + "logs"):
    print(
        "\nConverting input sam file into bam format...\n"
    )
    os.system("%ssamtools view -Sb %s  > %ssorted.bam" % (path_samtools,
                                                          input_file, out))
    bam_file = "%ssorted.bam" % (out)

    is_variant_file_OK(bam_file, "bam", "convert")

    os.system("touch  %slogs/sam2bam.log" % (out))

if format == "bam":
    bam_file = "%s" % (input_file)

# 10.Variant (snv and indel) calling

if variantcalling:
    if vcf:
        print(
            "WARNING: Using input vcf as variant file. Do not provide vcf file if you wish to perform variant calling\n"
        )

        variant_results_file = vcf

    else:

        if "VC.log" in os.listdir(out + "logs"):
            print(
                "WARNING: The presence of VC.log in logs is telling you that the variant calling was already peformed, please remove VC.log if you wish to perform this stage anyway\n"
            )

            variant_results_file = "%s%s_sorted.vcf.gz" % (out, sample_name)

        else:
            #10.1 Strelka snv and indel calling (in fast and normal mode)
            if "VC_strelka" not in os.listdir(out + "logs"):

                if paired == "1":
                    print("\nSNVs and indels are being called with Strelka...\n")
                    
                    if BED == "True":
                        os.system("bgzip -c %s  > %s/temp.bed.gz" % (path_bed, out))
                        os.system(
                        "%ssortBed -i %s/temp.bed.gz | bgzip -c > %s/sorted.bed.gz" %
                        (path_bedtools, out, out))

                        os.system("%stabix -p bed %s/sorted.bed.gz" % (path_tabix, out))
                        os.system("mkdir %sstrelka" % (out))
                        os.system(
                        "%sconfigureStrelkaGermlineWorkflow.py --bam %s --referenceFasta %s --runDir %sstrelka --callRegions %s/sorted.bed.gz"
                        % (path_strelka, bam_file, path_reference, out, out))

                    if exome:
                        os.system("mkdir %sstrelka" % (out))
                        os.system("%sconfigureStrelkaGermlineWorkflow.py --bam %s --referenceFasta %s --runDir %sstrelka --exome" %
                        (path_strelka, bam_file, path_reference, out))
                        
                    else:
                        os.system("mkdir %sstrelka" % (out))
                        os.system(
                        "%sconfigureStrelkaGermlineWorkflow.py --bam %s --referenceFasta %s --runDir %sstrelka"
                        % (path_strelka, bam_file, path_reference, out))

                    os.system("%sstrelka/runWorkflow.py -j %s -m local" % (out, num_cpu))
                    os.system(
                    "mv %s/strelka/results/variants/genome.S1.vcf.gz  %s/results/%s_strelka.vcf.gz"
                    % (out, out, sample_name))
                    os.system(
                    "mv %s/strelka/results/variants/genome.S1.vcf.gz.tbi  %s/results/%s_strelka.vcf.gz.tbi"
                    % (out, out, sample_name))

                    if mode != "intensive":

                        os.system("mv %s/results/%s_strelka.vcf.gz %s%s_sorted.vcf.gz" % (out, sample_name, out, sample_name))
                        os.system ("tabix -p bcf %s%s_sorted.vcf.gz" % (out, sample_name))

                        variant_results_file = "%s%s_sorted.vcf.gz" % (out, sample_name)

                        is_variant_file_OK(variant_results_file, "Vcf", "variantcalling")

                    os.system("touch  %slogs/VC_strelka.log" % (out))

                    print('\nSNV and indel calling with Strelka is complete.\n')

                    if mode == "intensive":
                        # 10.3 Strelka indel calling (only performed in intensive mode)
                        # In intensive mode, DNAscan calls snvs with Freebayes and indels
                        # with Strelka. Strelka is used only on those positions of the
                        # genome for which the alignment stage identifies one insertion or
                        # deletion in at least on read.
                        # 10.3.1 identification of potential indel sites
                        # Samtools mpileup is used to identify those positions of the genome for which the alignment stage identifies one insertion or deletion in at least on read.
                        # Strelka is used to call indels only on those positions
                        print(
                            "\nIdentifying potential indel sites with samtools mpileup...\n"
                        )

                        counter = 1
                        ps = []

                        while counter < int(num_cpu) + 1:
                            command = "%ssamtools mpileup --max-depth 10000 -AQ 0 -l %stemp%s.bed %s | awk '$5 ~/[\+,\-]/ {print $1\"\t\"$2-1\"\t\"$2}' > %smpileup_positions%s.bed" % (
                                path_samtools, out, str(counter), bam_file, out,
                                str(counter))

                            proc_mpileup = subprocess.Popen(command, shell=True)
                            ps.append(proc_mpileup)
                            counter += 1

                        for proc_mpileup in ps:
                            proc_mpileup.wait()

                        # 10.3.2 Strelka indel calling on selected positions (from
                        # previous step 10.3.1)
                        counter = 1
                        os.system("touch %smpileup_positions.bed" % (out))
                        while counter < int(num_cpu) + 1:
                            os.system(
                                "cat %smpileup_positions%s.bed >> %smpileup_positions.bed"
                                % (out, str(counter), out))

                            counter += 1

                        os.system(
                            "%sbedtools sort -i %smpileup_positions.bed > %stemp_sorted.bed ; mv  %stemp_sorted.bed %smpileup_positions.bed"
                            % (path_bedtools, out, out, out, out))

                        print(
                            "\nCalling indels with Strelka on indel positions identified by samtools mpileup...\n"
                        )

                        os.system(
                            "%svcftools  --gzvcf %s%s_sorted.vcf.gz --minGQ 30 --minDP 2 --bed %smpileup_positions.bed --recode --recode-INFO-all --out %sindels_only"
                            % (path_vcftools, out, sample_name, out, out))
                        os.system("touch  %slogs/indels_only.log" % (out))

                        print(
                            "\nCompleted indel calling with Strelka.\n"
                        )


                        # 10.3.3 If intensive mode variant calling was performed, DNAscan called snvs
                        # with Freebayes and indels with Strelka, resulting in two vcf files.
                        # For the annotation step these two files are merged together.

                        os.system(
                            "%svcftools  --gzvcf %s%s_sorted.vcf.gz --minGQ 30 --minDP 2 --exclude-bed %smpileup_positions.bed  --recode --recode-INFO-all --out %sSNPs_only"
                            % (path_vcftools, out, sample_name, out, out))
                        os.system("touch %slogs/SNPs_only.log" % (out))

                        print("\nCompleted SNV calling with Strelka.\n")

                        print("\nMerging the SNV and indel VCFs to enable downstream annotation to take place...\n")
                        os.system(
                            "bgzip  %sSNPs_only.recode.vcf ; bgzip %sindels_only.recode.vcf "
                            % (out, out))
                        os.system(
                            "%stabix -p vcf %sSNPs_only.recode.vcf.gz ; %stabix -p vcf %sindels_only.recode.vcf.gz "
                            % (path_tabix, out, path_tabix, out))
                        os.system(
                            "%sbcftools merge --merge all --force-samples %sSNPs_only.recode.vcf.gz %sindels_only.recode.vcf.gz -O v -o %s%s.vcf "
                            % (path_bcftools, out, out, out,
                            sample_name))
                        os.system(
                            "perl %svcf-sort.pl  %s%s.vcf | bgzip -c > %s%s_sorted.vcf.gz"
                            % (path_scripts, out, sample_name, out, sample_name))

                        variant_results_file = "%s%s_sorted.vcf.gz" % (out,
                                                                        sample_name)

                        is_variant_file_OK(variant_results_file, "Vcf", "variantcalling")

                        os.system("%stabix -p vcf %s%s_sorted.vcf.gz" % (path_tabix, out, sample_name))

                        if not debug:
                            os.system(
                                "rm %s%s.vcf* %slogs/SNPs_only.log %smpileup_positions* %slogs/indels_only.log"
                                % (out, sample_name, out, out, out))
                            
                        print("\nSuccessfully merged SNV and indel calls.\n")
                       
                    if not debug:
                        os.system("rm -r %sstrelka %s/results/%s_strelka.vcf.gz.tbi" % (out, out, sample_name))
                        
                else:
                    print('\nSmall variant calling with Strelka requires alignment from paired-end reads.\n')

            os.system("touch %slogs/VC.log" % (out))

            print("\nCompleted SNV and indel calling.\n")

# 11. Perform variant hard filtering

if filter_string and len(variant_results_file) != 0:
    print("\nHard filtering of SNV/indel variants is being performed...\n")
    # filters for those variants which have passed all Strelka filters plus custom string (according to strelka vcf file parameters)
    print('%sbcftools filter -i \" %s \" %s | bgzip -c > %s%s_sorted_filtered.vcf.gz ; %stabix -fp vcf %s%s_sorted_filtered.vcf.gz'
    % (path_bcftools, filter_string, variant_results_file, out, sample_name, path_tabix, out, sample_name))
    os.system(
    "%sbcftools filter -i \'%s\' %s | bgzip -c > %s%s_sorted_filtered.vcf.gz ; %stabix -fp vcf %s%s_sorted_filtered.vcf.gz " % (
     path_bcftools, filter_string, variant_results_file, out, sample_name, path_tabix, out, sample_name))

    variant_results_file = "%s%s_sorted_filtered.vcf.gz" % (out, sample_name)
    is_variant_file_OK(variant_results_file, "Vcf", "variantcalling")

    print("\nVariant hard filtering is complete.\n")

# 12. Performs known expansions search with ExpansionHunter
# Known expansions have to be described in a json file placed in the json
# folder (paths_configs.py).

if expansion:
    if "EH.log" in os.listdir(out):

        print(
            "WARNING: The presence of EH.log in logs is telling you that the expansion scan was already peformed, please remove SV.log if you wish to perform this stage anyway\n"
        )

    else:

        print("\nExpansionHunter is scanning the genome for repeat expansions specified in the variant_catalog json...\n")

        os.system(
            "%sExpansionHunter --reads %s --reference %s  --variant-catalog %s/variant_catalog_%s.json --output-prefix %s/temp_EH"
            % (path_expansionHunter, bam_file, path_reference,
               path_expansionHunter_jsons, reference, out))
        os.system(
            "mv %s/temp_EH.vcf %s/results/%s_expansions.vcf ; bgzip %s/results/%s_expansions.vcf ; %stabix -p vcf %s/results/%s_expansions.vcf.gz"
            % (out, out, sample_name, out, sample_name, path_tabix, out,
               sample_name))

        expansion_results_file = "%s/results/%s_expansions.vcf.gz" % (out, sample_name)

        is_variant_file_OK(expansion_results_file, "Vcf", "expansion")

        os.system("touch  %slogs/EH.log" % (out))

        print("\nRepeat expansion scanning is complete.\n")

        if not debug:
            os.system("rm %stemp_EH.json*  %stemp_EH.log" % (out, out))

# 13. Structural Variant calling

if SV:

    if "SV.log" in os.listdir(out + "logs"):

        print(
            "WARNING: The presence of SV.log in logs is telling you that structural variant calling was already peformed, please remove SV.log if you wish to perform this stage anyway\n"
        )

    else:
        if paired == "1":
            if BED == "True":
                print(
                    "\nStructural variants are being called with Manta...\n"
                )
                os.system("bgzip -c %s  > %s/temp.bed.gz" % (path_bed, out))
                os.system(
                    "%ssortBed -i %s/temp.bed.gz | bgzip -c > %s/sorted.bed.gz" %
                    (path_bedtools, out, out))

                os.system("%stabix -p bed %s/sorted.bed.gz" % (path_tabix, out))
                os.system("mkdir %smanta" % (out))
                os.system(
                    "%sconfigManta.py --bam %s --referenceFasta %s --runDir %smanta --callRegions %s/sorted.bed.gz"
                    % (path_manta, bam_file, path_reference, out, out))

            else:
                print(
                    "\nStructural variants are being called with Manta...\n")

                os.system("mkdir %smanta" % (out))
                os.system(
                    "%sconfigManta.py --bam %s --referenceFasta %s --runDir %smanta"
                    % (path_manta, bam_file, path_reference, out))

            os.system("%smanta/runWorkflow.py -j %s -m local" % (out, num_cpu))
            os.system("gzip -d %s/manta/results/variants/diploidSV.vcf.gz")
            os.system("%s/convertInversion.py %s/samtools %s %s/manta/results/variants/diploidSV.vcf | bgzip -c > %s/results/%s_manta_SV.vcf.gz" % (
            path_scripts, path_samtools, path_reference, out, out, sample_name))
            os.system(
                "tabix -p vcf %s/results/%s_manta_SV.vcf.gz"
                % (out, sample_name))

            if mode == "fast":
                #13.1 Manta is used to call all SVs in fast mode
                os.system("tabix -p vcf %s/results/%s_manta_SV.vcf.gz" % (out, sample_name))

                SV_results_file = "%s/results/%s_manta_SV.vcf.gz" % (out, sample_name)

            else:
                os.system ("gzip -d %s/results/%s_manta_SV.vcf.gz" % (out, sample_name))

                manta_SV_results_file = "%s/results/%s_manta_SV.vcf" % (out, sample_name)

            if not debug:
                os.system("rm -r %stemp.bed.gz  %ssorted.bed.gz %smanta" %
                          (out, out, out))

            is_variant_file_OK(manta_SV_results_file, "Vcf", "SV")

            print("\nStructural variant calling with Manta is complete.\n")

            if mode == "normal" or mode == "intensive":
                # 13.2 In normal mode, Delly is used to call inversions and deletions and Manta calls all SV types.
                # In intensive mode, both Delly and Manta call all SV events

                os.system("mkdir %sdelly" % (out))

                if mode == "normal":

                    print("\nDelly will additionally call inversion and deletion structural variants...\n")

                    os.system("%sdelly call -t DEL -g %s -o %sdelly/%s_delly_DEL.bcf -x %s %s" % (
                    path_delly, path_reference, out, sample_name, path_delly_exclude_regions, bam_file))
                    os.system("%sdelly call -t INV -g %s -o %sdelly/%s_delly_INV.bcf -x %s %s" % (
                    path_delly, path_reference, out, sample_name, path_delly_exclude_regions, bam_file))
                    os.system("%sbcftools view %sdelly/%s_delly_DEL.bcf > %sdelly/%s_delly_DEL.vcf" % (
                    path_bcftools, out, sample_name, out, sample_name))
                    os.system("%sbcftools view %sdelly/%s_delly_INV.bcf > %sdelly/%s_delly_INV.vcf" % (
                    path_bcftools, out, sample_name, out, sample_name))
                    os.system("grep '^#' %sdelly/%s_delly_DEL.vcf > %sdelly/%s_delly.vcf ; grep -v '^#' %sdelly/%s_delly_DEL.vcf %sdelly/%s_delly_INV.vcf >> %sdelly/%s_delly.vcf" % (
                    out, sample_name, out, sample_name, out, sample_name, out, sample_name, out, sample_name))

                else:
                    print("\nStructural variants are being called with Delly...\n")

                    os.system("%sdelly call -g %s -o %sdelly/%s_delly.bcf -x %s %s" % (
                    path_delly, path_reference, out, sample_name, path_delly_exclude_regions, bam_file))
                    os.system("%sbcftools view %sdelly/%s_delly.bcf > %sdelly/%s_delly.vcf" % (
                    path_bcftools, out, sample_name, out, sample_name))

                os.system("mv %sdelly/%s_delly_SV.vcf %s/results/" % (out, sample_name, out))

                delly_SV_results_file = "%s/results/%s_delly_SV.vcf" % (out, sample_name)

                is_variant_file_OK(delly_SV_results_file, "Vcf", "SV")

                print("\nStructural variant calling with Delly is complete.\n")

                # 13.3 In normal and intensive modes, SV calls with Manta and Delly are merged together using SURVIVOR to create a union set of structural variants.

                print("\nStructural variants called with Manta and Delly are being merged with SURVIVOR to create a union callset...\n")

                os.system("ls %s/results/*SV.vcf > %s/results/survivor_sample_files" % (out, out))
                os.system("%sSURVIVOR merge %s/results/survivor_sample_files 1000 1 1 1 0 30 > %s/results/%s_SV_merged.vcf" % (
                path_SURVIVOR, out, out, sample_name))
                os.system("bgzip %s/results/%s_SV_merged.vcf ; tabix -p vcf %s/results/%s_SV_merged.vcf.gz" % (out, sample_name, out_sample_name))

                SV_results_file = "%s/results/%s_SV_merged.vcf.gz"

                is_variant_file_OK(SV_results_file, "Vcf", "SV")

                print("\nMerging of structural variant calls is complete.\n")

                if not debug:
                    os.system("rm -r %sdelly %s/results/survivor_sample_files" % (out, out))

            os.system("touch %slogs/SV.log" % (out))

        else:
            print(
                "\n\nWARNING: Structural variant calling cannot be performed with single end reads. If you want this to be performed, please provide paired-end aligned reads to DNAscan.\n\n"
            )

# 14. Transposable element insertion detection with MELT
if MEI:
    if "mei.log" in os.listdir(out + "logs"):

        print(
            "WARNING: The presence of MEI.log in logs is telling you that transposable element calling was already peformed, please remove MEI.log if you wish to perform this stage anyway\n"
        )
    else:

        print("\nTransposable element insertions are being scanned for with MELT...\n")

        os.system("mkdir %smelt" % (out))

        if reference == "hg19":
            os.system("ls %sme_refs/1KGP_Hg19/*zip | sed 's/\*//g' > %smelt/transposon.list" % (path_melt, out))

            melt_bed = "%sadd_bed_files/1KGP_Hg19/hg19.genes.bed" % (path_melt)

        else:
            os.system("ls %sme_refs/Hg38/*zip | sed 's/\*//g' > %smelt/transposon.list" % (path_melt, out))

            melt_bed = "%sadd_bed_files/Hg38/Hg38.genes.bed" % (path_melt)

        os.system("%sjava -Xmx%sg -jar %sMELT.jar Single -bamfile %s -h %s -t %smelt/transposon.list -n %s -w %smelt" % (path_java, RAM_GB, path_melt, bam_file, path_reference, melt_bed, out))

        os.system("cat %melt/SVA.final_comp.vcf | grep '^#' > %smelt/%s.header.txt" % (out, out, sample_name))
        os.system("cat %melt/SVA.final_comp.vcf | grep -v '^#' > %smelt/%s.sva.vcf" % (out, out, sample_name))
        os.system("cat %melt/LINE1.final_comp.vcf | grep -v '^#' > %smelt/%s.line1.vcf" % (out, out, sample_name))
        os.system("cat %melt/ALU.final_comp.vcf | grep -v '^#' > %smelt/%s.alu.vcf" % (out, out, sample_name))
        os.system("cat %smelt/%s.header.txt %smelt/%s.sva.vcf %smelt.%s.line1.vcf %smelt/%s.alu.vcf | perl %svcf-sort.pl -c | bgzip -c > %s/results/%s_MEI.vcf.gz" % (out, sample_name, out, sample_name, out, sample_name, out, sample_name, path_scripts, out, sample_name))
        os.system("tabix -p vcf %s/results/%s_MEI.vcf.gz" % (out, sample_name))

        MEI_results_file = "%s/results/%s_MEI.vcf.gz" % (out, sample_name)

        is_variant_file_OK(MEI_results_file, "Vcf", "MEI")

        if not debug:
            os.system("rm -r %smelt" % (out))

        os.system("touch  %slogs/mei.log" % (out))

        print("\nTransposable element insertion scanning with MELT is complete.\n")

# 15. Annotation with Annovar with optional missense variant prioritisation according to ACMG guidelines (intervar_20180118 database)

if annotation:
    if "annovar.log" in os.listdir(out + "logs"):

        print(
            "WARNING: The presence of annovar.log in logs is telling you that annotation was already peformed, please remove annovar.log if you wish to perform this stage anyway\n"
        )

        variant_results_file = "%sresults/%s_annotated.vcf.gz" % (out,
                                                                  sample_name)
    else:

        print("\nAnnotation is being performed using Annovar, with the databases and respective operations defined in paths_configs.py...\n")

        os.system(
            "perl %stable_annovar.pl  --thread %s --vcfinput %s %s -buildver %s -remove -protocol %s -operation %s -nastring . --outfile %s/annovar.vcf"
            % (path_annovar, num_cpu, variant_results_file, path_annovar_db,
               reference, annovar_protocols, annovar_operations, out))

        if not debug and not alsgenescanner:
            os.system(
                "rm %sannovar.vcf.%s_multianno.txt %sannovar.vcf.avinput" %
                (out, reference, out))

        os.system(
            "mv %s/annovar.vcf.%s_multianno.vcf %sresults/%s_annotated.vcf ; bgzip -f %sresults/%s_annotated.vcf ; %stabix -fp vcf %sresults/%s_annotated.vcf.gz"
            % (out, reference, out, sample_name, out, sample_name, path_tabix, out,
               sample_name))

        os.system("mv %s %sresults/" % (variant_results_file, out))
        os.system("mv %s.tbi %sresults/" % (variant_results_file, out))

        variant_results_file = "%sresults/%s_annotated.vcf.gz" % (out,
                                                                  sample_name)
        is_variant_file_OK(variant_results_file, "Vcf", "annotation")

        os.system("touch  %slogs/annovar.log" % (out))

        print("\nAnnotation with ANNOVAR is complete.\n")

        if SV:
        #15.1 Structural variant annotation and prioritisation is carried out using AnnotSV
            print("\nStructural variant annotation is being performed with AnnotSV...\n")

            os.system("export ANNOTSV=%s" % (path_annotsv))

            if alsgenescanner or len(path_gene_list) != 0:
                candidate_gene_cmd = "-candidateGenesFile %s" % (path_gene_list)

            if reference == "hg19":
                genome_build = "GRCh37"

            else:
                genome_build = "GRCh38"

            os.system("AnnotSV -annotationsDir %s/share/AnnotSV/ -bcftools %sbcftools -bedtools %sbedtools -SvinputFile %s %s -genomeBuild %s -outputFile %s/results/%s_annotated_SV -SVminSize 30 %s" % (
            path_annotsv, path_bcftools, path_bedtools, SV_results_file, candidate_gene_cmd, genome_build, out, sample_name, annotsv_custom_options))

            SV_annotation_file = "%s/results/%s_annotated_SV.tsv" % (out, sample_name)

            if not debug:
                os.system("rm %s/results/*sorted.bed %s/results/*breakpoints.bed" % (out, out))

            print("\nStructural variant annotation is complete.\n")

        print("\nAnnotation and prioritisation is complete.\n")

else:
    if variant_results_file:
        os.system("mv %s* %sresults/" % (variant_results_file, out))

        variant_results_file = "%s%s_sorted.vcf.gz" % (out, sample_name)


# 16. Microbes screening

if virus or bacteria or custom_microbes:

    if "microbes.log" in os.listdir(out + "logs"):

        print(
            "WARNING: The presence of microbes.log in logs is telling you that microbes scanning was already peformed, please remove microbes.log if you wish to perform this stage anyway\n"
        )
    else:
        # 16.1 Exctract non human reads
        print("\nExtracting non-human reads...\n")
        os.system(
            "%ssamtools view -@ %s -hf 4 %s | %ssamtools bam2fq -s %ssingleton_reads.fastq -@ %s - > %sunaligned_reads.fastq ; cat %ssingleton_reads.fastq >> %sunaligned_reads.fastq ; gzip  %sunaligned_reads.fastq "
            % (path_samtools, num_cpu, bam_file, path_samtools, out, num_cpu,
               out, out, out, out))

        # 16.2 Identifies present non-human microbes
        if virus:

            # 16.2.1 Identifies present viruses
            print("\nIdentifying viruses present in sample...\n")
            os.system(
                "%shisat2 --no-spliced-alignment -p %s -x %s -U %sunaligned_reads.fastq.gz | %ssamtools view -@ %s -hSb -  | %ssamtools sort -@ %s -T %stemp.file -o %soutput_virus.bam -"
                % (path_hisat, num_cpu, path_virus_index, out, path_samtools,
                   num_cpu, path_samtools, num_cpu, out, out))
            os.system(
                "%ssamtools index -@ %s %soutput_virus.bam; %ssamtools idxstats %soutput_virus.bam > %svirus_stats.txt"
                % (path_samtools, num_cpu, out, path_samtools, out, out))

            # 16.2.2 Generates virus report
            print("\nGenerating virus coverage and stats report...\n")
            os.system(
                "awk \'{print $1}\' %svirus_stats.txt > %svirus_list.txt ; for i in $(cat %svirus_list.txt | grep C_); do printf \"$i \"; %ssamtools depth %soutput_virus.bam -r $i | awk \'$3>0 {print $0}\' | wc -l; done > %svirus_coverage_stats.txt"
                % (out, out, out, path_samtools, out, out))

            virus_coverage_stats = open('%svirus_coverage_stats.txt' % (out),
                                        'r')

            virus_coverage_stats_lines = virus_coverage_stats.readlines()
            virus_stats_file = open('%svirus_stats.txt' % (out), 'r')
            virus_stats_file_lines = virus_stats_file.readlines()

            i = 0

            virus_results = open('%sresults/virus_results.txt' % (out), 'w')
            virus_results.write("Id\tGenome_length\tNumber_of_reads\tCoverage\n")

            while i < len(virus_coverage_stats_lines):
                virus_results.write(
                    "%s\t%s\t%s\t%s\n" %
                    (virus_stats_file_lines[i].split('\t')[0],
                     virus_stats_file_lines[i].split('\t')[1],
                     virus_stats_file_lines[i].split('\t')[2],
                     virus_coverage_stats_lines[i].split(' ')[1].strip()))

                i += 1

            virus_results.close()

            if not debug:
                os.system("rm %svirus_stats.txt  %svirus_coverage_stats.txt" %
                          (out, out))

            print("\nVirus coverage and stats report is now available.\n"
                 )

        if bacteria:

            # 16.2.3 Identifies present bacteria
            print("\nIdentifying bacteria present in sample...\n")

            os.system(
                "%shisat2 --no-spliced-alignment -p %s -x %s -U %sunaligned_reads.fastq.gz | %ssamtools view -hSb -  | %ssamtools sort -T %stemp.file -o %soutput_bacteria.bam -"
                % (path_hisat, num_cpu, path_bacteria_index, out,
                   path_samtools, path_samtools, out, out))
            os.system(
                "%ssamtools index -@ %s %soutput_bacteria.bam; %ssamtools idxstats %soutput_bacteria.bam > %sbacteria_stats.txt"
                % (path_samtools, num_cpu, out, path_samtools, out, out))

            # 16.2.4 Generates bacteria report
            print("\nGenerating bacteria coverage and stats report...\n")

            os.system(
                "awk \'{print $1}\' %sbacteria_stats.txt > %sbacteria_list.txt ; for i in $(cat %sbacteria_list.txt| grep ref); do printf \"$i \"; %ssamtools depth %soutput_bacteria.bam -r $i | awk \'$3>0 {print $0}\' | wc -l; done > %sbacteria_coverage_stats.txt"
                % (out, out, out, path_samtools, out, out))

            bacteria_coverage_stats = open(
                '%sbacteria_coverage_stats.txt' % (out), 'r')

            bacteria_coverage_stats_lines = bacteria_coverage_stats.readlines()
            bacteria_stats_file = open('%sbacteria_stats.txt' % (out), 'r')
            bacteria_stats_file_lines = bacteria_stats_file.readlines()

            i = 0

            bacteria_results = open('%sresults/bacteria_results.txt' % (out),
                                    'w')
            bacteria_results.write(
                "Id\tGenome_length\tNumber_of_reads\tCoverage\n")

            while i < len(bacteria_coverage_stats_lines):
                bacteria_results.write(
                    "%s\t%s\t%s\t%s\n" %
                    (bacteria_stats_file_lines[i].split('|')[1],
                     bacteria_stats_file_lines[i].split('\t')[1],
                     bacteria_stats_file_lines[i].split('\t')[2],
                     bacteria_coverage_stats_lines[i].split(' ')[1].strip()))

                i += 1

            bacteria_results.close()

            if not debug:
                os.system(
                    "rm %sbacteria_stats.txt  %sbacteria_coverage_stats.txt" %
                    (out, out))

            print("\nBacteria coverage and stats report is now available.\n")

        if custom_microbes:

            # 16.2.5 Identifies present user-selected microbes
            print("\nIdentifying user-selected microbes that are present in sample...\n")

            os.system(
                "%shisat2 --no-spliced-alignment -p %s -x %s -U %sunaligned_reads.fastq.gz | %ssamtools view -hSb -  | %ssamtools sort -T %stemp.file -o %soutput_custom_microbes.bam -"
                % (path_hisat, num_cpu, path_custom_microbes_index, out,
                   path_samtools, path_samtools, out, out))
            os.system(
                "%ssamtools index -@ %s %soutput_custom_microbes.bam; %ssamtools idxstats %soutput_custom_microbes.bam > %scustom_microbes_stats.txt"
                % (path_samtools, num_cpu, out, path_samtools, out, out))

            # 16.2.6 Generates user-selected microbes report
            print("\nGenerating user-selected microbe coverage and stats report...\n")

            os.system(
                "awk \'{print $1}\' %scustom_microbes_stats.txt > %scustom_microbes_list.txt ; for i in $(cat %scustom_microbes_list.txt| grep ref); do printf \"$i \"; %ssamtools depth %soutput_custom_microbes.bam -r $i | awk \'$3>0 {print $0}\' | wc -l; done > %scustom_microbes_coverage_stats.txt"
                % (out, out, out, path_samtools, out, out))

            custom_microbes_coverage_stats = open(
                '%scustom_microbes_coverage_stats.txt' % (out), 'r')

            custom_microbes_coverage_stats_lines = custom_microbes_coverage_stats.readlines(
            )
            custom_microbes_stats_file = open(
                '%scustom_microbes_stats.txt' % (out), 'r')
            custom_microbes_stats_file_lines = custom_microbes_stats_file.readlines(
            )

            i = 0

            custom_microbes_results = open(
                '%sresults/custom_microbes_results.txt' % (out), 'w')
            custom_microbes_results.write(
                "Id\tGenome_length\tNumber_of_reads\tCoverage\n")

            while i < len(custom_microbes_coverage_stats_lines):
                custom_microbes_results.write(
                    "%s\t%s\t%s\t%s\n" %
                    (custom_microbes_stats_file_lines[i],
                     custom_microbes_stats_file_lines[i].split('\t')[1],
                     custom_microbes_stats_file_lines[i].split('\t')[2],
                     custom_microbes_coverage_stats_lines[i].split(
                         ' ')[1].strip()))

                i += 1

            custom_microbes_results.close()

            if not debug:
                os.system(
                    "rm %scustom_microbes_stats.txt  %scustom_microbes_coverage_stats.txt"
                    % (out, out))

            print("\nMicrobe coverage and stats report is now available.\n"
                 )

        os.system("touch  %slogs/microbes.log" % (out))

        print("\nMicrobe screening is complete.\n")

# 17. Alignment report generation ( samtools flagstat and stats )

if alignment_report:
    if "alignment_report.log" in os.listdir(out + "logs"):

        print(
            "WARNING: The presence of alignment_report.log in logs is telling you that the alignment report was already produced, please remove alignment_report.log if you wish to perform this stage anyway\n"
        )

    else:

        print("\nGenerating alignment report...\n")

        os.system("%ssamtools flagstat -@ %s %s > %sreports/%s_flagstat.txt" %
                  (path_samtools, num_cpu, bam_file, out, sample_name))
        os.system("%ssamtools stats -@ %s  %s > %sreports/%s_stats.txt" %
                  (path_samtools, num_cpu, bam_file, out, sample_name))
        os.system("touch  %slogs/alignment_report.log" % (out))

        print("\nAlignment report is now available.\n")


# 18. Sequencing data report generation ( fastqc )

if sequencing_report:
    if "sequencing_report.log" in os.listdir(out + "logs"):

        print(
            "WARNING: The presence of sequencing_report.log in logs is telling you that the sequence report was already produced, please remove sequencing_report.log if you wish to perform this stage anyway\n"
        )

    else:
        if format == "fastq":
            if path_java != "":
                java_option = "-j " + path_java + "java"
            else:
                java_option = ""

            print("\nGenerating sequencing report...\n")

            os.system("%sfastqc %s -o %sreports -f %s -t %s %s %s" %
                  (path_fastqc, java_option, out, format, num_cpu, input_file,
                   input_file2))
            os.system("touch  %slogs/sequencing_report.log" % (out))

            print("\nSequencing report is now available.\n")

        else:
            print(
                "WARNING: The sequencing report cannot be generated as the input data is not in fastq format.\n"
            )


# 19. Snv and indel calling report generation ( bcftools stats )

if calls_report:
    if "calls_report.log" in os.listdir(out + "logs"):

        print(
            "WARNING: The presence of calls_report.log in logs is telling you that the calls report was already produced, please remove calls_report.log if you wish to perform this stage anyway\n"
        )

    else:
        print("\nGenerating SNV and indel calls report...\n")

        os.system(
            "%sbcftools stats --threads %s %s > %sreports/%s_vcfstats.txt" %
            (path_bcftools, num_cpu, variant_results_file, out, sample_name))
        os.system("touch  %slogs/calls_report.log" % (out))

        print("\nCalls report for SNVs and indels are now available.\n")

# 20. Html report generation ( Multiqc )

if alignment_report or calls_report or sequencing_report:
    if "multiqc.log" in os.listdir(out + "logs"):

        print(
            "WARNING: The presence of multiqc.log in logs is telling you that the multiqc report was already produced, please remove multiqc.log if you wish to perform this stage anyway\n"
        )

    else:

        print("\nGenerating HTML report...\n")

        os.system(
            "%smultiqc -o %sreports %sreports" % (path_multiqc, out, out))
        os.system("touch  %slogs/multiqc.log" % (out))

        print("\nHTML report created.\n")

# 21. Annotated variants report generation

if results_report:
    if "annovar.log" not in os.listdir(out + "logs") or not path_gene_list:

        print(
            "WARNING: Either the annotation was not peformed or path_gene_list was not provided in paths_configs.py, please perform annotation using the -annotation flag and/or specify the a gene list in paths_configs.py if you wish to generate an annotation results report.\n"
        )

    else:
        if "results_report.log" in os.listdir(out + "logs"):
            print(
                "WARNING: The presence of results_report.log in logs is telling you that the results report was already produced, please remove results_report.log if you wish to perform this stage anyway\n"
            )

        else:
            print("\nGenerating report of annotated variant calls...\n")

            os.system("zcat %s > %stemp.vcf" % (variant_results_file, out))

            #vcf = open('%stemp.vcf' % (out), 'r')

            #vcf_lines = vcf.readlines()

            gene_list_file = open(path_gene_list)

            gene_list_lines = gene_list_file.readlines()

            gene_list = gene_list_lines

            out_file_all = open(
                '%sreports/%s_all_variants.txt' % (out, sample_name), 'w')

            counter = 0

            for i in gene_list:

                with open('%stemp.vcf' % (out)) as vcf:

                    for j in vcf:

                        check1 = re.search(
                            r'(^chr)|(^[0-9,X,Y,M]+\t)', j, flags=0)

                        check = re.search(
                            "=%s;" % (i.strip().upper()), j, flags=0)

                        if check and check1:

                            infos = j.split('ANNOVAR_DATE')[1][12:].split(
                                'ALLELE_END')[0].replace(";", "\t")

                            if counter == 0:

                                replaced_1 = re.sub(
                                    '=[a-z,A-Z,0-9,\.,\_,\-,:,>,<]+', '',
                                    infos)

                                out_file_all.write(
                                    'CHR\tPosition\tRef\tAlt\tGenotype\t%s\n' %
                                    (replaced_1))

                                counter = 1

                            replaced = re.sub('[a-z,A-Z,0-9,\.,\_,\-,:,>,<]+=',
                                              '', infos)

                            out_file_all.write(
                                '%s\t%s\t%s\t%s\t%s\t%s\n' %
                                (j.split('\t')[0], j.split('\t')[1],
                                 j.split('\t')[3], j.split('\t')[4],
                                 j.split('\t')[-1].split(':')[0], replaced))

            out_file_all.close()

        #21.1 knotAnnotSV SV report generation
            if SV:
                print("\nGenerating SV annotation HTML report...\n")

                os.system("mkdir %s%s_SVanno" % (out, sample_name))

                os.system("perl %sknotAnnotSV.pl --configFile %s/config_AnnotSV.yaml --annotSVfile %s --outDir %s%s_SVanno --genomeBuild %s " % (path_knotannotsv, path_knotannotsv, SV_annotation_file, out, reference))

                os.system("mv %s%s_SVanno/%s_SVanno.annotated.html %s/reports/%s_SVannotatedvariants.html" % (out, sample_name, sample_name, out))

                if not debug:
                    os.system("rm -r %s%s_SVanno" % (out, sample_name))

                print("\nSV HTML report created.\n")

            os.system("touch  %slogs/results_report.log" % (out))

            print("\nCalls report for annotated variants is now available.\n")


# 22. Starting iobio services

if iobio:
    if "iobio.log" in os.listdir(out + "logs"):

        print(
            "WARNING: The presence of iobio.log in logs is telling you that the iobio services were already started, please remove iobio.log if you wish to start them again\n"
        )

    else:
        os.system("touch  %slogs/iobio.log" % (out))

        print(
            "\n\nIobio services can now be started. Select the service by copying and pasting the relevant iobio URLs into your browser (https://vcf.iobio.io, https://bam.iobio.io and/or https://gene.iobio.io) and upload your data into the selected service\n\nIf you want to explore your variant calling results please copy and paste the following URL into your browser and upload the vcf file (../%s%s_sorted.vcf.gz):\n\n"
            % (out, sample_name),
            end='',
            flush=True)

        if "annovar.log" in os.listdir(out + "logs"):
            print(
                "https://gene.iobio.io/?species=Human&rel0=proband&rel1=mother&rel2=father&genes=",
                end='',
                flush=True)

            a = {}

            os.system("zcat %s > %stemp.vcf" % (variant_results_file, out))

            file = open('%stemp.vcf' % (out), 'r')
            file_lines = file.readlines()

            for i in file_lines:
                check = re.search(r'(^chr)|(^[0-9,X,Y,M]+\t)', i, flags=0)

                if check:
                    a[i.split('Gene.refGene=')[1].split(';')[0]] = []

            for i in a.keys():
                print('%s,' % (i.split(',')[0]), end='', flush=True)
            print(' \n\n')

        else:
            print("https://gene.iobio.io")

if alsgenescanner: 
    print("\n\nALSGeneScanner is running...\n\n")
    os.system(
        "python3 %s/alsgenescanner.py %s/annovar.vcf.%s_multianno.txt %s/results/%s_alsgenescanner_all.txt"
        % (path_scripts, out, reference, out, sample_name))
    os.system(
        "cat %s/results/%s_alsgenescanner_all.txt | head -1 > %s/results/%s_alsgenescanner_alsod.txt; cat %s/results/%s_alsgenescanner_all.txt | grep -iwf %s/list_genes_alsod.txt >> %s/results/%s_alsgenescanner_alsod.txt"
        % (out, sample_name, out, sample_name, out, sample_name, path_to_db,
           out, sample_name))
    os.system(
        "cat %s/results/%s_alsgenescanner_all.txt | head -1 > %s/results/%s_alsgenescanner_clinvar.txt ; cat %s/results/%s_alsgenescanner_all.txt | grep -iwf %s/list_genes_clinvar.txt >> %s/results/%s_alsgenescanner_clinvar.txt"
        % (out, sample_name, out, sample_name, out, sample_name, path_to_db,
           out, sample_name))
    os.system(
        "cat %s/results/%s_alsgenescanner_all.txt | head -1 > %s/results/%s_alsgenescanner_manual_review.txt ; cat %s/results/%s_alsgenescanner_all.txt | grep -iwf %s/list_genes_manual_review.txt >> %s/results/%s_alsgenescanner_manual_review.txt"
        % (out, sample_name, out, sample_name, out, sample_name, path_to_db,
           out, sample_name))
    os.system(
        "cat %s/results/%s_alsgenescanner_all.txt | head -1 > %s/results/%s_alsgenescanner_all_ranked.txt ; cat %s/results/%s_alsgenescanner_all.txt | grep -i pathog | sed 's/ /_/g'| sort -k10nr >> %s/results/%s_alsgenescanner_all_ranked.txt ; cat %s/results/%s_alsgenescanner_all.txt | grep '^chr'  | grep -iv pathog | sed 's/ /_/g'| sort -k10nr >> %s/results/%s_alsgenescanner_all_ranked.txt  "
        % (out, sample_name, out, sample_name, out, sample_name, out, sample_name,
           out, sample_name, out, sample_name))
    print("\n\nALSGeneScanner is complete.\n\n")
