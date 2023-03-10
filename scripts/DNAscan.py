#!/usr/bin/env python3

################################################################
# Program: DNAscan
# Version 2.0
# Author: Heather Marriott (heather.marriott@kcl.ac.uk) and Alfredo Iacoangeli (alfredo.iacoangeli@kcl.ac.uk)
#################################################################

#################################################################
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
#       8.1.1 Alignment for SNP/indel use only
#       8.1.2 Alignment (general)
#   8.2 Aligns single end reads
#       8.2.1 Alignment for SNP/indel use only
#       8.2.2 Alignment (general)
# 9. If input file is a sam file, it converts it into bam
# 10. Variant (snv and indel) calling with Strelka2
# 11. Perform variant hard filtering
# 12. Perform known expansions search with ExpansionHunter
# 13. Compute genome-wide short tandem repeat profiles with ExpansionHunter Denovo
#   13.1 Convert the novel/non-reference loci identified with ExpansionHunter Denovo to variant catalog format
#   13.2 Genotype the novel/non reference variant catalog with ExpansionHunter
# 14. Structural Variant calling
#   14.1 Manta (all SV)
#   14.2 Delly (general SV)
#   14.3 SV calls are merged to create union callset
# 15. Transposable element insertion detection with MELT
# 16. Annotation with Annovar -  with optional missense variant prioritisation according to ACMG guidelines (intervar_20180118 database)
#   16.1 AnnotSV for structural variant/transposable element annotation and prioritisation
# 17. Microbes screening
#   17.1 Extract non human reads
#   17.2 Identifies present non-human microbes
#       17.2.1 Identifies present viruses
#       17.2.2 Generates virus report
#       17.2.3 Identifies present bacteria
#       17.2.4 Generates bacteria report
#       17.2.5 Identifies present user-selected microbes
#       17.2.6 Generates user-selected microbes report
# 18. Alignment report generation ( samtools flagstat and stats )
# 19. Sequencing data report generation ( fastqc )
# 20. Snv and indel calling report generation ( bcftools stats )
# 21. Html report generation ( Multiqc )
# 22. Annotated variants results report generation
#   22.1 knotAnnotSV annotated structural variants report generation
#   22.2 Concise results report for all variants called (SV, MEI, expansions, SNVs and indels)
# 23. Starting iobio services
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
import create_variant_catalog

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
alsgenescanner = args.alsgenescanner
exome = args.exome
format = args.format
paired = args.paired
reference = args.reference
input_file = args.input_file
input_file2 = args.input_file2
expansion = args.expansion
STR = args.STR
genotypeSTR = args.genotypeSTR
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
fast_mode = args.fast_mode

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
        sys.exit("WARNING: %s does not exist - DNAscan will now terminate.\n" % file)

# 6. Bed splitting: splitting the analysis region into subsets of equal length.
#

if alsgenescanner:

    alignment = True
    variantcalling = True
    annotation = True
    SV = True
    expansion = True
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
    # reference genome index or the default exome bed.
    if exome == "True":
        os.system("zcat %sdb/exome_%s.bed.gz > %stmp/exome_%s.bed" %
                  (dnascan_dir, reference, out, reference))
        path_bed = "%stmp/exome_%s.bed" % (out, reference)

        if len(path_bed) != 0:
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
        "\nDuplicates will be removed.\n"
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

            # 8.1.1 Alignment for SNP/indel use only
            if variantcalling:
                if not (SV or expansion or MEI or STR or genotypeSTR):
                    print(
                        "\nPerforming paired read alignment with HISAT2...\n"
                        )

                    os.system(
                        "%shisat2 %s --no-spliced-alignment -p %s -x %s -1 %s -2 %s | %s %ssamtools view -@ %s -Sb -  | %ssambamba sort -t %s  --tmpdir=%s -o %ssorted.bam  /dev/stdin"
                        % (path_hisat,hisat_custom_options, num_cpu, path_hisat_index, input_file,
                           input_file2, samblaster_cmq, path_samtools, num_cpu,
                           path_sambamba, num_cpu, tmp_dir, out))

                    bam_file = "%ssorted.bam" % (out)

                    is_variant_file_OK(bam_file, "bam", "alignment")

                    os.system("touch  %slogs/alignment.log" % (out))

                    print("\nCompleted paired read alignment with HISAT2.\n")

            # 8.1.2 Alignment (general)
            # use HISAT2 to align all reads,
            # then soft-clipped and unaligned reads are realigned with BWA mem
                else:
                    if (SV or expansion or MEI or STR or genotypeSTR):
                        if RG:
                            rg_option_hisat2 = " --rg-id %s --rg LB:%s --rg PL:%s  --rg PU:%s --rg SM:%s " % (
                                RG_ID, RG_LB, RG_PL, RG_PU, RG_SM)
                            rg_option_bwa = " -R '@RG\\tID:%s\\tLB:%s\\tPL:%s\\tPU:%s\\tSM:%s' " % (
                                RG_ID, RG_LB, RG_PL, RG_PU, RG_SM)
                        else:
                            rg_option_hisat2 = ""
                            rg_option_bwa = ""

                        print(
                            "\nPerforming paired read alignment with HISAT2...\n"
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
                            "\nPerforming paired read alignment of soft-clipped and unaligned HISAT2 reads with BWA-MEM...\n"
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

            # 8.2.1 Alignment for SNP/indel use only
            if variantcalling:
                if not (SV or expansion or MEI or STR or genotypeSTR):
                    print(
                        "\nPerforming single-end read alignment with HISAT2...\n"
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

            # 8.1.2 Alignment (general)
            # use HISAT2 to align all reads,
            # then soft-clipped and unaligned reads are realigned with BWA mem
                else:
                    if (SV or expansion or MEI or STR or genotypeSTR):
                        if RG:
                            rg_option_hisat2 = " --rg-id %s --rg LB:%s --rg PL:%s  --rg PU:%s --rg SM:%s " % (
                                RG_ID, RG_LB, RG_PL, RG_PU, RG_SM)
                            rg_option_bwa = " -R '@RG\\tID:%s\\tLB:%s\\tPL:%s\\tPU:%s\\tSM:%s' " % (
                                RG_ID, RG_LB, RG_PL, RG_PU, RG_SM)
                        else:
                            rg_option_hisat2 = ""
                            rg_option_bwa = ""

                        print(
                            "\nPerforming single-end read alignment with HISAT2...\n"
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
                            "\nPerforming single-end read alignment of soft-clipped and unaligned HISAT2 reads with BWA-MEM...\n"
                        )
                        os.system(
                            "%sbwa mem %s %s -t %s %s %sunaligned_reads.fq | %s %ssamtools view -Sb -  | %ssambamba sort -t %s --tmpdir=%s -o %ssorted_bwa.bam /dev/stdin; %ssamtools index -@ %s %ssorted_bwa.bam "
                            % (path_bwa, bwa_custom_options, rg_option_bwa, num_cpu, path_bwa_index, out,
                               samblaster_cmq, path_samtools, path_sambamba, num_cpu, tmp_dir,
                               out, path_samtools, num_cpu, out))
                        os.system("%ssamtools view -H %ssorted.bam > %sheader.txt" % (path_samtools, out, out))
                        print("\nMerging HISAT2 and BWA-MEM aligned reads...\n")
                        os.system(
                            "%ssamtools merge -c -@ %s -f -h %sheader.txt %ssorted_merged.bam %ssorted.bam  %ssorted_bwa.bam"
                            % (path_samtools, num_cpu, out, out, out, out))
                        os.system("%ssamtools index -@ %s %ssorted_merged.bam " % (path_samtools, num_cpu, out))

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

            if (SV or expansion or MEI or STR or genotypeSTR):
                bam_file = "%ssorted_merged.bam" % (out)

            if variantcalling:
                if not (SV or expansion or MEI or STR or genotypeSTR):
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

if format == "cram":
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
            #10.1 Strelka snv and indel calling
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
                        "%s %sconfigureStrelkaGermlineWorkflow.py --bam %s --referenceFasta %s --runDir %sstrelka --callRegions %s/sorted.bed.gz"
                        % (path_python2, path_strelka, bam_file, path_reference, out, out))

                        if not debug:
                            os.system("rm %s/temp.bed.gz %s/sorted.bed.gz" % (out, out))

                    if exome:
                        os.system("mkdir %sstrelka" % (out))
                        os.system("%s %sconfigureStrelkaGermlineWorkflow.py --bam %s --referenceFasta %s --runDir %sstrelka --exome" %
                        (path_python2, path_strelka, bam_file, path_reference, out))

                    else:
                        os.system("mkdir %sstrelka" % (out))
                        os.system(
                        "%s %sconfigureStrelkaGermlineWorkflow.py --bam %s --referenceFasta %s --runDir %sstrelka"
                        % (path_python2, path_strelka, bam_file, path_reference, out))

                    os.system("%s %sstrelka/runWorkflow.py -j %s -m local" % (path_python2, out, num_cpu))
                    os.system(
                    "mv %s/strelka/results/variants/genome.S1.vcf.gz  %s/results/%s_strelka.vcf.gz"
                    % (out, out, sample_name))
                    os.system(
                    "mv %s/strelka/results/variants/genome.S1.vcf.gz.tbi  %s/results/%s_strelka.vcf.gz.tbi"
                    % (out, out, sample_name))

                    os.system("mv %s/results/%s_strelka.vcf.gz %s%s_sorted.vcf.gz" % (out, sample_name, out, sample_name))
                    os.system ("%stabix -p vcf %s%s_sorted.vcf.gz" % (path_tabix, out, sample_name))

                    variant_results_file = "%s%s_sorted.vcf.gz" % (out, sample_name)

                    is_variant_file_OK(variant_results_file, "Vcf", "variantcalling")

                    os.system("touch  %slogs/VC_strelka.log" % (out))

                    print('\nSNV and indel calling with Strelka is complete.\n')

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

        print("\nRepeat expansion scanning is complete.\n")

        os.system("touch  %slogs/EH.log" % (out))

        if not debug:
            os.system("rm %stemp_EH.json* %stemp_EH_realigned.bam" % (out, out))

# 13. Compute genome-wide short tandem repeat profiles with ExpansionHunter Denovo
if STR:
    if "STR.log" in os.listdir(out):
        print(
            "WARNING: The presence of EH.log in logs is telling you that the expansion scan was already peformed, please remove SV.log if you wish to perform this stage anyway\n"
        )
    else:
        print("\nExpansionHunter Denovo is scanning the genome to construct a catalog-free short tandem repeat profile...\n")

        os.system("%s/bin/ExpansionHunterDenovo profile --reads %s --reference %s --output-prefix %s/results/%s --min-anchor-mapq 50 --max-irr-mapq 40 --log-reads" % (
            path_expansionHunterDenovo_dir, bam_file, path_reference, out, sample_name))

        STR_profile = "%s/results/%s.str_profile.json" % (out, sample_name)

        if len(STR_profile) != 0:
            print("\nShort tandem repeat profiling is complete.\n")

            if not fast_mode and genotypeSTR == "True":
                #13.1 Convert the novel/non-reference loci identified with ExpansionHunter Denovo to variant catalog format
                print("\nConverting loci identified with ExpansionHunter Denovo into ExpansionHunter variant catalog format...\n")

                os.system("cat %s/results/%s.locus.tsv | sed 's/contig/chr/g' | cut -f1-4 > %s/results/%s_EHDNinput.txt" % (out, sample_name, out, sample_name))

                EHDN_input = "%s/results/%s_EHDNinput.txt" % (out, sample_name)

                EHDN_variant_catalog = "%s/results/%s_EHDN_variant_catalog.json" % (out, sample_name)

                EHDN_excluded = "%s/results/%s_EHDN_excluded.csv" % (out, sample_name)

                EHDN_unmatched = "%s/results/%s_EHDN_unmatched.csv" % (out, sample_name)

                create_variant_catalog.transform_format(EHDN_input, path_reference, EHDN_variant_catalog, EHDN_unmatched, EHDN_excluded)

                #13.2 Genotype the novel/non reference variant catalog with ExpansionHunter
                if len(EHDN_variant_catalog) != 0:
                    print("\nGenotyping denovo loci with ExpansionHunter...\n")

                    os.system("%sExpansionHunter --reads %s --reference %s  --variant-catalog %s --output-prefix %s/temp_EHDN" % (
                    path_expansionHunter, bam_file, path_reference, EHDN_variant_catalog, out))
                    os.system("mv %s/temp_EHDN.vcf %s/results/%s_EHDNexpansions.vcf ; bgzip %s/results/%s_EHDNexpansions.vcf ; %stabix -p vcf %s/results/%s_EHDNexpansions.vcf.gz" % (
                    out, out, sample_name, out, sample_name, path_tabix, out, sample_name))

                    EHDNexpansion_results_file = "%s/results/%s_EHDNexpansions.vcf.gz" % (out, sample_name)

                    is_variant_file_OK(EHDNexpansion_results_file, "Vcf", "EHDNexpansion")

                    print("\nDenovo expansion loci genotyping is complete.\n")

                    if not debug:
                        os.system("rm %stemp_EHDN.json* %stemp_EHDN_realigned.bam, %s, %s, %s" % (out, out, EHDN_excluded, EHDN_unmatched, EHDN_input))

            else:
                print("\nWARNING: You have indicated that you want to run DNAscan2 in fast mode, which does not include the genotyping of STR loci (--genotypeSTR). If you want to run STR genotyping, please remove the fast mode flag.\n")

        else:
            print("\nWARNING: %s is empty, please run again and make sure all paths are correct if you want to perform genome-wide short tandem repeat profiling.\n" % STR_profile)

        os.system("touch  %slogs/STR.log" % (out))

# 14. Structural Variant calling
if SV or MEI:
    if SV:

        if "SV.log" in os.listdir(out + "logs"):

            print(
                "WARNING: The presence of SV.log in logs is telling you that structural variant calling was already peformed, please remove SV.log if you wish to perform this stage anyway\n"
            )

        else:
            #14.1 Manta (all SV)
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
                        "%s %sconfigManta.py --bam %s --referenceFasta %s --runDir %smanta --callRegions %s/sorted.bed.gz"
                        % (path_python2, path_manta, bam_file, path_reference, out, out))

                    if not debug:
                        os.system("rm %s/temp.bed.gz %s/sorted.bed.gz" % (out,out))

                else:
                    print(
                        "\nStructural variants are being called with Manta...\n")

                    os.system("mkdir %smanta" % (out))
                    os.system(
                        "%s %sconfigManta.py --bam %s --referenceFasta %s --runDir %smanta"
                        % (path_python2, path_manta, bam_file, path_reference, out))

                os.system("%s %smanta/runWorkflow.py -j %s -m local" % (path_python2, out, num_cpu))
                os.system("%s %s/convertInversion.py %ssamtools %s %s/manta/results/variants/diploidSV.vcf.gz > %s/results/%s_manta_SV.vcf" % (
                path_python2, path_scripts, path_samtools, path_reference, out, out, sample_name))
                os.system("bgzip -c %s/results/%s_manta_SV.vcf > %s/results/%s_manta_SV.vcf.gz" % (out, sample_name, out, sample_name))
                os.system("%stabix -p vcf %s/results/%s_manta_SV.vcf.gz" % (path_tabix, out, sample_name))

                manta_SV_results_file = "%s/results/%s_manta_SV.vcf" % (out, sample_name)

                if not debug:
                    os.system("rm -r %smanta" % (out))

                print("\nStructural variant calling with Manta is complete.\n")

                # 14.2 Delly (all SV)
                if not fast_mode:
                    os.system("mkdir %sdelly" % (out))

                    if format == "cram":
                        print("\nConverting input cram file to bam format for structural variant detection with Delly...\n")
                        os.system("samtools view -b -h -@ %s -T %s -o %s/%s.bam %s" % (num_cpu, path_reference, out, sample_name, bam_file))
                        os.system("samtools index -@ %s %s/%s.bam" % (num_cpu, out, sample_name))
                        delly_bam = "%s/%s.bam" % (out, sample_name)

                    if format == "fastq":
                        delly_bam = "%ssorted_merged.bam" % (out)

                    if format == "bam":
                        delly_bam = bam_file

                    print("\nStructural variants are being called with Delly...\n")

                    os.system("%sdelly call -g %s -o %sdelly/%s_delly.bcf -x %s %s" % (path_delly, path_reference, out, sample_name, path_delly_exclude_regions, delly_bam))
                    os.system("%sbcftools view %sdelly/%s_delly.bcf > %sdelly/%s_delly_SV.vcf" % (path_bcftools, out, sample_name, out, sample_name))

                    os.system("mv %sdelly/%s_delly_SV.vcf %s/results/" % (out, sample_name, out))
                    os.system("bgzip -c %s/results/%s_delly_SV.vcf > %s/results/%s_delly_SV.vcf.gz" % (out, sample_name, out, sample_name))
                    os.system("%stabix -p vcf %s/results/%s_delly_SV.vcf.gz" % (path_tabix, out, sample_name))

                    delly_SV_results_file = "%s/results/%s_delly_SV.vcf.gz" % (out, sample_name)

                    is_variant_file_OK(delly_SV_results_file, "Vcf", "SV")

                    print("\nStructural variant calling with Delly is complete.\n")

                # 14.3 SV calls with Manta and Delly are merged together using SURVIVOR to create a union set of structural variants.

                    print("\nStructural variants called with Manta and Delly are being merged with SURVIVOR to create a union callset...\n")

                    os.system("ls %s/results/*SV.vcf > %s/results/survivor_sample_files" % (out, out))
                    os.system("%sSURVIVOR merge %s/results/survivor_sample_files 1000 1 1 1 0 30 %s/results/%s_SV_merged.vcf" % (path_SURVIVOR, out, out, sample_name))
                    os.system("perl %svcf-sort.pl %s/results/%s_SV_merged.vcf | bgzip -c > %s/results/%s_SV_merged.vcf.gz" % (path_scripts, out, sample_name, out, sample_name))
                    os.system("%stabix -p vcf %s/results/%s_SV_merged.vcf.gz" % (path_tabix, out, sample_name))

                    SV_results_file = "%s/results/%s_SV_merged.vcf.gz" % (out, sample_name)

                    is_variant_file_OK(SV_results_file, "Vcf", "SV")

                    print("\nMerging of structural variant calls is complete.\n")

                    if not debug:
                        os.system("rm -r %sdelly %s/results/survivor_sample_files %s/results/*.vcf %s/results/*SV.vcf.gz*" % (out, out, out, out))

                        if format == "cram":
                            os.system("rm %s/%s.bam*" % (out, sample_name))

                    os.system("touch %slogs/SV.log" % (out))

                else:
                    SV_results_file = "%s/results/%s_manta_SV.vcf.gz" % (out, sample_name)

                    is_variant_file_OK(SV_results_file, "Vcf", "SV")

                    os.system("touch %slogs/SV.log" % (out))

            else:
                print(
                    "\n\nWARNING: Structural variant calling cannot be performed with single end reads. If you want this to be performed, please provide paired-end aligned reads to DNAscan.\n\n"
                )

    # 15. Transposable element insertion detection with MELT
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

            if reference == "hg38":
                os.system("ls %sme_refs/Hg38/*zip | sed 's/\*//g' > %smelt/transposon.list" % (path_melt, out))

                melt_bed = "%sadd_bed_files/Hg38/Hg38.genes.bed" % (path_melt)

            if exome == "True":
                os.system("%sjava -Xmx%sg -jar %sMELT.jar Single -bamfile %s -h %s -t %smelt/transposon.list -n %s -w %smelt -exome %s" % (
                    path_java, RAM_GB, path_melt, bam_file, path_reference, out, melt_bed, out, melt_custom_options))
            else:
                os.system("%sjava -Xmx%sg -jar %sMELT.jar Single -bamfile %s -h %s -t %smelt/transposon.list -n %s -w %smelt %s" % (
                    path_java, RAM_GB, path_melt, bam_file, path_reference, out, melt_bed, out, melt_custom_options))

            os.system("cat %smelt/SVA.final_comp.vcf | grep '^#' > %smelt/%s.header.txt" % (out, out, sample_name))
            os.system("cat %smelt/SVA.final_comp.vcf | grep -v '^#' > %smelt/%s.sva.vcf" % (out, out, sample_name))
            os.system("cat %smelt/LINE1.final_comp.vcf | grep -v '^#' > %smelt/%s.line1.vcf" % (out, out, sample_name))
            os.system("cat %smelt/ALU.final_comp.vcf | grep -v '^#' > %smelt/%s.alu.vcf" % (out, out, sample_name))
            os.system("cat %smelt/%s.header.txt %smelt/%s.sva.vcf %smelt/%s.line1.vcf %smelt/%s.alu.vcf | perl %svcf-sort.pl -c | bgzip -c > %s/results/%s_MEI.vcf.gz" % (out, sample_name, out, sample_name, out, sample_name, out, sample_name, path_scripts, out, sample_name))
            os.system("%stabix -p vcf %s/results/%s_MEI.vcf.gz" % (path_tabix, out, sample_name))

            MEI_results_file = "%s/results/%s_MEI.vcf.gz" % (out, sample_name)

            is_variant_file_OK(MEI_results_file, "Vcf", "MEI")

            if not debug:
                os.system("rm -r %smelt" % (out))

                if variantcalling:
                    os.system("rm %s/*bam.disc %s/*bam.disc.bai %s/*bam.fq" % (out, out, out))

            os.system("touch  %slogs/mei.log" % (out))

            print("\nTransposable element insertion scanning with MELT is complete.\n")

    if 'SV_results_file' in locals() and os.path.isfile("%s/results/%s_MEI.vcf.gz" % (out, sample_name)) == True:
        print("\nMerging SV and MEI callsets together with SURVIVOR to create a union callset...\n")
        os.system("bgzip -d %s" % (SV_results_file))
        os.system("bgzip -d %s" % (MEI_results_file))
        os.system("ls %s/results/*.vcf > %s/results/survivor_sample_files" % (out, out))
        os.system("%sSURVIVOR merge %s/results/survivor_sample_files 1000 1 1 1 0 30 %s/results/%s_SV_MEI_merged.vcf" % (path_SURVIVOR, out, out, sample_name))
        os.system("perl %svcf-sort.pl %s/results/%s_SV_MEI_merged.vcf | bgzip -c > %s/results/%s_SV_MEI_merged.vcf.gz" % (path_scripts, out, sample_name, out, sample_name))
        os.system("%stabix -p vcf %s/results/%s_SV_MEI_merged.vcf.gz" % (path_tabix, out, sample_name))

        SV_MEI_results_file = "%s/results/%s_SV_MEI_merged.vcf.gz" % (out, sample_name)

        is_variant_file_OK(SV_MEI_results_file, "Vcf", "SVMEI")

        print("\nMerging of SV and MEI variant calls is complete.\n")

        if not debug:
            os.system("rm %s/results/survivor_sample_files %s/results/*.vcf" % (out, out, out, sample_name))

# 16. Annotation with Annovar with optional missense variant prioritisation according to ACMG guidelines (intervar_20180118 database)

if annotation:
    if "annovar.log" in os.listdir(out + "logs"):

        print(
            "WARNING: The presence of annovar.log in logs is telling you that annotation was already peformed, please remove annovar.log if you wish to perform this stage anyway\n"
        )

        variant_results_file = "%s/results/%s_SNPindel_annotated.vcf.gz" % (out, sample_name)

        expansion_variant_results_file = "%sresults/%s_expansions_annotated.vcf.gz" % (out, sample_name)

    else:
        if variantcalling:
            print("\nAnnotation of SNP and indels are being performed using ANNOVAR, with the databases and respective operations defined in paths_configs.py...\n")

            os.system(
                "perl %stable_annovar.pl  --thread %s --vcfinput %s %s -buildver %s -remove -protocol %s -operation %s -nastring . --outfile %s/annovar_SNPindel.vcf"
                % (path_annovar, num_cpu, variant_results_file, path_annovar_db,
                   reference, annovar_protocols, annovar_operations, out))

            if not debug and not alsgenescanner:
                os.system(
                    "rm %sannovar_SNPindel.vcf.%s_multianno.txt %sannovar_SNPindel.vcf.avinput" %
                    (out, reference, out))

            os.system(
                "mv %s/annovar_SNPindel.vcf.%s_multianno.vcf %sresults/%s_SNPindel_annotated.vcf" % (out, reference, out, sample_name))

            os.system("bgzip -f %sresults/%s_SNPindel_annotated.vcf ; %stabix -fp vcf %sresults/%s_SNPindel_annotated.vcf.gz" % (
                out, sample_name, path_tabix, out, sample_name))

            os.system("mv %s %sresults/" % (variant_results_file, out))
            os.system("mv %s.tbi %sresults/" % (variant_results_file, out))

            variant_results_file = "%sresults/%s_SNPindel_annotated.vcf.gz" % (out, sample_name)
            is_variant_file_OK(variant_results_file, "Vcf", "annotation")

            print("\nSNP and indel annotation is complete.\n")

        if expansion:
            if os.path.isfile(EHDNexpansion_results_file) == True:
                print('Repeat expansions and short tandem repeats called with ExpansionHunter software are being merged with SURVIVOR to create a union callset...\n')
                os.system("gzip -dk %s/results/*expansions.vcf.gz" % (out, sample_name))
                os.system("ls %s/results/*expansions.vcf > %s/results/survivor_expansion_sample_files" % (out, out))
                os.system("%sSURVIVOR merge %s/results/survivor_expansion_sample_files 1000 1 1 1 0 30 %s/results/%s_expansions_merged.vcf" % (path_SURVIVOR, out, out, sample_name))
                os.system("perl %svcf-sort.pl %s/results/%s_expansions_merged.vcf | bgzip -c > %s/results/%s_expansions_merged.vcf.gz" % (path_scripts, out, sample_name, out, sample_name))
                os.system("%stabix -p vcf %s/results/%s_expansions_merged.vcf.gz" % (path_tabix, out, sample_name))

                expansion_variant_results_file = "%s/results/%s_expansions_merged.vcf.gz" % (out, sample_name)

                is_variant_file_OK(expansion_variant_results_file, "Vcf", "expansion")

                if not debug:
                    os.system("rm %s/results/survivor_expansion_sample_files %s/results/%s_expansions_merged.vcf %s/results/*expansions.vcf" % (out, out, sample_name, out))

            else:
                expansion_variant_results_file = "%s/results/%s_expansions.vcf.gz" % (out, sample_name)

            print("\nAnnotation of repeat expansions and/or short tandem repeats are being performed using ANNOVAR, with the databases and respective operations defined in paths_configs.py...\n")
            os.system(
                "perl %stable_annovar.pl  --thread %s --vcfinput %s %s -buildver %s -remove -protocol %s -operation %s -nastring . --outfile %s/annovar_expansions.vcf"
                % (path_annovar, num_cpu, expansion_variant_results_file, path_annovar_db,
                   reference, annovar_protocols, annovar_operations, out))

            if not debug:
                os.system(
                    "rm %sannovar_expansions.vcf.%s_multianno.txt %sannovar_expansions.vcf.avinput" %
                    (out, reference, out))

            os.system(
                "mv %s/annovar_expansions.vcf.%s_multianno.vcf %sresults/%s_expansions_annotated.vcf" % (out, reference, out, sample_name))

            os.system("bgzip -f %sresults/%s_expansions_annotated.vcf ; %stabix -fp vcf %sresults/%s_expansions_annotated.vcf.gz" % (
                out, sample_name, path_tabix, out, sample_name))

            os.system("mv %s %sresults/" % (expansion_variant_results_file, out))
            os.system("mv %s.tbi %sresults/" % (expansion_variant_results_file, out))

            expansion_variant_results_file = "%sresults/%s_expansions_annotated.vcf.gz" % (out, sample_name)

            is_variant_file_OK(expansion_variant_results_file, "Vcf", "expansion")

            print("\nRepeat expansion and/or short tandem repeat annotation is complete.\n")

        os.system("touch %slogs/annovar.log" % (out))

        if SV or MEI:
        #16.1 Structural variant annotation and prioritisation is carried out using AnnotSV
            if "annotsv.log" in os.listdir(out + "logs"):
                print(
                    "WARNING: The presence of annotsv.log in logs is telling you that SV annotation was already peformed, please remove annotsv.log if you wish to perform this stage anyway\n"
                )
                if SV:
                    SV_annotation_file = "%s/results/%s_annotated_SV.tsv" % (out, sample_name)
                if MEI:
                    MEI_annotation_file = "%s/results/%s_annotated_MEI.tsv" % (out, sample_name)
                if os.path.isfile("%s/results/%s_SV_MEI_merged.vcf.gz" % (out, sample_name)) == True:
                    SV_MEI_annotation_file = "%s/results/%s_annotated_SV_MEI.tsv" % (out, sample_name)
            else:
                os.environ["ANNOTSV"] = "%s" % (path_annotsv)

                if alsgenescanner or len(path_gene_list) != 0:
                    candidate_gene_cmd = "-candidateGenesFile %s" % (path_gene_list)

                if reference == "hg19":
                    genome_build = "GRCh37"

                else:
                    genome_build = "GRCh38"

                if os.path.isfile("%s/results/%s_SV_MEI_merged.vcf.gz" % (out, sample_name)) == True:
                    print("\nStructural variant and transposable element annotation is being performed with AnnotSV...\n")
                    os.system("%s/bin/AnnotSV -annotationsDir %s/share/AnnotSV/ -bcftools %sbcftools -bedtools %sbedtools -SvinputFile %s %s -genomeBuild %s -outputFile %s/results/%s_annotated_SV_MEI -SVminSize 30 %s" % (
                        path_annotsv, path_annotsv, path_bcftools, path_bedtools, SV_MEI_results_file, candidate_gene_cmd, genome_build, out, sample_name, annotsv_custom_options))

                    SV_MEI_annotation_file = "%s/results/%s_annotated_SV_MEI.tsv" % (out, sample_name)

                    print("\nStructural variant and transposable element annotation is complete.\n")

                else:
                    if SV:

                        print("\nStructural variant annotation is being performed with AnnotSV...\n")
                        os.system("%s/bin/AnnotSV -annotationsDir %s/share/AnnotSV/ -bcftools %sbcftools -bedtools %sbedtools -SvinputFile %s %s -genomeBuild %s -outputFile %s/results/%s_annotated_SV -SVminSize 30 %s" % (
                            path_annotsv, path_annotsv, path_bcftools, path_bedtools, SV_results_file, candidate_gene_cmd, genome_build, out, sample_name, annotsv_custom_options))

                        SV_annotation_file = "%s/results/%s_annotated_SV.tsv" % (out, sample_name)

                        print("\nStructural variant annotation is complete.\n")

                    if MEI:
                        print("\nTransposable element annotation is being performed with AnnotSV...\n")
                        os.system("%s/bin/AnnotSV -annotationsDir %s/share/AnnotSV/ -bcftools %sbcftools -bedtools %sbedtools -SvinputFile %s %s -genomeBuild %s -outputFile %s/results/%s_annotated_MEI -SVminSize 30 %s" % (
                            path_annotsv, path_annotsv, path_bcftools, path_bedtools, MEI_results_file, candidate_gene_cmd, genome_build, out, sample_name, annotsv_custom_options))

                        MEI_annotation_file = "%s/results/%s_annotated_MEI.tsv" % (out, sample_name)

                        print("\nTransposable element annotation is complete.\n")

                os.system("touch  %slogs/annotsv.log" % (out))

        print("\nAnnotation and prioritisation is complete.\n")

else:
    if variant_results_file:
        os.system("mv %s* %sresults/" % (variant_results_file, out))

        variant_results_file = "%s%s_sorted.vcf.gz" % (out, sample_name)


# 17. Microbes screening

if virus or bacteria or custom_microbes:

    if "microbes.log" in os.listdir(out + "logs"):

        print(
            "WARNING: The presence of microbes.log in logs is telling you that microbes scanning was already peformed, please remove microbes.log if you wish to perform this stage anyway\n"
        )
    else:
        # 17.1 Exctract non human reads
        print("\nExtracting non-human reads...\n")
        os.system(
            "%ssamtools view -@ %s -hf 4 %s | %ssamtools bam2fq -s %ssingleton_reads.fastq -@ %s - > %sunaligned_reads.fastq ; cat %ssingleton_reads.fastq >> %sunaligned_reads.fastq ; gzip  %sunaligned_reads.fastq "
            % (path_samtools, num_cpu, bam_file, path_samtools, out, num_cpu,
               out, out, out, out))

        # 17.2 Identifies present non-human microbes
        if virus:

            # 17.2.1 Identifies present viruses
            print("\nIdentifying viruses present in sample...\n")
            os.system(
                "%shisat2 --no-spliced-alignment -p %s -x %s -U %sunaligned_reads.fastq.gz | %ssamtools view -@ %s -hSb -  | %ssamtools sort -@ %s -T %stemp.file -o %soutput_virus.bam -"
                % (path_hisat, num_cpu, path_virus_index, out, path_samtools,
                   num_cpu, path_samtools, num_cpu, out, out))
            os.system(
                "%ssamtools index -@ %s %soutput_virus.bam; %ssamtools idxstats %soutput_virus.bam > %svirus_stats.txt"
                % (path_samtools, num_cpu, out, path_samtools, out, out))

            # 17.2.2 Generates virus report
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

            # 17.2.3 Identifies present bacteria
            print("\nIdentifying bacteria present in sample...\n")

            os.system(
                "%shisat2 --no-spliced-alignment -p %s -x %s -U %sunaligned_reads.fastq.gz | %ssamtools view -hSb -  | %ssamtools sort -T %stemp.file -o %soutput_bacteria.bam -"
                % (path_hisat, num_cpu, path_bacteria_index, out,
                   path_samtools, path_samtools, out, out))
            os.system(
                "%ssamtools index -@ %s %soutput_bacteria.bam; %ssamtools idxstats %soutput_bacteria.bam > %sbacteria_stats.txt"
                % (path_samtools, num_cpu, out, path_samtools, out, out))

            # 17.2.4 Generates bacteria report
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

            # 17.2.5 Identifies present user-selected microbes
            print("\nIdentifying user-selected microbes that are present in sample...\n")

            os.system(
                "%shisat2 --no-spliced-alignment -p %s -x %s -U %sunaligned_reads.fastq.gz | %ssamtools view -hSb -  | %ssamtools sort -T %stemp.file -o %soutput_custom_microbes.bam -"
                % (path_hisat, num_cpu, path_custom_microbes_index, out,
                   path_samtools, path_samtools, out, out))
            os.system(
                "%ssamtools index -@ %s %soutput_custom_microbes.bam; %ssamtools idxstats %soutput_custom_microbes.bam > %scustom_microbes_stats.txt"
                % (path_samtools, num_cpu, out, path_samtools, out, out))

            # 17.2.6 Generates user-selected microbes report
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

# 18. Alignment report generation ( samtools flagstat and stats )

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


# 19. Sequencing data report generation ( fastqc )

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


# 20. Snv and indel calling report generation ( bcftools stats )

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

# 21. Html report generation ( Multiqc )

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

# 22. Annotated variants report generation

if results_report:
    if "annovar.log" not in os.listdir(out + "logs"):

        print(
            "WARNING: Variant annotation was not peformed - please perform variant calling and annotation using the -variantcalling/-expansion and -annotation flags if you wish to generate an ANNOVAR results report.\n"
        )

    else:
        if "results_report.log" in os.listdir(out + "logs"):
            print(
                "WARNING: The presence of results_report.log in logs is telling you that the results report was already produced, please remove results_report.log if you wish to perform this stage anyway\n"
            )

        else:
            print("\nGenerating report of annotated SNP and indel variant calls...\n")

            os.system("zcat %s > %stemp_SNPindel.vcf" % (variant_results_file, out))

            #vcf = open('%stemp_SNPindel.vcf' % (out), 'r')

            #vcf_lines = vcf.readlines()

            gene_list_file = open(path_gene_list)

            gene_list_lines = gene_list_file.readlines()

            gene_list = gene_list_lines

            out_file_all = open(
                '%sreports/%s_annovar_SNPindelvariants.txt' % (out, sample_name), 'w')

            counter = 0

            for i in gene_list:

                with open('%stemp_SNPindel.vcf' % (out)) as vcf:

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

            if not debug:
                os.system("rm %stemp_SNPindel.vcf" % (out))

            if expansion:
                print("\nGenerating report of annotated expansion/short tandem repeat calls...\n")

                os.system("zcat %s > %stemp_expansions.vcf" % (expansion_variant_results_file, out))

                #vcf = open('%stemp_expansions.vcf' % (out), 'r')

                #vcf_lines = vcf.readlines()

                gene_list_file = open(path_gene_list)

                gene_list_lines = gene_list_file.readlines()

                gene_list = gene_list_lines

                out_file_all = open(
                '%sreports/%s_annovar_expansionvariants.txt' % (out, sample_name), 'w')

                counter = 0

                for i in gene_list:

                    with open('%stemp_expansions.vcf' % (out)) as vcf:

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

                if not debug:
                    os.system("rm %stemp_expansions.vcf" % (out))

        #22.1 knotAnnotSV SV report generation
            if SV or MEI:
                    if "annotsv.log" not in os.listdir(out + "logs"):
                        print("WARNING: Structural variant and/or mobile element insertion annotation was not peformed - please perform structural variant calling and/or mobile element insertion calling and annotation using the -SV and/or -MEI and -annotation flags if you wish to generate an AnnotSV results report.\n")

                    else:
                        if os.path.isfile("%s/results/%s_SV_MEI_merged.vcf.gz" % (out, sample_name)) == True:
                            print("\nGenerating structural variant and transposable element annotation HTML report...\n")

                            os.system("mkdir %s%s_SVMEIanno" % (out, sample_name))

                            os.system("perl %sknotAnnotSV.pl --configFile %s/config_AnnotSV.yaml --annotSVfile %s --outDir %s%s_SVMEIanno --genomeBuild %s" % (
                                path_knotannotsv, path_knotannotsv, SV_MEI_annotation_file, out, sample_name, reference))

                            os.system("mv %s%s_SVMEIanno/%s_annotated_SV_MEI.html %s/reports/%s_SV_MEIannotatedvariants.html" % (
                                out, sample_name, sample_name, out, sample_name))

                            if not debug:
                                os.system("rm -r %s%s_SVMEIanno" % (out, sample_name))

                            print("\nStructural element and transposable element HTML report created.\n")
                        else:
                            if SV:
                                print("\nGenerating SV annotation HTML report...\n")

                                os.system("mkdir %s%s_SVanno" % (out, sample_name))

                                os.system("perl %sknotAnnotSV.pl --configFile %s/config_AnnotSV.yaml --annotSVfile %s --outDir %s%s_SVanno --genomeBuild %s" % (
                                    path_knotannotsv, path_knotannotsv, SV_annotation_file, out, sample_name, reference))

                                os.system("mv %s%s_SVanno/%s_annotated_SV.html %s/reports/%s_SVannotatedvariants.html" % (
                                    out, sample_name, sample_name, out, sample_name))

                                if not debug:
                                    os.system("rm -r %s%s_SVanno" % (out, sample_name))

                                print("\nSV HTML report created.\n")

                            if MEI:
                                print("\nGenerating transposable element annotation HTML report...\n")

                                os.system("mkdir %s%s_MEIanno" % (out, sample_name))

                                os.system("perl %sknotAnnotSV.pl --configFile %s/config_AnnotSV.yaml --annotSVfile %s --outDir %s%s_MEIanno --genomeBuild %s" % (
                                    path_knotannotsv, path_knotannotsv, MEI_annotation_file, out, sample_name, reference))

                                os.system("mv %s%s_MEIanno/%s_annotated_MEI.html %s/reports/%s_MEIannotatedvariants.html" % (
                                    out, sample_name, sample_name, out, sample_name))

                                if not debug:
                                    os.system("rm -r %s%s_MEIanno" % (out, sample_name))

                                print("\nTransposable element HTML report created.\n")

            #22.2 Concise results report for all variants called (SV, MEI, expansion, SNVs and indels)
            if os.path.isfile("%s/reports/%s_annovar_SNPindelvariants.txt" % (out, sample_name)) == True:
                if SV_results_file == "%s/results/%s_SV_merged.vcf.gz" % (out, sample_name) or os.path.isfile("%s/results/%s_MEI.vcf.gz" % (out, sample_name)) == True:
                    if os.path.isfile("%s/results/%s_annotated_SV.tsv" % (out, sample_name)) == True and SV:
                        annotsv_file = SV_annotation_file

                    if not SV and os.path.isfile("%s/results/%s_annotated_MEI.tsv" % (out, sample_name)) == True:
                        annotsv_file = MEI_annotation_file

                    print("\nGenerating a concise results report for all annotated variants (SNVs, indels, expansion, SV and/or MEI)...\n")

                    os.system("cat %s | cut -f 2,3,6,15,18,28,29,34,73,87,88 | awk -v OFS='\t' '{split($4,a,/:/);$4=a[1]}1' | awk -v OFS='\t' ' {NR==1?$11=\"Clinvar_ID\t\"$11:$11=\"\t\"$11 } 1 ' | awk  -v OFS='\t' ' {NR==1?$13=\"Clinvar_Phenotype\t\"$13:$13=\"\t\"$13 } 1 ' | awk  -v OFS='\t' ' {NR==1?$14=\"Variant_Frequency_ExAC\t\"$14:$14=\"\t\"$14 } 1 ' | awk  -v OFS='\t' ' {NR==1?$14=\"Variant_Frequency_1000g\t\"$14:$14=\"\t\"$14 } 1 ' | awk  -v OFS='\t' ' {NR==1?$14=\"Variant_Frequency_gnomAD\t\"$14:$14=\"\t\"$14 } 1 ' | awk -F '\t' 'NR>1 {print \"chr\"$0}' > %s/reports/temp_%s_SV_variants.tsv" % (
                        annotsv_file, out, sample_name))

                else:
                    if os.path.isfile("%s/results/%s_annotated_SV_MEI.tsv" % (out, sample_name)) == True:
                        annotsv_file = SV_MEI_annotation_file

                    print("\nGenerating a concise results report for all annotated variants (SNVs, indels, expansion, SV and/or MEI)...\n")

                    os.system("cat %s | cut -f 2,3,6,15,18,28,29,34,73,87,88 | awk -v OFS='\t' '{split($4,a,/:/);$4=a[1]}1' | awk -v OFS='\t' ' {NR==1?$11=\"Clinvar_ID\t\"$11:$11=\"\t\"$11 } 1 ' | awk  -v OFS='\t' ' {NR==1?$12=\"Clinvar_Phenotype\t\"$13:$13=\"\t\"$13 } 1 ' | awk  -v OFS='\t' ' {NR==1?$14=\"Variant_Frequency_ExAC\t\"$14:$14=\"\t\"$14 } 1 ' | awk  -v OFS='\t' ' {NR==1?$14=\"Variant_Frequency_1000g\t\"$14:$14=\"\t\"$14 } 1 ' | awk  -v OFS='\t' ' {NR==1?$14=\"Variant_Frequency_gnomAD\t\"$14:$14=\"\t\"$14 } 1 ' | awk -F '\t' 'NR>1 {print \"chr\"$0}' > %s/reports/temp_%s_SV_variants.tsv" % (
                        annotsv_file, out, sample_name))

                os.system("cat %s/reports/%s_annovar_SNPindelvariants.txt | awk '{print $1 \"\t\" $2 \"\t\" $9 \"\t\" $5 \"\t\" $7 \"\t\" $6 \"\t\" $10 \"\t\" $91 \"\t\" $79 \"\t\" $77 \"\t\" $78 \"\t\" $83 \"\t\" $120 \"\t\" $121}' | awk -v OFS='\t' '{$3=\"small_variant\" ; print ;}' | awk -v OFS='\t' '{split($9,a,/:/);$9=a[5]}1' | awk -v OFS='\t' ' {NR==1?$8=\"Overlapping_Genes\t\"$8:$8=\"\t\"$8 } 1 ' | awk -v OFS='\t' ' {NR==1?$11=\"OMIM_Phenotype\t\"$11:$11=\"\t\"$11 } 1 ' > %s/reports/temp_%s_snvindel_variants.tsv" % (
                    out, sample_name, out, sample_name))

                os.system("cp %sall_variants_report_header.txt %s/reports/%s_all_variants.tsv" % (path_scripts, out, sample_name))
                os.system("tail -n +2 %s/reports/temp_%s_snvindel_variants.tsv >> %s/reports/%s_all_variants.tsv" % (out, sample_name, out, sample_name))
                os.system("tail -n +2 %s/reports/temp_%s_SV_variants.tsv >> %s/reports/%s_all_variants.tsv" % (out, sample_name, out, sample_name))

                if os.path.isfile("%s/reports/%s_annovar_expansionvariants.txt" % (out, sample_name)) == True:
                    os.system("cat %s/reports/%s_annovar_expansionvariants.txt | awk '{print $1 \"\t\" $2 \"\t\" $9 \"\t\" $5 \"\t\" $7 \"\t\" $6 \"\t\" $10 \"\t\" $91 \"\t\" $79 \"\t\" $77 \"\t\" $78 \"\t\" $83 \"\t\" $120 \"\t\" $121}' | awk -v OFS='\t' '{$3=\"STR\" ; print ;}' | awk -v OFS='\t' '{split($9,a,/:/);$9=a[5]}1' | awk -v OFS='\t' ' {NR==1?$8=\"Overlapping_Genes\t\"$8:$8=\"\t\"$8 } 1 ' | awk -v OFS='\t' ' {NR==1?$11=\"OMIM_Phenotype\t\"$11:$11=\"\t\"$11 } 1 ' > %s/reports/temp_%s_expansion_variants.tsv" % (out, sample_name, out, sample_name))

                    os.system("tail -n +2 %s/reports/temp_%s_expansion_variants.tsv >> %s/reports/%s_all_variants.tsv" % (out, sample_name, out, sample_name))

                if not debug:
                    os.system("rm %s/reports/temp*" % (out))

                print("\nConcise results report for all annotated variants (SNVs, indels, expansion, SV and or MEI) is now available.\n")

            os.system("touch %slogs/results_report.log" % (out))

            print("\nResults report for annotated variants is now available.\n")

# 23. Starting iobio services

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
