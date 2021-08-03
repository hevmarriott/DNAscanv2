################################################################
# Program: DNAscan - script to run DNAscan on a list of samples
# Version 1.0
# Author: Alfredo Iacoangeli (alfredo.iacoangeli@kcl.ac.uk)
#################################################################

################################################################
# Script structure:
# 1. Import python modules
# 2. Define paths variables
# 3. Define options from command line
# 4. Parse options from command line
# 5. Run DNAscan for each line in the input sample list
#   5.1 Create DNAscan input file option string per line in the input list
#   5.2 Create working dir tree
#   5.3 Run DNAscan for one sample
#   5.4 Create multisample results files
#   5.5 Run ExpansionHunterDenovo to detect unknown/non-reference repeat expansions
################################################################

import argparse, os, paths_configs, loadInNS, os.path

from argparse import RawTextHelpFormatter

# 2. Define paths variables

var = loadInNS.load('scripts/paths_configs.py')
locals().update(var)

# 3. Define options from command line

parser = argparse.ArgumentParser(prog='python3 analyse_list_of_samples.py ', usage='%(prog)s -format "string" -paired "string" -sample_list "string" -out_dir "string" -option_string "string"', description = '############Help Message############ \n\nThis is a script to run DNAscan on a list of samples. Each line of the list must contain the path to one sample. If samples are in paired reads in fastq format and you have two files per sample, these will have to be on the same line spaced bt a tab.\n\n E.g. sample.1.fq.gz  sample.2.fq.gz\n\nDNAscan uses the file paths.py to locate the needed tools and files. Please make sure your paths.py file is properly filled \n\nUsage example: \n\npython analyse_list_of_samples.py -option_string "-format fastq -mode intensive -reference hg19 -alignment -variantcalling -annotation" -out_dir /path/to/dir -sample_list list.txt -format bam\n\nPlease check the following list of required options\n\n################################################', formatter_class=RawTextHelpFormatter)

requiredNamed = parser.add_argument_group('required named arguments')

requiredNamed.add_argument( '-option_string' , required=True , action = "store" , dest = "option_string" , default = "" , help = 'string of option to pass to the main script. Do not include neither -out option nor -in option . E.g. "-format fastq -mode intensive -reference hg19 -alignment -variantcalling -annotation"' )

requiredNamed.add_argument( '-out_dir' , required=True , action = "store" , dest = "out_dir" , default = "" , help = 'working directory [string]' )

requiredNamed.add_argument( '-sample_list' , required=True , action = "store" , dest = "sample_list" , default = "" , help = 'file containing the list of samples to Analyse [string]' )

requiredNamed.add_argument( '-format' , required=True , action = "store" , dest = "format" , help = 'options are bam, sam, fastq, vcf [string]' )

requiredNamed.add_argument( '-paired' , required=True , action = "store" , dest = "paired" , default = "1" , help = 'options are 1 for paired end reads and 0 for single end reads [string]' )


# 4. Parse options from command line

args = parser.parse_args()

option_string = args.option_string

out_dir = args.out_dir

sample_list = args.sample_list

format = args.format

paired = args.paired

list_file = open( "%s" %(sample_list) , 'r' )

list_file_lines = list_file.readlines()

# 5. Run DNAscan for each line in the input sample list

for sample in list_file_lines:

    # 5.1 Create DNAscan input file option string per line in the input list

    if paired == "1" and format == "fastq" :

        input_file_string = "-in %s -in2 %s" %( sample.split('\t')[0] , sample.split('\t')[1].strip() )

        sample_name = sample.split('\t')[0].split("/")[-1].split("1.f")[-2]

    else :

           input_file_string = "-in %s" %( sample.strip() )

           sample_name = sample.split('.')[-2].split("/")[-1]

    # 5.2 Create working dir tree

    os.system( "mkdir %s ; mkdir %s/%s ; touch %s/multisample.txt" %( out_dir , out_dir , sample_name, out_dir ) )

    # 5.3 Run DNAscan for one sample

    os.system( "python3 %s/scripts/DNAscan.py %s -sample_name %s %s -out %s/%s/ " %( dnascan_dir , option_string , sample_name , input_file_string , out_dir , sample_name) )
    
    sample_names = open("%s/multisample_list.txt" % (out_dir), "a")
    sample_names.write('%s\n' % (sample.split('.')[0]))
    sample_names.close()
    
#5.4 Create multisample results files #need to do this so comes with a list of samples, one per line - this is not right

samples = open("%s/multisample_list.txt" % (out_dir) , 'r' )
samples_lines = sample_names.readlines()

if "-variantcalling" in option_string:
    for sample in samples_lines:
        os.system("cp %s/%s/results/%s_sorted_filtered.vcf.gz %s" % (out_dir, sample_name, sample_name, out_dir))

    os.system("%sbcftools merge --merge all --force-samples %s/*filtered.vcf.gz -O v -o %smultisample_filtered.vcf" % (path_bcftools, out_dir, out_dir))

    os.system("perl %svcf-sort.pl  %smultisample_filtered.vcf | bgzip -c > %smultisample_sorted_filtered.vcf.gz"
    % (path_scripts, out_dir, out_dir))

    os.system("tabix -p vcf %smultisample_sorted_filtered.vcf.gz" % (out_dir))

    multisample_results_file = "%smultisample_sorted_filtered.vcf.gz" % (out_dir)

    if "-annotation" in option_string:
        os.system("cp %s/%s/results/%s_annotated.vcf.gz %s" % (out_dir, sample_name, sample_name, out_dir))
        os.system("%sbcftools merge --merge none %s/*annotated.vcf.gz -O v -o %smultisample_annotated.vcf" % (path_bcftools, out_dir, out_dir))
        os.system("perl %svcf-sort.pl  %smultisample_annotated.vcf | bgzip -c > %smultisample_annotated.vcf.gz"
        % (path_scripts, out_dir, out_dir))

        os.system("tabix -p vcf %smultisample_annotated.vcf.gz" % (out_dir))

        multisample_results_file = "%smultisample_annotated.vcf.gz" % (out_dir)

        if "-results_report" in option_string:
                if "annovar.log" not in os.listdir(out_dir + "logs"):
                    print(
                    "WARNING: Annotation was not peformed - please perform annotation using the -annotation flag if you wish to generate an annotation results report.\n"
                    )

                else:
                    if "results_report.log" in os.listdir(out_dir + "logs"):
                        print(
                        "WARNING: The presence of results_report.log in logs is telling you that the results report was already produced, please remove results_report.log if you wish to perform this stage anyway\n"
                        )

                    else:
                        print("\nGenerating multisample report of annotated variant calls...\n")

                        os.system("zcat %s > %stemp.vcf" % (multisample_results_file, out_dir))

                        #vcf = open('%stemp.vcf' % (out), 'r')

                        #vcf_lines = vcf.readlines()

                        gene_list_file = open(path_gene_list)

                        gene_list_lines = gene_list_file.readlines()

                        gene_list = gene_list_lines

                        out_file_all = open('%sreports/multisample_all_variants.txt' % (out_dir), 'w')

                        counter = 0

                        for i in gene_list:

                            with open('%stemp.vcf' % (out_dir)) as vcf:

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
                        
                        print("\nMultisample report of annotated variant calls is now available.\n")


if "-SV" or "-MEI" in option_string:
    if "-SV" in option_string:
        print("\nMerging structural variant files to create multisample SV file...\n")
        for sample in samples_lines:
            if "-mode fast" in option_string:
                os.system("cp %s/%s/results/%s_manta_SV.vcf.gz %s" % (out_dir, sample_name, sample_name, out_dir))
                os.system("gzip -d %s/%s_manta_SV.vcf.gz" % (out_dir, sample_name))
                os.system("ls %s/*SV.vcf > %s/survivor_sample_files" % (out_dir, out_dir))

            os.system("cp %s/%s/results/%s_SV.merged.vcf.gz %s" % (out_dir, sample_name, sample_name, out_dir))
            os.system("gzip -d %s/%s_SV.merged.vcf.gz" % (out_dir, sample_name))
            os.system("ls %s/*SV.merged.vcf > %s/survivor_SV_sample_files" % (out_dir, out_dir))

        os.system("%sSURVIVOR merge %s/survivor_SV_sample_files 1000 1 1 1 0 30 %s/multisample_SV_merged.vcf" % (
        path_SURVIVOR, out_dir, out_dir))
        os.system("bgzip %s/multisample_SV_merged.vcf ; %stabix -p vcf %s/multisample_SV_merged.vcf.gz" % (out_dir, path_tabix, out_dir))

        multisample_SV_results_file = "%s/multisample_SV_merged.vcf.gz" % (out_dir)
        
        print("\nMultisample SV file has been created.\n")

    if "-MEI" in option_string:
        print("\nMerging transposable element files to create multisample MEI file...\n")
        for sample in samples_lines:
            os.system("cp %s/%s/results/%s_MEI.vcf.gz %s" % (out_dir, sample_name, sample_name, out_dir))
            os.system("gzip -d %s/%s_MEI.vcf.gz" % (out_dir, sample_name))

        os.system("ls %s/*MEI.vcf > %s/survivor_MEI_sample_files" % (out_dir, out_dir))
        os.system("%sSURVIVOR merge %s/survivor_MEI_sample_files 1000 1 1 1 0 30 %s/multisample_MEI_merged.vcf" % (
        path_SURVIVOR, out_dir, out_dir))
        os.system("bgzip %s/multisample_MEI_merged.vcf ; %stabix -p vcf %s/multisample_MEI_merged.vcf.gz" % (out_dir, path_tabix, out_dir))

        multisample_MEI_results_file = "%s/multisample_MEI_merged.vcf.gz" % (out_dir)
        
        print("\nMultisample MEI file has been created.\n")

    if "-annotation" in option_string:
        os.environ["ANNOTSV"] = "%s" % (path_annotsv)

        if alsgenescanner or len(path_gene_list) != 0:
            candidate_gene_cmd = "-candidateGenesFile %s" % (path_gene_list)

        if reference == "hg19":
            genome_build = "GRCh37"

        else:
            genome_build = "GRCh38"

        if "-SV" in option_string:
            print("\nMultisample structural variant annotation is being performed with AnnotSV...\n")
            os.system("%s/bin/AnnotSV -annotationsDir %s/share/AnnotSV/ -bcftools %sbcftools -bedtools %sbedtools -SvinputFile %s %s -genomeBuild %s -outputFile %s/multisample_annotated_SV -SVminSize 30 %s" % (
                path_annotsv, path_annotsv, path_bcftools, path_bedtools, multisample_SV_results_file, candidate_gene_cmd, genome_build, out_dir, annotsv_custom_options))

            SV_annotation_file = "%s/multisample_annotated_SV.tsv" % (out_dir)

            print("\nMultisample structural variant annotation is complete.\n")

        if "-MEI" in option_string:
            print("\nMultisample transposable element annotation is being performed with AnnotSV...\n")
            os.system("%s/bin/AnnotSV -annotationsDir %s/share/AnnotSV/ -bcftools %sbcftools -bedtools %sbedtools -SvinputFile %s %s -genomeBuild %s -outputFile %s/multisample_annotated_MEI -SVminSize 30 %s" % (
                path_annotsv, path_annotsv, path_bcftools, path_bedtools, multisample_MEI_results_file, candidate_gene_cmd, genome_build, out_dir, annotsv_custom_options))

            MEI_annotation_file = "%s/multisample_annotated_MEI.tsv" % (out_dir)

            print("\nMultisample transposable element annotation is complete.\n")

        if not debug:
            os.system("rm %s/results/multisample_unannotated*" % (out))

        if "-results_report" in option_string:
            print("\nCreating multisample HTML annotation report...\n")
            if "-SV" in option_string:
                os.system("mkdir %smultisample_SVannoreport" % (out_dir))

                annotation_dir = "%smultisample_SVannoreport" % (out_dir)
                            
                os.system("perl %sknotAnnotSV.pl --configFile %s/config_AnnotSV.yaml --annotSVfile %s --outDir %s --genomeBuild %s" % (
                    path_knotannotsv, path_knotannotsv, SV_annotation_file, annotation_dir, reference))
                            
                os.system("mv %s/multisample_annotated_SV.html %s/reports/multisample_SVannotatedvariants.html" % (
                    annotation_dir, out_dir))

            if "-MEI" in option_string:
                os.system("mkdir %smultisample_MEIannoreport" % (out_dir))

                annotation_dir = "%smultisample_MEIannoreport" % (out_dir)
                            
                os.system("perl %sknotAnnotSV.pl --configFile %s/config_AnnotSV.yaml --annotSVfile %s --outDir %s --genomeBuild %s" % (
                    path_knotannotsv, path_knotannotsv, MEI_annotation_file, annotation_dir, reference))
                
                os.system("mv %s/multisample_annotated_MEI.html %s/reports/multisample_MEIannotatedvariants.html" % (
                    annotation_dir, out_dir))
                            
            print("\nMultisample HTML annotation report is now available.\n")

# 5.5 Run ExpansionHunterDenovo to detect unknown/non-reference repeat expansions
if "-expansion" in option_string:
    print("\nPerforming unknown and non-reference repeat expansion analysis (motif and locus outlier modes) for all analysed samples using ExpansionHunter Denovo...\n")
    manifest_file = open("%s/multisample_manifest.txt" % (out_dir), 'w')

    for sample in samples_lines:
        if "-mode fast" in option_string:
            os.system("%s/bin/ExpansionHunterDenovo profile --reads %s/%s/sorted.bam --reference %s --output-prefix %s/%s --min-anchor-mapq 50 --max-irr-mapq 40" %
            (path_expansionHunterDenovo_dir, out_dir, sample_name, path_reference, out_dir, sample_name))

        os.system("%s/bin/ExpansionHunterDenovo profile --reads %s/%s/sorted_merged.bam --reference %s --output-prefix %s/%s --min-anchor-mapq 50 --max-irr-mapq 40" %
        (path_expansionHunterDenovo_dir, out_dir, sample_name, path_reference, out_dir, sample_name))

        STR_profile = "%s/%s.str_profile.json" % (out_dir, sample_name)

        manifest_file.write("%s\tcase\t%s\n" % (sample_name, STR_profile))

    manifest_file.close()

    os.system("%s/bin/ExpansionHunterDenovo merge --reference %s --manifest %s/multisample_manifest.txt --output-prefix %s/multisample" % (path_expansionHunterDenovo_dir, path_reference, out_dir, out_dir))

    os.system("%s/scripts/outlier.py locus --manifest %s/multisample_manifest.txt--multisample-profile %s/multisample.multisample_profile.json --output %s/multisample.outlier_locus.tsv" % (path_expansionHunterDenovo_dir, out_dir, out_dir, out_dir))

    os.system("%s/scripts/outlier.py motif --manifest %s/multisample_manifest.txt --multisample_profile %s/multisample.multisample_profile.json --output %s/multisample.outlier_motif.tsv" % (path_expansionHunterDenovo_dir, out_dir, out_dir, out_dir))

    print("\nRepeat expansion analysis with ExpansionHunter Denovo is complete\n")

    if "-annotation" in option_string:
        print("\nAnnotating Expansion Hunter Denovo outlier locus results...\n")
        os.system("%s/scripts/annotate_ehdn.sh --ehdn-results %s/multisample.outlier_locus.tsv --ehdn-annotated-results %s/multisample.outlier_locus_annotated.tsv --annovar-annotate-variation %s/annotate_variation.pl --annovar-humandb %s --annovar-buildver %s" %
        (path_expansionHunterDenovo_dir, out_dir, out_dir, path_annovar, path_annovar_db, reference))

        print("\nRepeat expansion annotation is complete.\n")
