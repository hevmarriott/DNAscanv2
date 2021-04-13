If you are reading this file you have already installed all DNAscan dependencies except GATK 4.1.9.0 and Annovar.
Please download GATK-4.1.9.0 at https://github.com/broadinstitute/gatk/releases/download/4.1.9.0/gatk-4.1.9.0.zip and unzip in the desired directory. Annovar requires registration to be used - please register and download Annovar at http://www.openbioinformatics.org/annovar/annovar_download_form.php

After successfully downloading the two softwares you can either change the appropiate paths in paths_configs.py, follow the installation instructions on the software websites and download the needed Annovar databases (see below) or extract annovar.tar download and run the following command:

bash /path/to/DNAscan/scripts/install_annovar_and_gatk.sh $path_to_setup_dir $path_to_DNASCAN_dir $path_to_annovar $path_to_gatk_directory

To download the Annovar databases use the following commands replacing the appropiate paths:

$ANNOVAR_DIR/annotate_variation.pl -buildver hg19 -downdb -webfrom annovar cadd $DNAscan_INSTALL_DIR/humandb/ &

$ANNOVAR_DIR/annotate_variation.pl -buildver hg19 -downdb -webfrom annovar refGene $DNAscan_INSTALL_DIR/humandb/

$ANNOVAR_DIR/annotate_variation.pl -buildver hg19 -downdb -webfrom annovar exac03 $DNAscan_INSTALL_DIR/humandb/

$ANNOVAR_DIR/annotate_variation.pl -buildver hg19 -downdb -webfrom annovar dbnsfp30a $DNAscan_INSTALL_DIR/humandb/

$ANNOVAR_DIR/annotate_variation.pl -buildver hg19 -downdb -webfrom annovar clinvar_20210123 $DNAscan_INSTALL_DIR/humandb/

$ANNOVAR_DIR/annotate_variation.pl -buildver hg19 -downdb -webfrom annovar avsnp147 $DNAscan_INSTALL_DIR/humandb/

IMPORTANT: if you are using Docker to deploy DNAscan, after downloading Annovar and GATK you have to mirror their folders inside the Docker container. You can do this either by using the "docker cp" or the -v flag to mount the volume when starting the container. For more detailed instructions about how to do this please read the $DNAscan_INSTALL_DIR/README.md file.

