#Usage: bash install_dependencies_hg38.sh $path_to_setup_dir $path_to_DNASCAN_dir $path_to_ANNOVAR $path_to_MELT $num_cpu 
#Example: bash install_dependencies_hg38.sh /home/local/ /home/DNAscan /home/annovar.tar.gz /home/MELTv2.2.2.tar.gz 4

INSTALL_DIR=$1

DNASCAN_DIR=$2

ANNOVAR_EXEC=$3

MELT_EXEC=$4

NUM_CPUS=$5

sudo apt-get update

sudo apt-get install -y vim

sudo apt-get install -y python2

sudo apt-get install -y perl

sudo apt-get install -y ttf-dejavu

sudo apt-get install -y wget bzip2 gzip make git tcl http json tar csv

mkdir $INSTALL_DIR

mkdir $INSTALL_DIR/humandb

tar -zxf $ANNOVAR_EXEC --directory $INSTALL_DIR

$ANNOVAR_DIR=$INSTALL_DIR/annovar

cd $DNASCAN_DIR

chmod +x $ANNOVAR_DIR/*

#nohup $ANNOVAR_DIR/annotate_variation.pl -buildver hg38 -downdb -webfrom annovar cadd $INSTALL_DIR/humandb/ &

$ANNOVAR_DIR/annotate_variation.pl -buildver hg38 -downdb -webfrom annovar refGene $INSTALL_DIR/humandb/

$ANNOVAR_DIR/annotate_variation.pl -buildver hg38 -downdb -webfrom annovar exac03 $INSTALL_DIR/humandb/

$ANNOVAR_DIR/annotate_variation.pl -buildver hg38 -downdb -webfrom annovar dbnsfp33a $INSTALL_DIR/humandb/

$ANNOVAR_DIR/annotate_variation.pl -buildver hg38 -downdb -webfrom annovar clinvar_20210501 $INSTALL_DIR/humandb/

$ANNOVAR_DIR/annotate_variation.pl -buildver hg38 -downdb -webfrom annovar avsnp147 $INSTALL_DIR/humandb/

$ANNOVAR_DIR/annotate_variation.pl -buildver hg38 -downdb -webfrom annovar intervar_20180118 $INSTALL_DIR/humandb/

cd $INSTALL_DIR

tar -zxf $MELT_EXEC --directory $INSTALL_DIR

$MELT_DIR=$INSTALL_DIR/MELTv2.2.2

wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

chmod +x Miniconda3-latest-Linux-x86_64.sh

bash Miniconda3-latest-Linux-x86_64.sh -b -p $INSTALL_DIR/Miniconda3/

export PATH=$INSTALL_DIR/Miniconda3/bin:$PATH

echo export PATH=$INSTALL_DIR/Miniconda3/bin:$PATH >> ~/.bashrc

conda config --add channels conda-forge

conda config --add channels defaults

conda config --add channels r

conda config --add channels bioconda

conda install -y python=3.8.0

conda install -y pysam=0.16.0.1

conda install -y samtools=1.9

conda install -y freebayes=1.3.2

conda install -y bedtools=2.25.0

conda install -y vcftools=0.1.16

conda install -y bcftools=1.9

conda install -y hisat2=2.2.1

conda install -y bwa=0.7.17

conda install -y rtg-tools=3.12

conda install -y multiqc=1.10.1

conda install -y fastqc=0.11.9

conda install -y expansionhunter=3.2.2

conda install -y sambamba=0.7.1

conda install -y samblaster=0.1.26

conda install -y delly=0.8.7

conda install -y pysimplegui=4.40.0

conda install -y survivor=1.0.7

conda install -y perl=5.26.2=h470a237_0

conda install -y perl-app-cpanminus

cd $DNASCAN_DIR

mkdir hg38

cd hg38

wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz

gzip -d hg38.fa.gz

samtools faidx hg38.fa

nohup bwa index hg38.fa &

nohup hisat2-build -p $NUM_CPUS hg38.fa hg38 &

sudo apt-get update -qq

sudo apt-get install -y -qq bzip2 gcc g++ make python zlib1g-dev

cd $INSTALL_DIR

wget https://github.com/Illumina/strelka/releases/download/v2.9.10/strelka-2.9.10.centos6_x86_64.tar.bz2

tar -xjf strelka-2.9.10.centos6_x86_64.tar.bz2

export PATH=$INSTALL_DIR/strelka-2.9.10.centos6_x86_64/bin:$PATH

echo export PATH=$INSTALL_DIR/strelka-2.9.10.centos6_x86_64/bin:$PATH >> ~/.bashrc

wget https://github.com/Illumina/manta/releases/download/v1.6.0/manta-1.6.0.centos6_x86_64.tar.bz2

tar -xjf manta-1.6.0.centos6_x86_64.tar.bz2

export PATH=$INSTALL_DIR/manta-1.6.0.centos6_x86_64/bin:$PATH

echo export PATH=$INSTALL_DIR/manta-1.6.0.centos6_x86_64/bin:$PATH >> ~/.bashrc

wget https://github.com/Illumina/ExpansionHunterDenovo/releases/download/v0.9.0/ExpansionHunterDenovo-v0.9.0-linux_x86_64.tar.gz

tar -xzf ExpansionHunterDenovo-v0.9.0-linux_x86_64.tar.gz

export PATH=$INSTALL_DIR/ExpansionHunterDenovo-v0.9.0-linux_x86_64/bin:$PATH

echo export PATH=$INSTALL_DIR/ExpansionHunterDenovo-v0.9.0-linux_x86_64/bin:$PATH >> ~/.bashrc

cd AnnotSV

make PREFIX=. install

make PREFIX=. install-human-annotation

export ANNOTSV=$INSTALL_DIR/AnnotSV/

cd $ANNOTSV/share/AnnotSV/Annotations_Human/Genes/GRCh38

wget https://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz

cd $INSTALL_DIR

git clone https://github.com/mobidic/knotAnnotSV

cpan YAML::XS

cpan Sort::Key::Natural

export PATH=$INSTALL_DIR/knotAnnotSV/:$PATH

echo export PATH=$INSTALL_DIR/knotAnnotSV/:$PATH >> ~/.bashrc

cd $DNASCAN_DIR

sed "s|path_reference = \"\"|path_reference = \"$DNASCAN_DIR\/hg38\/hg38.fa\"|" scripts/paths_configs.py > scripts/paths_configs.py_temp

sed "s|path_hisat_index = \"\"|path_hisat_index = \"$DNASCAN_DIR\/hg38\/hg38\"|" scripts/paths_configs.py_temp > scripts/paths_configs.py

sed "s|path_bwa_index = \"\"|path_bwa_index = \"$DNASCAN_DIR\/hg38\/hg38.fa\"|" scripts/paths_configs.py > scripts/paths_configs.py_temp

sed "s|path_annovar = \"\"|path_annovar = \"$ANNOVAR_DIR\/\"|" scripts/paths_configs.py_temp > scripts/paths_configs.py

sed "s|path_annovar_db = \"\"|path_annovar_db = \"$INSTALL_DIR\/humandb\/\"|" scripts/paths_configs.py > scripts/paths_configs.py_temp

sed "s|dnascan_dir = \"\"|dnascan_dir = \"$DNASCAN_DIR\/\"|" scripts/paths_configs.py_temp > scripts/paths_configs.py

sed "s|path_manta = \"\"|path_manta = \"$INSTALL_DIR\/manta-1.6.0.centos6_x86_64\/bin\/\"|" scripts/paths_configs.py > scripts/paths_configs.py_temp

sed "s|path_expansionHunter_jsons = \"\"|path_expansionHunter_jsons = \"$DNASCAN_DIR\/repeats\/hg38\/\"|" scripts/paths_configs.py_temp > scripts/paths_configs.py

sed "s|path_to_db = \"\"|path_to_db = \"$DNASCAN_DIR\/db\/\"|" scripts/paths_configs.py > scripts/paths_configs.py_temp

sed "s|path_annotsv = \"\"|path_annotsv = \"$INSTALL_DIR\/AnnotSV\/\"|" scripts/paths_configs.py_temp > scripts/paths_configs.py

sed "s|path_knotannotsv = \"\"|path_knotannotsv = \"$INSTALL_DIR\/knotAnnotSV\/\"|" scripts/paths_configs.py > scripts/paths_configs.py_temp

sed "s|path_strelka = \"\"|path_strelka = \"$INSTALL_DIR\/strelka-2.9.10.centos6_x86_64\/bin\/\"|" scripts/paths_configs.py_temp > scripts/paths_configs.py

sed "s|path_delly_exclude_regions = \"\"|path_delly_exclude_regions = \"$DNASCAN_DIR\/db\delly_hg38.excl.tsv\"|"  scripts/paths_configs.py > scripts/paths_configs.py_temp

sed "s|path_expansionHunterDenovo_dir = \"\"|path_expansionHunterDenovo_dir = \"$INSTALL_DIR\/ExpansionHunterDenovo-v0.9.0-linux_x86_64\/\"|" scripts/paths_configs.py_temp > scripts/paths_configs.py

sed "s|path_melt = \"\"|path_melt = \"$INSTALL_DIR\/MELTv2.2.2\/\"|" scripts/paths_configs.py > scripts/paths_configs.py_temp

mv scripts/paths_configs.py_temp scripts/paths_configs.py

chmod +x scripts/*

export PATH=$DNASCAN_DIR/scripts/:$PATH

echo export PATH=$DNASCAN_DIR/scripts/:$PATH >> ~/.bashrc

echo "###########################################IMPORTANT######################################################"
echo "Hisat2-build and bwa-index are still creating their indexes. Please wait untill they complete their task."
echo "You can check whether or not they are still running using the 'top' command"
echo "##########################################################################################################"
