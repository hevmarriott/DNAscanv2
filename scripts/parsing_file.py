import argparse
from argparse import RawTextHelpFormatter

def create_parser():
  parser = argparse.ArgumentParser(
      prog='python DNAscan.py',
      usage=
      '%(prog)s [options] -format "string" -reference "string" -in "string" \n please use the --help flag to print further information\n ',
      description=
      '############DNAscan help Message############ \n\nDNAscan uses the file paths_configs.py to locate the needed tools and files. Please make sure your paths_configs.py file is properly filled \n\nUsage example: \n\npython DNAscan.py -format fastq -out /home/user/test/ -in sample_sorted_1.fq.gz -reference hg19 -alignment -variantcalling -results_report\n\nPlease check the following list of optional and required options\n\n################################################',
      formatter_class=RawTextHelpFormatter)

  parser.add_argument(
      '-RG',
      action="store_true",
      dest="RG",
      default=False,
      help=
      'if this flag is set the alignment stage will use the read group provided in paths_configs.py (Default = "False")'
  )

  parser.add_argument(
      '-tmp',
      action="store_true",
      dest="tmp",
      default=False,
      help=
      'if this flag is set the alignment stage will use the directory provided in paths_configs.py to save the temporary files instead of /tmp  (Default = "False")'
  )


  parser.add_argument(
      '-alsgenescanner',
      action="store_true",
      dest="alsgenescanner",
      default=False,
      help=
      'if this flag is set DNAscan will perform the ALSgeneScanner analyses and reports (Default = "False")'
  )

  parser.add_argument(
      '-filter_string',
      action="store",
      dest="filter_string",
      default='FORMAT/FT == "PASS" && FORMAT/DP > 10 && MQ > 40 && GQ > 20 && ID/SB < 2 && ADF > 0 && ADR > 0',
      help='bcftools filter string for strelka, eg "GQ>20 & DP>10" (Default = \'FORMAT/FT == "PASS" && FORMAT/DP > 10 && MQ > 40 && GQ > 20 && ID/SB < 2 && ADF > 0 && ADR > 0\')')

  parser.add_argument(
      '-paired',
      action="store",
      dest="paired",
      default="1",
      help=
      'options are 1 for paired end reads and 0 for single end reads (Default = "1")'
  )

  parser.add_argument(
      '-vcf', action="store", dest="vcf", help='complementary vcf file')

  parser.add_argument(
      '-in2',
      action="store",
      dest="input_file2",
      help='input file 2, for paired end reads only (fastq file)')

  parser.add_argument(
      '-iobio',
      action="store_true",
      dest="iobio",
      help=
      'if this flag is set the iobio service links will be provided at the end of the analysis (Default = "False")',
      default=False)

  parser.add_argument(
      '-fast_mode',
      action="store_true",
      dest="fast_mode",
      help=
      'if this flag is set DNAscan2 will run without SV calling with Delly and the genotyping of STR loci identified with ExpansionHunter Denovo if those flags are set (Default = "False")',
      default=False)

  parser.add_argument(
      '-alignment',
      action="store_true",
      dest="alignment",
      help=
      'if this flag is set the alignment stage will be performed (Default = "False")',
      default=False)

  parser.add_argument(
      '-expansion',
      action="store_true",
      dest="expansion",
      help=
      'if this flag is set DNAscan will look for the expansions described in the json folder in paths_configs.py with ExpansionHunter and non-reference/novel tandem repeats with ExpansionHunter Denovo  (Default = "False") ',
      default=False)

  parser.add_argument(
      '-STR',
      action="store_true",
      dest="STR",
      help=
      'if this flag is set DNAscan will look for genome-wide novel/non-reference tandem repeats with ExpansionHunter Denovo  (Default = "False") ',
      default=False)

  parser.add_argument(
      '-genotypeSTR',
      action="store_true",
      dest="genotypeSTR",
      help=
      'if this flag is set in addition to STR, DNAscan will genotype the identified STR loci with ExpansionHunter (Default = "False") ',
      default=False)

  parser.add_argument(
      '-SV',
      action="store_true",
      dest="SV",
      help=
      'if this flag is set the structural variant calling stage will be performed with Manta (in all modes) and Delly (normal mode for inversions and deletions; all variants in intensive mode) (Default = "False") ',
      default=False)

  parser.add_argument(
      '-MEI',
      action="store_true",
      dest="MEI",
      help=
      ' if this flag is set, mobile element insertion/transposable element calling will be performed with MELT (Default = "False") ',
      default=False)

  parser.add_argument(
      '-BED',
      action="store_true",
      dest="BED",
      help=
      'restrict the analysis to the regions in the bed file (Default = "False") ',
      default=False)

  parser.add_argument(
      '-virus',
      action="store_true",
      dest="virus",
      help=
      'if this flag is set DNAscan will perform viral scanning (Default = "False")  ',
      default=False)

  parser.add_argument(
      '-bacteria',
      action="store_true",
      dest="bacteria",
      help=
      'if this flag is set DNAscan will perform bacteria scanning (Default = "False") ',
      default=False)

  parser.add_argument(
      '-custom_microbes',
      action="store_true",
      dest="custom_microbes",
      help=
      'if this flag is set DNAscan will perform a customized microbe scanning according to the provided microbe data base in paths_configs.py (Default = "False")  ',
      default=False)

  parser.add_argument(
      '-variantcalling',
      action="store_true",
      dest="variantcalling",
      help=
      'if this flag is set DNAscan will perform snv and indel calling (Default = "False")  ',
      default=False)

  parser.add_argument(
      '-annotation',
      action="store_true",
      dest="annotation",
      help=
      'if this flag is set DNAscan will annotate the found SNV and indel variants with Annovar and AnnotSV for SV and MEI variants (Default = "False")  ',
      default=False)

  parser.add_argument(
      '-results_report',
      action="store_true",
      dest="results_report",
      help=
      'if this flag is set DNAscan will generate a results report describing the annotated variants (Default = "False")  ',
      default=False)

  parser.add_argument(
      '-alignment_report',
      action="store_true",
      dest="alignment_report",
      help=
      'if this flag is set DNAscan will generate an alignment report (Default = "False") ',
      default=False)

  parser.add_argument(
      '-sequencing_report',
      action="store_true",
      dest="sequencing_report",
      help=
      'if this flag is set DNAscan will generate a report describing the input sequencing data (Default = "False") ',
      default=False)

  parser.add_argument(
      '-calls_report',
      action="store_true",
      dest="calls_report",
      help=
      'if this flag is set DNAscan will generate a report describing the found snvs and small indels (Default = "False")',
      default=False)

  parser.add_argument(
      '-sample_name',
      action="store",
      dest="sample_name",
      default="sample",
      help='specify sample name [string] (default = "sample")')

  parser.add_argument(
      '-rm_dup',
      action="store",
      dest="rm_dup",
      help=
      'do you want to mark duplicates while aligning the reads to the reference genome? Options are True or False [string] (Default = "True")',
      default="True")

  parser.add_argument(
      '-debug',
      action="store_true",
      dest="debug",
      help=
      'if this flag is set DNAscan will not delete intermediete and temporary files (Default = "False")',
      default=False)

  parser.add_argument(
      '-exome',
      action="store_true",
      dest="exome",
      help=
      'if this flag is set DNAscan will only look for variants in the whole exome (Default = "False")',
      default=False)

  parser.add_argument(
      '-format',
      action="store",
      dest="format",
      default="fastq",
      help='options are bam, sam, cram, fastq, vcf [string] ')

  parser.add_argument(
      '-reference',
      action="store",
      dest="reference",
      default="hg19",
      help=
      'options are sus_scrofa_99, hg19, hg38, grch37 and grch38 the path to the reference fasta file must be specified in paths_configs.py [string]'
  )

  parser.add_argument(
      '-ref_file',
      action="store",
      dest="ref_file",
      default= "",
      help=
      'path to the reference file in fasta format [string]'
  )

  parser.add_argument(
      '-dnascan_dir',
      action="store",
      dest="dnascan_main_dir",
      default= "",
      help=
      'path to the DNAscan main dir [string]'
  )



  return parser
