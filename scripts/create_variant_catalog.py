import pysam
import transform as t
import read as r
import write as w
from shared_datatypes import RefGenome


# Version for EHdn output in format for Sarah
def transform_format_sarah(ehdn_file: str, refgen_file: str, varcat_out_file: str, unmatched_out_file: str, excluded_out_file: str):
    ehdn_variants = r.read_ehdn_variants_sarah(ehdn_file)
    ref_genome = RefGenome(pysam.FastaFile(refgen_file))

    # If a variant has no basis in the reference genome, it gets excluded.
    transformed_variants, unmatched_variants = t.transform_variants(ehdn_variants, ref_genome, mode='sarah')
    print("Excluded: {}, transformed: {} variants with no basis in the ref.".format(len(unmatched_variants), len(transformed_variants)))

    # If a variant has too many Ns in the surrounding region, it gives an error when running EH
    # The default value is a margin of 1000, corresponding to --region-extension-length arg (=1000) optional argument
    # for EH, to define how far from on/off-target regions to search for informative reads
    retained_variants, excluded_variants = t.filter_out_too_many_Ns(transformed_variants, ref_genome, margin=1000)
    print("Excluded: {}, retained: {} variants after filtering out too many Ns.".format(len(excluded_variants), len(retained_variants)))

    # Output the actual variant catalog to file
    w.write_variant_catalog(varcat_out_file, retained_variants)
    # Output the excluded regions that have no basis in the ref to file
    w.write_variants_to_csv_format_sarah(unmatched_out_file, unmatched_variants)
    # Output the regions excluded due to too many Ns to file
    w.write_variants_to_csv_format_sarah(excluded_out_file, excluded_variants)

    print("Successfully created a variant catalog with {} entries from {} inputted regions.".format(len(retained_variants), len(ehdn_variants)))
    return retained_variants


# Version for EHdn casecontrol output format
def transform_casecontrol_format(ehdn_file: str, refgen_file: str, varcat_out_file: str, unmatched_out_file: str, excluded_out_file: str):
    ehdn_variants = r.read_ehdn_variants_casecontrol(ehdn_file, max_p_value=0.05, min_motifsize=7.0)
    ref_genome = RefGenome(pysam.FastaFile(refgen_file))

    transformed_variants, unmatched_variants = t.transform_variants(ehdn_variants, ref_genome, mode='casecontrol')

    retained_variants, excluded_variants = t.filter_out_too_many_Ns(transformed_variants, ref_genome)
    print("Excluded: {}, retained: {} after filtering out too many Ns".format(len(excluded_variants), len(retained_variants)))

    w.write_variant_catalog(varcat_out_file, retained_variants)
    w.write_variants_to_csv_outlier_format(unmatched_out_file, unmatched_variants)
    w.write_variants_to_csv_outlier_format(excluded_out_file, excluded_variants)
    print("Successfully created a variant catalog of {} variants with {} entries".format(len(ehdn_variants), len(retained_variants)))
    return retained_variants


# Version for EHdn outlier output format
def transform_outlier_format(ehdn_file: str, refgen_file: str, varcat_out_file: str, unmatched_out_file: str, excluded_out_file: str):
    ehdn_variants = r.read_ehdn_variants_outlier(ehdn_file, min_z_value=10.0, min_motifsize=7.0, min_cases_support=0)
    ref_genome = RefGenome(pysam.FastaFile(refgen_file))

    transformed_variants, unmatched_variants = t.transform_variants(ehdn_variants, ref_genome, mode='outlier')
    print("{} variants had no basis in the ref and were excluded.".format(len(unmatched_variants)))

    retained_variants, excluded_variants = t.filter_out_too_many_Ns(transformed_variants, ref_genome, margin=1000)
    print("Excluded: {}, retained: {} after filtering out too many Ns".format(len(excluded_variants), len(retained_variants)))

    w.write_variant_catalog(varcat_out_file, retained_variants)
    w.write_variants_to_csv_outlier_format(unmatched_out_file, unmatched_variants)
    w.write_variants_to_csv_outlier_format(excluded_out_file, excluded_variants)

    print("Successfully created a variant catalog of {} variants with {} entries".format(len(ehdn_variants), len(retained_variants)))
    return retained_variants


# Helper to filter out duplicate entries from variant catalogs when running multiple variant catalogs successively
def retain_unique_entries(varcat_file: str, filter_out_file: str, filtered_varcat_file: str):
    varcat = r.read_in_varcat(varcat_file)
    filter_out = r.read_in_varcat(filter_out_file)
    print("size of orig varcat: {}, size of filter_out_varcat: {}".format(len(varcat), len(filter_out)))
    retain, duplicates = [], []
    for item in varcat:
        if item not in filter_out:
            retain.append(item)
        else:
            duplicates.append(item)

    w.dump_varcat_to_file(filtered_varcat_file, retain)
    print("Filtered out {} entries from the original varcat with {} entries".format(len(duplicates), len(varcat)))


refgenfile = r'example/input/GRCh38_full_analysis_set_plus_decoy_hla.fa'
ehdnfile = r'example/input/Example_InputFile_EH4.txt'
varcat_output = r'example/output/variant_catalog.json'
unmatched_output = r'example/output/unmatched_regions.csv'
excluded_output = r'example/output/automagically_excluded.csv'

transform_format_sarah(ehdnfile, refgenfile, varcat_output, unmatched_output, excluded_output)
