from typing import Dict, List
import json
from shared_datatypes import Variant
import csv


def write_variant_catalog(filename: str, data: List[Variant]):
    output = []

    # optional: 'OfftargetRegions' for 'RareRepeat'
    # ReferenceRegion with format 'contignumber:start-end'
    # If just one region, string, otherwise a list
    # If ReferenceRegion a list, VariantType should become a list too of same size, specifying the type of
    # repeat for each reference region

    for d in data:
        region = _format_region(d)
        motif = "(" + (d['motif']) + ")*"
        output.append({
            'LocusId': region,
            'LocusStructure': motif,
            'ReferenceRegion': region,
            'VariantType': 'Repeat',
        })

    with open(filename, 'w+') as f:
        json.dump(output, f)


def _format_region(variant: Variant) -> str:
    return variant['contig'].lstrip('chr') + ":" + str(variant['start']) + "-" + str(variant['end'])


def write_variants_to_csv_format_sarah(filename: str, data: List[Variant]):
    with open(filename, 'w+') as f:
        fields = ['contig', 'start', 'end', 'motif']
        writer = csv.DictWriter(f, fieldnames=fields, dialect="excel", lineterminator='\r')
        writer.writeheader()
        for variant in data:
            writer.writerow(variant)


def write_variants_to_csv_casecontrol_format(filename: str, data: List[Variant]):
    with open(filename, 'w+') as f:
        fields = ['contig', 'start', 'end', 'motif', 'p', 'bonfp', 'counts']
        writer = csv.DictWriter(f, fieldnames=fields, dialect="excel", lineterminator='\r')
        writer.writeheader()
        for variant in data:
            writer.writerow(variant)


def write_variants_to_csv_outlier_format(filename: str, data: List[Variant]):
    with open(filename, 'w+') as f:
        fields = ['contig', 'start', 'end', 'motif', 'top_z', 'high_case_counts']
        writer = csv.DictWriter(f, fieldnames=fields, dialect="excel", lineterminator='\r')
        writer.writeheader()
        for variant in data:
            writer.writerow(variant)


def dump_varcat_to_file(filename: str, data: List[Dict]):
    with open(filename, 'w+') as f:
        json.dump(data, f)
