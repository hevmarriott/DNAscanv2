import csv
from typing import List, Dict
import json
from shared_datatypes import Variant


def read_ehdn_variants_sarah(filename: str) -> List[Variant]:
    variant_list = []
    with open(filename) as f:
        data = csv.reader(f, delimiter="\t")
        for d in data:

            contig, start, end, motif = d[0], d[1], d[2], d[3]

            # Excluding the header row
            if contig == 'chr':
                continue

            # Excluding the sex and decoy chromosomes from our analysis
            if any(c in contig for c in ['Un', 'KI', 'EBV', 'random', 'X', 'Y']):
                continue

            try:
                start = int(start)
                end = int(end)
            except ValueError:
                print("Something went wrong converting coordinates to numbers.")

            variant_list.append(Variant({
                'contig': contig,
                'start': start,
                'end': end,
                'motif': motif,
            }))
    return variant_list


def read_ehdn_variants_casecontrol(filename: str, max_p_value=1.0, min_motifsize=0, bonfp_filter=False) -> List[Variant]:
    variant_list = []

    with open(filename) as f:
        data = csv.reader(f, delimiter="\t")
        for d in data:

            contig, start, end, motif, p, bonfp, counts = d[0], d[1], d[2], d[3], d[4], d[5], d[6]

            # Excluding the sex and decoy chromosomes from our analysis
            if any(c in contig for c in ['contig', 'Un', 'KI', 'EBV', 'random', 'X', 'Y']):
                continue

            try:
                start = int(start)
                end = int(end)
                p = float(p)
                bonfp = float(bonfp)
            except ValueError:
                print("Something went wrong converting ehdn data to numbers.")

            if p > max_p_value:
                continue

            if len(motif) < min_motifsize:
                continue

            if bonfp_filter and bonfp == 1:
                continue

            variant_list.append(Variant({
                'contig': contig,
                'start': start,
                'end': end,
                'motif': motif,
                'p': p,
                'bonfp': bonfp,
                'counts': counts,
            }))

    return variant_list


def read_ehdn_variants_outlier(filename: str, min_z_value=1.0, min_motifsize=0, min_cases_support=0) -> List[Variant]:
    variant_list = []
    with open(filename) as f:
        data = csv.reader(f, delimiter="\t")
        for d in data:

            contig, start, end, motif, top_z, high_case_counts, counts = d[0], d[1], d[2], d[3], d[4], d[5], d[6]

            # Excluding the sex and decoy chromosomes from our analysis
            if any(c in contig for c in ['contig', 'Un', 'KI', 'EBV', 'random', 'X', 'Y']):
                continue

            try:
                start = int(start)
                end = int(end)
                top_z = float(top_z)
            except ValueError:
                print("Something went wrong converting ehdn data to numbers.")

            if top_z < min_z_value:
                continue

            if len(motif) < min_motifsize:
                continue

            high_case_counts_list = high_case_counts.split(',')
            if len(high_case_counts_list) < min_cases_support:
                continue

            variant_list.append(Variant({
                'contig': contig,
                'start': start,
                'end': end,
                'motif': motif,
                'top_z': top_z,
                'high_case_counts': high_case_counts_list,
            }))

    return variant_list


def json_to_dict(filename: str) -> List[Dict]:
    res = []

    with open(filename) as f:
        data = json.load(f)
        for d in data:
            contig = d['ehdn']['contig']
            start = d['ehdn']['start']
            end = d['ehdn']['end']
            motif = d['ehdn']['motif']
            if isinstance(motif, list):
                motif = motif[0]

            res.append({
                'contig': contig,
                'start': start,
                'end': end,
                'motif': motif,
            })

    return res


def read_in_varcat(filename: str):
    res = []
    with open(filename) as f:
        data = json.load(f)
        for d in data:
            res.append({
                'LocusId': d['LocusId'],
                'LocusStructure': d['LocusStructure'],
                'ReferenceRegion': d['ReferenceRegion'],
                'VariantType': d['VariantType'],
            })
    return res


def read_gff_file(filename: str):
    with open(filename) as f:
        data = f.readlines()
        for d in data[500:505]:
            print(d)
