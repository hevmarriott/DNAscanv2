
from Bio.Seq import Seq
from typing import List, Tuple
from shared_datatypes import Variant, RefGenome, RefSeq, StartPos


def transform_variants(variants: List[Variant], refgen: RefGenome, margin=500, mode='casecontrol') -> Tuple[List[Variant], List[Variant]]:
    transformed_variants, unmatched_variants = [], []

    for v in variants:
        approx_start = v['start'] - margin
        approx_end = v['end'] + margin

        refseq = RefSeq(refgen.fetch(v['contig'], approx_start, approx_end))

        longest_regions = _find_longest_regions(v['motif'], refseq)

        for region in longest_regions:
            motif, positions = region

            if len(positions) == 0:
                unmatched_variants.append(v)

            for p in positions:
                start, end = _update_positions(approx_start, p, len(motif))

                if mode == 'casecontrol':
                    transformed_variants.append(Variant({
                        'contig': v['contig'],
                        'start': start,
                        'end': end,
                        'motif': motif,
                        'p': v['p'],
                        'bonfp': v['bonfp'],
                        'counts': v['counts'],
                    }))
                elif mode == 'outlier':
                    transformed_variants.append(Variant({
                        'contig': v['contig'],
                        'start': start,
                        'end': end,
                        'motif': motif,
                        'top_z': v['top_z'],
                        'high_case_counts': v['high_case_counts'],
                    }))
                elif mode == 'sarah':
                    transformed_variants.append(Variant({
                        'contig': v['contig'],
                        'start': start,
                        'end': end,
                        'motif': motif,
                    }))

    _sort_variants(transformed_variants)
    _sort_variants(unmatched_variants)
    return transformed_variants, unmatched_variants


def _sort_variants(data: List[Variant]):
    data.sort(key=lambda p: p['start'])
    data.sort(key=lambda p: p['contig'])


def _find_longest_regions(motif: str, refseq: RefSeq) -> List[Tuple[str, List[StartPos]]]:
    motif_reversed = str(Seq(motif).reverse_complement())

    positions = _find_longest_perfect_stretch(motif, refseq)
    positions_revmotif = _find_longest_perfect_stretch(motif_reversed, refseq)

    extended_positions = _extend_regions(positions, refseq, motif)
    extended_positions_revmotif = _extend_regions(positions_revmotif, refseq, motif_reversed)

    selected = _select_motif((motif, extended_positions), (motif_reversed, extended_positions_revmotif))

    res = []
    for selection in selected:
        selected_motif, selected_positions = selection
        longest_regions = _select_longest_stretches(selected_positions, selected_motif)
        res.append((selected_motif, longest_regions))

    return res


def _find_longest_perfect_stretch(motif: str, refseq: RefSeq) -> List[StartPos]:
    start_pos, start, end = [], 0, len(refseq)

    index = refseq.find(motif, start, end)

    while index != -1 and start < end:
        if start_pos and _part_of_same_stretch(index, start_pos, len(motif)):
            start_pos[-1]['copies'] = start_pos[-1]['copies'] + 1
        else:
            start_pos.append(StartPos({'index': index, 'copies': 1}))
        index = refseq.find(motif, index+1, end)

    return start_pos


def _part_of_same_stretch(index: int, start_pos: List[StartPos], motifsize: int) -> bool:
    start_pos_diff = index - start_pos[-1]['index']
    size_of_stretch = start_pos[-1]['copies'] * motifsize
    return start_pos_diff == size_of_stretch


def _extend_regions(positions: List[StartPos], refseq: RefSeq, motif: str) -> List[StartPos]:
    for pos in positions:
        pos = _extend_region_left(pos, refseq, motif)
        pos = _extend_region_right(pos, refseq, motif)

    return positions


def _extend_region_left(startpos: StartPos, refseq: RefSeq, motif: str) -> StartPos:
    lseq, lstart = _get_seq_left(refseq, startpos['index'], len(motif))

    if _levenshtein(lseq, motif) == 1:
        llseq, llstart = _get_seq_left(refseq, lstart, len(motif))

        if llseq == motif:
            # region can be extended!
            startpos['index'] = llstart
            startpos['copies'] += 2
            return _extend_region_left(startpos, refseq, motif)

    elif _levenshtein(lseq, motif) == 0:
        startpos['index'] = lstart
        startpos['copies'] += 1
        return _extend_region_left(startpos, refseq, motif)

    return startpos


def _extend_region_right(startpos: StartPos, refseq: RefSeq, motif: str) -> Tuple[int, int]:
    rseq, rstart = _get_seq_right(refseq, startpos['index'], len(motif), startpos['copies'])

    if _levenshtein(rseq, motif) == 1:
        rrseq, _ = _get_seq_right(refseq, rstart, len(motif), 1)
        if rrseq == motif:
            # region can be extended!
            startpos['copies'] += 2
            return _extend_region_right(startpos, refseq, motif)

    elif _levenshtein(rseq, motif) == 0:
        startpos['copies'] += 1
        return _extend_region_right(startpos, refseq, motif)

    return startpos


def _get_seq_left(refseq: RefSeq, start: int, size: str) -> Tuple[str, int]:
    l_end = start
    l_start = l_end - size
    left = refseq[l_start:l_end]
    return left, l_start


def _get_seq_right(refseq: RefSeq, start: int, size: int, copies: int) -> Tuple[str, int]:
    r_start = start + (copies * size)
    r_end = r_start + size
    right = refseq[r_start:r_end]
    return right, r_start


def _levenshtein(s1: str, s2: str) -> int:
    if len(s1) < len(s2):
        return _levenshtein(s2, s1)

    if len(s2) == 0:
        return len(s1)

    previous_row = range(len(s2) + 1)
    for i, c1 in enumerate(s1):
        current_row = [i + 1]
        for j, c2 in enumerate(s2):
            insertions = previous_row[
                             j + 1] + 1  # j+1 instead of j since previous_row and current_row are one character longer
            deletions = current_row[j] + 1  # than s2
            substitutions = previous_row[j] + (c1 != c2)
            current_row.append(min(insertions, deletions, substitutions))
        previous_row = current_row

    return previous_row[-1]


def _select_motif(regular_motif: Tuple[str, List[StartPos]], reversed_motif: Tuple[str, List[StartPos]]) -> List[Tuple[str, List[StartPos]]]:
    motif, positions = regular_motif
    revmotif, positions_revmotif = reversed_motif

    if revmotif == motif:
        return [regular_motif]

    positions.sort(key=lambda p: p['copies'], reverse=True)
    positions_revmotif.sort(key=lambda p: p['copies'], reverse=True)

    if not positions:
        return [reversed_motif]
    if not positions_revmotif:
        return [regular_motif]

    max_copies = positions[0]['copies']
    max_copies_revmotif = positions_revmotif[0]['copies']

    if max_copies > max_copies_revmotif:
        return [regular_motif]
    elif max_copies == max_copies_revmotif:
        return [regular_motif, reversed_motif]
    return [reversed_motif]


def _pick_between_same_maxcopies(regular_motif: Tuple[str, List[StartPos]], reversed_motif: Tuple[str, List[StartPos]]) -> Tuple[str, List[StartPos]]:
    motif, positions = regular_motif
    revmotif, positions_revmotif = reversed_motif

    mindist_motif, pos_motif = _select_pos_closest_to_ehdnpos(positions)
    mindist_revmotif, pos_revmotif = _select_pos_closest_to_ehdnpos(positions_revmotif)

    if mindist_motif <= mindist_revmotif:
        return motif, pos_motif
    return revmotif, positions_revmotif


def _select_pos_closest_to_ehdnpos(positions: List[StartPos], margin=500) -> Tuple[int, List[StartPos]]:
    high_intialisation_value = 100000
    closest_pos = (high_intialisation_value, StartPos(0))
    for p in positions:
        dist = abs(p['index'] - margin)
        if dist < closest_pos[0]:
            closest_pos = (dist, p)
    return closest_pos[0], [closest_pos[1]]


def _select_longest_stretches(positions: List[StartPos], motif: str) -> List[StartPos]:
    longest = []

    if positions:
        positions.sort(key=lambda p: p['copies'], reverse=True)
        max_copies = positions[0]['copies']

        for p in positions:
            if p['copies'] < max_copies:
                break
            elif p not in longest:
                longest.append(p)

        if longest[0]['copies'] == 1 and len(motif) < 10:
            _, longest = _select_pos_closest_to_ehdnpos(longest)

    return longest


def _update_positions(origstart: int, startpos: StartPos, motifsize: int) -> Tuple[int, int]:
    start, copies = startpos['index'], startpos['copies']
    newstart = origstart + start
    newend = newstart + (copies * motifsize)
    return newstart, newend


def inspect(refgen, contig, start, end, motif):
    print(motif)
    print(refgen.fetch(contig, start, end))


def _more_than_5_Ns(refgen: RefGenome, contig: str, start: int, end: int) -> bool:
    seq = refgen.fetch(contig, start, end)
    if seq.count('N') > 5:
        return True
    return False


def filter_out_too_many_Ns(variants: List[Variant], refgen: RefGenome, margin=1000) -> Tuple[List[Variant], List[Variant]]:
    retained, excluded = [], []
    for v in variants:
        if _more_than_5_Ns(refgen, v['contig'], v['start'] - margin, v['end'] + margin):
            excluded.append(v)
        else:
            retained.append(v)
    return retained, excluded
