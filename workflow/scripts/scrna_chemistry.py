"""Evidence-based recognition of barcode-bearing ends of physical ONT inserts."""
from collections import Counter
from dataclasses import dataclass
from itertools import product
from typing import Optional

from library_sequences import edit_distance, rc
from scrna_parse_diagnostics import diagnose as diagnose_unmatched_parse

PARSE = 'parse-wt-58'
C4 = 'mgi-c4-linked'
C4_COMPACT = 'mgi-c4-contiguous'
TENX = {
    '10x-3p-v2': ('v2', 10, 'polyT'),
    '10x-3p-v3': ('3p-v3', 12, 'polyT'),
    '10x-3p-v4': ('3p-v4', 12, 'polyT'),
    '10x-5p-v1-v2': ('v2', 10, 'TSO'),
    '10x-5p-v3': ('5p-v3', 12, 'TSO'),
}
LABELS = {
    PARSE: 'Parse Evercode WT · 58-base barcode cassette on Read 2',
    C4: 'MGI DNBelab C4 · linked barcode pair',
    C4_COMPACT: 'MGI DNBelab C4 · contiguous barcode pair',
    **{name: name.replace('10x-', '10x Chromium ').replace('3p-', '3′ ').replace('5p-', '5′ ')
       for name in TENX},
}
TENX_TSO = 'TTTCTTATATGGG'
SMART_HANDLE = 'AAGCAGTGGTATCAACGCAGAGT'
MGI_R1_END = 'GCTCACAGAACGACATGGCTACGATCCGACTT'[-18:]
MGI_R2_START = 'AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA'[:18]
MIN_PARENTS = 10
MIN_COMBINATIONS = 3
PROTOCOL_CATEGORIES = ('parse', '10x', 'mgi', 'multiple', 'none')


def scope_label(scope):
    return {'all': 'Whole library', 'unassigned': 'Unassigned index',
            'native-unindexed': 'Native unindexed'}.get(scope, scope)


def expected_cdna_strand(kit):
    return 'sense' if kit == PARSE or kit.startswith('10x-5p') else 'antisense'


@dataclass
class Unit:
    name: str
    parent: str
    scope: str
    seq: str  # Physical Read 1 -> Read 2, including for raw reverse-complements.
    qual: str
    raw_bases_beyond_r2: Optional[int] = None


@dataclass(frozen=True)
class Tag:
    barcode: str
    umi: str
    quality: str
    cut: int
    edits: int = 0


def poly_t(seq, start):
    end = start
    while end < len(seq) and seq[end] == 'T':
        end += 1
    n = end - start
    return n if n >= 6 else 0


def similar(seq, motif, mismatches=1):
    return len(seq) == len(motif) and sum(a != b for a, b in zip(seq, motif)) <= mismatches


def trim_smart_handle(seq):
    """Accept truncated terminal handles with at least 13 bases of evidence."""
    best, cut = 0, 0
    for template in (SMART_HANDLE, SMART_HANDLE + 'GAATGGG'):
        for offset in range(len(template) - 12):
            motif = template[offset:]
            allowance = 2 if len(motif) >= 20 else 1
            if len(motif) < best:
                break
            # An alignment with <=2 errors contains an intact 4-base seed.
            if not any(motif[i:i + 4] in seq[max(0, i - 2):i + 6] for i in (0, 4, 8)):
                continue
            for n in range(len(motif) - allowance, min(len(seq), len(motif) + allowance) + 1):
                errors = edit_distance(seq[:n], motif, allowance)
                score = len(motif) - 3 * errors
                if errors <= allowance and score > best:
                    best, cut = score, n
    return cut


def round_match(seq, qual, onlist, rescue):
    if seq in onlist:
        return seq, 0
    if not rescue or len(seq) != len(qual) or set(seq) - set('ACGTN'):
        return None
    matches = []
    for i, (base, q) in enumerate(zip(seq, qual)):
        for alt in 'ACGT':
            variant = seq[:i] + alt + seq[i + 1:]
            if base != alt and variant in onlist:
                matches.append((variant, ord(q) - 33))
    if len(matches) == 1 and matches[0][1] <= 20:
        return matches[0][0], 1
    return None


def tag_from_slices(seq, qual, bounds, lists, umi, cut, rescue):
    matches = [round_match(seq[a:b], qual[a:b], refs, rescue)
               for (a, b), refs in zip(bounds, lists)]
    if not all(matches):
        return None
    edits = sum(m[1] for m in matches)
    if edits > 1:
        return None
    a, b = umi
    return Tag('-'.join(m[0] for m in matches), seq[a:b], qual[a:b], cut, edits)


def parse_offsets(length1=12, length2=12):
    second = 18 + length1
    third = second + 8 + length2
    return ((10, 18), (second, second + 8), (third, third + 8))


def learn_parse_linkers(group, lists):
    votes = [Counter(), Counter()]
    parents = set()
    combinations = set()
    for unit in group:
        seq = rc(unit.seq)
        parts = tuple(seq[a:b] for a, b in parse_offsets())
        if all(part in refs for part, refs in zip(parts, lists)):
            combinations.add(parts)
            if unit.parent not in parents:
                votes[0][seq[18:30]] += 1
                votes[1][seq[38:50]] += 1
                parents.add(unit.parent)
    if len(parents) < MIN_PARENTS or len(combinations) < MIN_COMBINATIONS:
        return None
    winners = [counter.most_common(1)[0] for counter in votes]
    if all(count >= 0.60 * len(parents) and len(set(seq)) >= 3 and 'N' not in seq for seq, count in winners):
        return tuple(seq for seq, _ in winners)
    return None


def parse_tag(seq, qual, lists, linkers, rescue):
    if not linkers:
        return None
    options = [12] if not rescue else [12, 11, 13]
    candidates = []
    for length1, length2 in product(options, repeat=2):
        bounds = parse_offsets(length1, length2)
        middle, end = bounds[1], bounds[2]
        left, right = seq[18:middle[0]], seq[middle[1]:end[0]]
        if length1 == length2 == 12 and not rescue:
            if not (similar(left, linkers[0]) and similar(right, linkers[1])):
                continue
            distance = 0
        else:
            distance1 = edit_distance(left, linkers[0], 1)
            distance2 = edit_distance(right, linkers[1], 1)
            if max(distance1, distance2) > 1:
                continue
            distance = distance1 + distance2
        cut = end[1]
        tag = tag_from_slices(seq, qual, bounds, lists, (0, 10), cut + poly_t(seq, cut), rescue)
        if tag:
            candidates.append((distance + tag.edits, tag))
    if not candidates:
        return None
    best = min(c[0] for c in candidates)
    tags = {tag for score, tag in candidates if score == best}
    return tags.pop() if len(tags) == 1 else None


def barcode_exact(kit, seq, refs):
    if kit == PARSE:
        return all(seq[a:b] in ref for (a, b), ref in zip(parse_offsets(), refs['parse']))
    if kit in TENX:
        return seq[:16] in refs['tenx'][TENX[kit][0]]
    bounds = ((6, 16), (22, 32)) if kit == C4 else ((0, 10), (10, 20))
    return all(seq[a:b] in ref for (a, b), ref in zip(bounds, refs['mgi']))


def tag_for(kit, seq, qual, refs, linkers=None, rescue=False):
    if kit == PARSE:
        return parse_tag(seq, qual, refs['parse'], linkers, rescue)
    if kit in TENX:
        book, umi_len, anchor = TENX[kit]
        if seq[:16] not in refs['tenx'][book]:
            return None
        end = 16 + umi_len
        if anchor == 'polyT':
            tail = poly_t(seq, end)
        else:
            tail = len(TENX_TSO) if similar(seq[end:end + len(TENX_TSO)], TENX_TSO) else 0
        return Tag(seq[:16], seq[16:end], qual[16:end], end + tail) if tail else None
    if kit == C4:
        if not all(similar(seq[i:i + len(s)], s) for i, s in ((0, 'TCTGCG'), (16, 'CCTTCC'), (32, 'CGATG'))):
            return None
        bounds, start = ((6, 16), (22, 32)), 37
    else:
        bounds, start = ((0, 10), (10, 20)), 20
    end = start + 10
    tail = poly_t(seq, end)
    if not tail:
        return None
    return tag_from_slices(seq, qual, bounds, refs['mgi'], (start, end), end + tail, rescue)


def detect(group, refs):
    linkers = learn_parse_linkers(group, refs['parse'])
    parents = {u.parent for u in group}
    evidence = {}
    unit_protocols = ['none'] * len(group)
    for kit in LABELS:
        family = kit.partition('-')[0]
        support, combos = set(), set()
        exact = supported = wrong_end = 0
        for index, unit in enumerate(group):
            seq, qual = (rc(unit.seq), unit.qual[::-1]) if kit == PARSE else (unit.seq, unit.qual)
            opposite = unit.seq if kit == PARSE else rc(unit.seq[-58:])
            wrong_end += barcode_exact(kit, opposite, refs)
            if not barcode_exact(kit, seq, refs):
                continue
            exact += 1
            tag = tag_for(kit, seq, qual, refs, linkers)
            if tag:
                supported += 1
                support.add(unit.parent)
                combos.add(tag.barcode)
                # One unit can match sibling kit versions or distinct protocol families.
                previous = unit_protocols[index]
                if previous == 'none':
                    unit_protocols[index] = family
                elif previous != family:
                    unit_protocols[index] = 'multiple'
        evidence[kit] = {'exact_barcode_units': exact, 'opposite_end_exact': wrong_end, 'supported_units': supported,
                         'parents': len(support), 'distinct_barcodes': len(combos)}
    supported = sorted((kit for kit, e in evidence.items() if e['parents'] >= MIN_PARENTS and
                        e['distinct_barcodes'] >= MIN_COMBINATIONS and e['parents'] >= 0.05 * len(parents)
                        and e['exact_barcode_units'] >= 3 * e['opposite_end_exact']),
                       key=lambda k: -evidence[k]['parents'])
    if len(supported) > 1:
        status, kit = 'mixed_or_ambiguous', None
    elif supported:
        status, kit = 'identified', supported[0]
    else:
        status, kit = ('insufficient_evidence' if len(parents) < MIN_PARENTS else 'unresolved'), None
    counts = Counter(unit_protocols)
    unmatched = (unit for unit, family in zip(group, unit_protocols) if family == 'none')
    return {'status': status, 'kit': kit, 'units': len(group), 'parents': len(parents),
            'evidence': evidence, 'protocol_support': {key: counts[key] for key in PROTOCOL_CATEGORIES},
            'unmatched_parse': diagnose_unmatched_parse(unmatched, refs['parse'], linkers)}, linkers


def extract(unit, kit, refs, linkers):
    barcode_seq, barcode_qual = (rc(unit.seq), unit.qual[::-1]) if kit == PARSE else (unit.seq, unit.qual)
    tag = tag_for(kit, barcode_seq, barcode_qual, refs, linkers, rescue=True)
    if not tag:
        return None
    if kit == PARSE:
        # Read 1 carries the transcript. Remove a terminal template-switch handle
        # only when its sequence is present on that end.
        start = trim_smart_handle(unit.seq)
        end = len(unit.seq) - tag.cut
        seq, qual = unit.seq[start:end], unit.qual[start:end]
    else:
        seq, qual = unit.seq[tag.cut:], unit.qual[tag.cut:]
    return tag, seq, qual


def native_mgi_units(parent, sequence, quality):
    """Only accept complete, ordered MGI primer pairs in unmatched raw reads."""
    found = []
    for reverse, (seq, qual) in enumerate(((sequence, quality), (rc(sequence), quality[::-1]))):
        starts = []
        pos = seq.find(MGI_R1_END)
        while pos != -1:
            starts.append(pos + len(MGI_R1_END))
            pos = seq.find(MGI_R1_END, pos + len(MGI_R1_END))
        for start, bound in zip(starts, starts[1:] + [len(seq)]):
            end = seq.find(MGI_R2_START, start, bound)
            if end < 0 or seq.find(MGI_R2_START, end + len(MGI_R2_START), bound) >= 0:
                continue
            if not 50 <= end - start <= 5000:
                continue
            a, b = (len(seq) - end, len(seq) - start) if reverse else (start, end)
            found.append((a, b, Unit(f'{parent}|mgi:{a}-{b}', parent, 'native-unindexed',
                                     seq[start:end], qual[start:end], len(seq) - end)))
    return [unit for a, b, unit in found if all(other is unit or e <= a or s >= b for s, e, other in found)]
