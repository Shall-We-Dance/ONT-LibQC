"""Read-2 cassette evidence in units left unmatched by every strict protocol test."""
from collections import Counter
from dataclasses import dataclass

from library_sequences import rc


CATEGORIES = ('barcode', 'linker', 'shifted', 'short_end', 'mixed', 'unresolved')
SITES = ('bc3', 'linker1', 'bc2', 'linker2', 'bc1')
STATES = ('exact', 'near', 'mismatch', 'incomplete')


def substitutions(reference):
    for motif in reference:
        for index, old in enumerate(motif):
            for base in 'ACGTN':
                if base != old:
                    yield motif[:index] + base + motif[index + 1:]


def unique_neighbors(reference):
    return {word for word, n in Counter(substitutions(reference)).items() if n == 1}


def linker_words(motif):
    lengths = {12: {word: 1 for word in substitutions([motif])},
               11: {motif[:i] + motif[i + 1:]: 1 for i in range(12)},
               13: {motif[:i] + b + motif[i:]: 1 for i in range(13) for b in 'ACGTN'}}
    lengths[12][motif] = 0
    return lengths


@dataclass(frozen=True)
class Hit:
    start: int
    end: int
    edits: int


def linker_hits(seq, words, last_start):
    # Search only around the reconstructed Read 2 boundary. Allow prefix loss
    # and up to 12 extra bases before the expected position.
    hits = []
    for length, variants in words.items():
        for start in range(min(last_start + 1, len(seq) - length + 1)):
            edits = variants.get(seq[start:start + length])
            if edits is not None:
                hits.append(Hit(start, start + length, edits))
    return hits


class Cassette:
    def __init__(self, barcodes, linkers):
        self.barcodes = barcodes
        self.words = [linker_words(motif) for motif in linkers]
        self.sites = [(10, 18, barcodes[0], unique_neighbors(barcodes[0])),
                      (18, 30, {linkers[0]}, self.words[0][12]),
                      (30, 38, barcodes[1], unique_neighbors(barcodes[1])),
                      (38, 50, {linkers[1]}, self.words[1][12]),
                      (50, 58, barcodes[2], unique_neighbors(barcodes[2]))]

    def canonical(self, seq):
        states = []
        for start, end, onlist, near in self.sites:
            word = seq[start:end]
            if end > len(seq):
                state = 'incomplete'
            elif word in onlist:
                state = 'exact'
            elif word in near:
                state = 'near'
            else:
                state = 'mismatch'
            states.append(state)
        return states

    def classify(self, seq, states):
        bc = states[::2]
        exact = bc.count('exact')
        linkers = sum(state in ('exact', 'near') for state in states[1::2])
        if len(seq) < 58 and (linkers == 2 or (linkers and exact)):
            return 'short_end', 'short_insert'
        if linkers == 2 and len(seq) >= 58 and exact < 3:
            return 'barcode', None
        if exact == 3 and linkers < 2:
            return 'linker', None

        lefts = linker_hits(seq, self.words[0], 30)
        rights = linker_hits(seq, self.words[1], 50)
        # Prefer directly on-list neighbors and the smallest edit explanation.
        # A single short motif found somewhere near the boundary is insufficient.
        candidates = []
        for left in lefts:
            for right in rights:
                gap = right.start - left.end
                if not 7 <= gap <= 9:
                    continue
                hits = sum(word in book for word, book in zip(
                    (seq[max(0, left.start - 8):left.start], seq[left.end:right.start],
                     seq[right.end:right.end + 8]), self.barcodes))
                score = (hits, -(left.edits + right.edits), -abs(gap - 8),
                         -abs(left.start - 18), -abs(right.start - 38))
                if hits and left.start <= 15 and right.start <= 35:
                    reason = ('short_end', 'early_anchors')
                elif hits and right.end + 8 > len(seq) and len(seq) < 58:
                    reason = ('short_end', 'short_insert')
                elif hits == 3 and left.start == 18 and gap == 8:
                    reason = ('linker', None)
                elif hits >= 2 and left.start == 18 and left.end == 30 and gap != 8 and right.end - right.start == 12:
                    reason = ('barcode', None)
                elif hits and (left.start != 18 or right.start != 38 or gap != 8):
                    reason = ('shifted', None)
                else:
                    reason = ('mixed', None)
                candidates.append((score, reason))
        if candidates:
            score = max(s for s, _ in candidates)
            reasons = {reason for s, reason in candidates if s == score}
            return reasons.pop() if len(reasons) == 1 else ('mixed', None)

        for right in rights:
            # Where the first linker would cross the Read 2 boundary, require
            # both barcode neighbors flanking linker 2 to be directly on-list.
            if (right.start < 20 and seq[max(0, right.start - 8):right.start] in self.barcodes[1]
                    and seq[right.end:right.end + 8] in self.barcodes[2]):
                return 'short_end', 'early_anchors'
        if exact >= 2 or (linkers and exact):
            return 'mixed', None
        for hits, before, after in ((lefts, self.barcodes[0], self.barcodes[1]),
                                    (rights, self.barcodes[1], self.barcodes[2])):
            for hit in hits:
                if (seq[max(0, hit.start - 8):hit.start] in before
                        or seq[hit.end:hit.end + 8] in after):
                    return 'mixed', None
        return 'unresolved', None


def diagnose(units, barcodes, linkers):
    """Count physical units, retaining only aggregate evidence and the opposite-end check."""
    result = {'units': 0, 'assessed': 0, 'partial': 0, 'shorter_than_58': 0, 'opposite_end_partial': 0,
              'categories': dict.fromkeys(CATEGORIES, 0),
              'short_end_patterns': {'short_insert': 0, 'early_anchors': 0},
              'raw_r2_boundary': {'measured': 0, 'at_least_1': 0, 'at_least_12': 0, 'min_bp': None},
              'sites': {name: dict.fromkeys(STATES, 0) for name in SITES}}
    model = Cassette(barcodes, linkers) if linkers else None
    boundary = result['raw_r2_boundary']
    for unit in units:
        result['units'] += 1
        result['shorter_than_58'] += len(unit.seq) < 58
        distance = unit.raw_bases_beyond_r2
        if distance is not None:
            boundary['measured'] += 1
            boundary['at_least_1'] += distance >= 1
            boundary['at_least_12'] += distance >= 12
            if boundary['min_bp'] is None or distance < boundary['min_bp']:
                boundary['min_bp'] = distance
        if model is None:
            continue
        result['assessed'] += 1
        seq = rc(unit.seq[-72:])
        states = model.canonical(seq)
        category, short_pattern = model.classify(seq, states)
        result['categories'][category] += 1
        result['partial'] += category != 'unresolved'
        if short_pattern:
            result['short_end_patterns'][short_pattern] += 1
        for name, state in zip(SITES, states):
            result['sites'][name][state] += 1
        opposite = unit.seq[:72]
        category, _ = model.classify(opposite, model.canonical(opposite))
        result['opposite_end_partial'] += category != 'unresolved'
    return result
