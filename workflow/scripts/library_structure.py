"""Fixed-length UMI layouts in the physical P5-to-P7 strand orientation."""
import re

P5 = 'AATGATACGGCGACCACCGAGATCTACAC'
P7_FORWARD = 'ATCTCGTATGCCGTCTTCTGCTTG'
POSITIONS = ('before_i5', 'after_i5', 'before_insert', 'after_insert', 'before_i7', 'after_i7')
DEFAULT = dict(mode='auto', adapter_family='auto', p5_sequence=None, p7_sequence=None,
               i5_length=None, i7_length=None, adapter_i5_sequence=None, adapter_i7_sequence=None, umis=[])


def normalize_structure(value=None):
    """Validate CLI/YAML input with the same bounds as the public API."""
    if value is None:
        value = {}
    if not isinstance(value, dict) or set(value) - set(DEFAULT):
        raise ValueError('Unsupported library structure settings')
    result = dict(DEFAULT, **value)
    if result['mode'] not in ('auto', 'custom'):
        raise ValueError('Library structure mode must be auto or custom')
    if result['adapter_family'] not in ('auto', 'truseq', 'nextera'):
        raise ValueError('Choose Auto, ligation or tagmentation inner adapters')
    for name in ('p5_sequence', 'p7_sequence', 'adapter_i5_sequence', 'adapter_i7_sequence'):
        sequence = result[name]
        if sequence is not None and (not isinstance(sequence, str) or not re.fullmatch('[ACGT]{16,80}', sequence)):
            raise ValueError(name + ' must contain 16-80 A/C/G/T bases in P5-to-P7 orientation')
    for name in ('i5_length', 'i7_length'):
        length = result[name]
        if length is not None and (type(length) is not int or not 0 <= length <= 24):
            raise ValueError(name + ' must be Auto or an integer from 0 to 24')
    if not isinstance(result['umis'], list) or len(result['umis']) > len(POSITIONS):
        raise ValueError('Use at most one UMI at each of the six positions')
    seen = set()
    for umi in result['umis']:
        if not isinstance(umi, dict) or set(umi) != {'position', 'length'}:
            raise ValueError('Each UMI requires a position and a fixed length')
        position, length = umi['position'], umi['length']
        if position not in POSITIONS or position in seen:
            raise ValueError('UMI positions must be supported and unique')
        if type(length) is not int or not 1 <= length <= 64:
            raise ValueError('UMI length must be an integer from 1 to 64 bp')
        seen.add(position)
    result['umis'] = sorted((dict(u) for u in result['umis']), key=lambda u: POSITIONS.index(u['position']))
    if result['mode'] == 'auto' and any(result[k] != DEFAULT[k] for k in DEFAULT if k != 'mode'):
        raise ValueError('Choose Customize to set a library structure or add UMIs')
    return result


def umi_lengths(structure):
    return {u['position']: u['length'] for u in structure['umis']}


def handles(structure):
    return structure['p5_sequence'] or P5, structure['p7_sequence'] or P7_FORWARD


def extract_umis(sequence, unit, structure, reverse_complement):
    """Separate anchored outer UMIs from indexes and inner UMIs from the insert.

    Index bases and UMI bases are reported on the physical P5-to-P7 strand.
    Insert coordinates remain on the original FASTQ read, preserving qualities.
    Missing outer anchors leave those UMI values unknown. Fixed UMI extraction
    uses observed bases; it cannot correct indels inside a random UMI sequence.
    """
    lengths = umi_lengths(structure)
    observed = {position: None for position in lengths}
    before, after = lengths.get('before_insert', 0), lengths.get('after_insert', 0)
    if unit['end'] - unit['start'] < before + after:
        return False
    reverse = unit['orientation'] == 'p7_to_p5'
    central = sequence[unit['start']:unit['end']]
    central = reverse_complement(central) if reverse else central
    if before:
        observed['before_insert'] = central[:before]
    if after:
        observed['after_insert'] = central[-after:]
    unit['start'] += after if reverse else before
    unit['end'] -= before if reverse else after
    for index in ('i5', 'i7'):
        raw = unit[index]
        prefix, suffix = lengths.get('before_' + index, 0), lengths.get('after_' + index, 0)
        if raw is None:
            continue
        expected = structure[index + '_length']
        index_size = len(raw) - prefix - suffix
        # Retain the existing two-base ONT index-indel rescue allowance.
        if index_size < 0 or (expected is not None and abs(index_size - expected) > 2):
            unit[index] = unit[index + '_observed'] = None
            continue
        stop = len(raw) - suffix
        unit[index] = raw[prefix:stop]
        unit[index + '_observed'] = reverse_complement(unit[index]) if reverse else unit[index]
        if prefix:
            observed['before_' + index] = raw[:prefix]
        if suffix:
            observed['after_' + index] = raw[stop:]
    unit['umis'] = observed
    unit['umi_status'] = ('complete' if all(v is not None for v in observed.values()) else 'partial') if lengths else 'not_configured'
    unit['umi_sequence'] = '-'.join(observed[p] for p in POSITIONS if p in observed) if lengths and unit['umi_status'] == 'complete' else None
    return True


def adapter_families(structure, templates):
    """Resolve physical custom adapters into the matcher's forward-primer pairs.

    Both inputs follow the displayed P5-to-P7 strand. The matcher stores its
    right primer in the opposite orientation, so reverse-complement only here.
    A blank side retains the selected preset(s); explicit sides override them.
    """
    structure = structure or {}
    left, right = structure.get('adapter_i5_sequence'), structure.get('adapter_i7_sequence')
    reverse = right.translate(str.maketrans('ACGT', 'TGCA'))[::-1] if right else None
    if left and right:
        return {'custom': (left, reverse)}
    selected = structure.get('adapter_family', 'auto')
    return {('custom_' + name if left or right else name): (left or r1, reverse or r2)
            for name, (r1, r2) in templates.items() if selected == 'auto' or name == selected}
