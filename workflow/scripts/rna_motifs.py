"""Bounded HOMER motif discovery on oriented, technically trimmed RNA inserts."""
import hashlib
import io
import math
import os
from pathlib import Path
import random
import re
import shutil
import subprocess
import tempfile


def parse_motifs(text, target_count, background_count):
    result, current = [], None
    for line in text.splitlines():
        if line.startswith('>'):
            fields = line[1:].split('\t')
            if len(fields) < 6:
                raise ValueError('Incomplete HOMER motif header')
            counts = re.search(r'T:([\d.]+)\([\d.]+%\),B:([\d.]+)\([\d.]+%\)', fields[5])
            if not counts:
                raise ValueError('Missing HOMER motif support')
            observed, background = map(float, counts.groups())
            logp = float(fields[3])
            if not all(math.isfinite(x) for x in (observed, background, logp)) or not 0 <= observed <= target_count or not 0 <= background <= background_count:
                raise ValueError('Invalid HOMER counts')
            expected = background * target_count / background_count
            current = {'motif': fields[0].replace('T', 'U'), 'dna_kmer': fields[0],
                       'observed_units': observed, 'expected_units': expected,
                       'target_fraction': observed / target_count, 'background_fraction': background / background_count,
                       'enrichment': (observed + .5) / (expected + .5),
                       'homer_log10_p': logp / math.log(10), 'pwm': []}
            result.append(current)
        elif line.strip():
            values = [float(x) for x in line.split()]
            if current is None or len(values) != 4 or not all(math.isfinite(x) and x >= 0 for x in values) or sum(values) <= 0:
                raise ValueError('Invalid HOMER motif matrix')
            current['pwm'].append([x / sum(values) for x in values])
    if any(not row['pwm'] or len(row['pwm']) != len(row['dna_kmer']) for row in result):
        raise ValueError('Incomplete HOMER motif matrix')
    return result


def adapter_match(consensus):
    """Flag exact adapter-sequence matches without labelling them biological motifs."""
    from ont_qc import FAMILIES, P5, P7
    word = consensus.replace('U', 'T')
    if len(word) < 6 or set(word) - set('ACGT'):
        return False
    adapters = [P5, P7, *(seq for pair in FAMILIES.values() for seq in pair)]
    return any(word in seq or word in seq.translate(str.maketrans('ACGT', 'TGCA'))[::-1] for seq in adapters)


def discover_homer(sequences, shuffle, minimum=100, seed=1729, shuffles=8, maximum=2000, executable=None, timeout=90):
    distinct = sorted({seq.upper().replace('U', 'T') for _, seq in sequences
                       if len(seq) >= 8 and set(seq.upper()) <= set('ACGTU')}, key=lambda s: hashlib.sha256(s.encode()).digest())
    selected, bases = [], 0
    for seq in distinct:
        if len(selected) == maximum:
            break
        if bases + len(seq) <= 2_000_000:
            selected.append(seq)
            bases += len(seq)
    output = {'engine': 'HOMER', 'status': 'insufficient_units', 'eligible_units': len(sequences),
              'distinct_sequences': len(distinct), 'sampled_units': len(selected), 'minimum_units': minimum,
              'shuffle_replicates': shuffles, 'widths': [5, 6, 8], 'threads': 1, 'maximum_input_bases': 2_000_000,
              'background': 'Dinucleotide-preserving shuffled RNA inserts; identical sequences collapsed before discovery',
              'interpretation': 'HOMER de novo motifs on aligned query sequence in annotated RNA strand; terminal soft clips excluded; enrichment against shuffled sequence, not IP/input specificity. HOMER scores are discovery statistics, not motif FDR. Adapter matches flag exact consensus matches of at least six bases to a standard library adapter on either strand',
              'motifs': []}
    if len(selected) < minimum:
        return output
    tool = shutil.which(executable or os.environ.get('HOMER2', 'homer2'))
    if not tool:
        output.update(status='engine_unavailable')
        return output
    output['executable_sha256'] = hashlib.sha256(Path(tool).read_bytes()).hexdigest()
    rng = random.Random(seed)
    rows = []
    try:
        # FASTA and native HOMER files stay outside the portable report and are removed afterwards.
        with tempfile.TemporaryDirectory(prefix='ont-rna-motifs-') as tmp:
            folder = Path(tmp)
            foreground, background = folder/'target.fa', folder/'background.fa'
            with foreground.open('w') as target, background.open('w') as control:
                for i, seq in enumerate(selected):
                    target.write(f'>t{i}\n{seq}\n')
                    for j in range(shuffles):
                        control.write(f'>b{i}_{j}\n{shuffle(seq, rng)}\n')
            for width in output['widths']:
                motifs = folder/f'{width}.motifs'
                command = [tool, 'denovo', '-i', str(foreground), '-b', str(background), '-o', str(motifs),
                           '-len', str(width), '-S', '3', '-strand', '+', '-p', '1']
                subprocess.run(command, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, timeout=timeout)
                rows.extend(parse_motifs(motifs.read_text(), len(selected), len(selected)*shuffles))
    except subprocess.TimeoutExpired:
        output.update(status='timed_out')
        return output
    except (OSError, subprocess.CalledProcessError, ValueError):
        output.update(status='engine_failed')
        return output
    motifs, seen = [], set()
    for row in sorted(rows, key=lambda x: (x['homer_log10_p'], -x['enrichment'])):
        if row['observed_units'] < 5 or row['target_fraction'] <= row['background_fraction'] or row['motif'] in seen:
            continue
        row['adapter_sequence_match'] = adapter_match(row['motif'])
        seen.add(row['motif']); motifs.append(row)
    output.update(status='complete', motifs=motifs[:5])
    return output


def logo_svg(pwm):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.textpath import TextPath
    from matplotlib.font_manager import FontProperties
    from matplotlib.patches import PathPatch
    from matplotlib.transforms import Affine2D
    colors = dict(A='#348a65', C='#3276b5', G='#c49a35', U='#c45c5c')
    with plt.rc_context({'svg.fonttype': 'none'}):
        fig, ax = plt.subplots(figsize=(2.3, .65))
        for index, column in enumerate(pwm):
            information = max(0., 2 + sum(p * math.log2(p) for p in column if p))
            bottom = 0
            for p, base in sorted(zip(column, 'ACGU')):
                height = p * information
                if height < .005:
                    continue
                glyph = TextPath((0, 0), base, prop=FontProperties(family='DejaVu Sans', weight='bold'), size=1)
                box = glyph.get_extents(); sx, sy = .85 / box.width, height / box.height
                transform = Affine2D().scale(sx, sy).translate(index-box.x0*sx, bottom-box.y0*sy)
                ax.add_patch(PathPatch(glyph, transform=transform+ax.transData, color=colors[base], linewidth=0))
                bottom += height
        ax.set(xlim=(-.1,len(pwm)), ylim=(0,2)); ax.axis('off'); fig.subplots_adjust(left=0,right=1,bottom=0,top=1)
        output=io.StringIO();fig.savefig(output,format='svg',transparent=True);plt.close(fig)
    svg=output.getvalue();return svg[svg.index('<svg'):].replace('<svg ', '<svg style="width:100%;height:auto;display:block" ', 1)
