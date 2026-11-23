"""Transcript-region composition of uniquely aligned cDNA bases."""
from collections import Counter
from dataclasses import replace
import io

import rna_qc as rna

REGIONS = (('five_prime_UTR', '5′ UTR', '#4787a7'),
           ('CDS', 'CDS', '#226e63'),
           ('three_prime_UTR', '3′ UTR', '#b3693e'))
OTHERS = (('ambiguous_UTR', 'Overlapping 5′ and 3′ UTR'),
          ('other_exon', 'Other exon'), ('intron', 'Intron'),
          ('genic_unannotated', 'Other genic'), ('intergenic', 'Intergenic'))


def transcript_features(gene):
    """Resolve each transcript on its own strand before merging isoforms."""
    cds, five, three = [], [], []
    for tx in gene.transcripts.values():
        exons = rna.merge_intervals(tx['exon'])
        segments = rna.transcript_segments(tx, exons, gene.strand)
        cds.extend(segments['CDS'])
        five.extend(segments['five_prime_UTR'])
        three.extend(segments['three_prime_UTR'])
    return replace(gene, cds=rna.merge_intervals(cds), utr5=rna.merge_intervals(five), utr3=rna.merge_intervals(three))


def summarize(rows, alignments, annotation):
    genes = {}
    counts, units, scopes = Counter(), Counter(), {}
    for row in rows:
        if row['mapping_status'] != 'unique':
            continue
        primary, = (r for r in alignments[row['unit_id']] if not (r['supplementary'] or r['secondary']))
        blocks = primary['blocks']
        overlaps = []
        for gene in annotation.query(primary['chrom'], blocks):
            if gene.gene_id not in genes:
                genes[gene.gene_id] = transcript_features(gene)
            overlaps.append(genes[gene.gene_id])
        assigned = rna.region_bases(blocks, overlaps)
        if assigned['five_prime_UTR']:
            left = [x for g in overlaps for x in g.utr5]
            right = [x for g in overlaps for x in g.utr3]
            conflict = rna.intersect_intervals(left, right)
            if conflict:
                conflict = rna.subtract_intervals(conflict, [x for g in overlaps for x in g.cds])
                assigned['ambiguous_UTR'] = rna.overlap_length(conflict, blocks)
                assigned['five_prime_UTR'] -= assigned['ambiguous_UTR']
        scope = row['scope_id'] or 'unassigned'
        scopes.setdefault(scope, Counter()).update(assigned)
        units[scope] += 1
        counts.update(assigned)

    def totals(counter, n):
        return {'unique_units': n, 'aligned_bases': sum(counter.values()),
                'coding_region_bases': sum(counter[k] for k, _, _ in REGIONS),
                'counts': {key: counter[key] for key in (*rna.REGIONS, 'ambiguous_UTR')}}
    return {**totals(counts, sum(units.values())), 'scopes': {k: totals(v, units[k]) for k, v in sorted(scopes.items())},
            'annotation_sha256': annotation.metadata['sha256']}


def figure(metrics):
    from matplotlib.backends.backend_svg import FigureCanvasSVG
    from matplotlib.figure import Figure

    den = metrics['coding_region_bases']
    if not den:
        return None
    fig = Figure(figsize=(7.6, 1.65))
    ax = fig.add_axes([.022, .08, .956, .87])
    ax.set(xlim=(0, 1), ylim=(0, 1))
    ax.axis('off')
    left = 0
    for i, (key, label, color) in enumerate(REGIONS):
        n = metrics['counts'][key]
        fraction = n / den
        ax.barh(.79, fraction, left=left, height=.21, color=color)
        left += fraction
        x = i / 3
        ax.text(x, .50, label, color=color, fontsize=12, weight='bold', va='top')
        ax.text(x, .25, f'{fraction:.1%}   ·   {n:,} bp', color='#263b50', fontsize=11, va='top')
    buf = io.BytesIO()
    FigureCanvasSVG(fig).print_svg(buf, metadata={'Title': 'Coding transcript regions', 'Date': None})
    return buf.getvalue()
