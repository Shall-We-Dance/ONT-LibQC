"""Transcript coverage figures for the shallow single-cell report."""
import io
from xml.sax.saxutils import escape

from report_palette import POOLED, scope_color
from rna_qc import profile_svg
from scrna_chemistry import scope_label


def support_counts(body):
    genes, fragments = body['contributing_genes'], body['supporting_units']
    return (f'{genes:,} gene' + ('s' if genes != 1 else ''),
            f'{fragments:,} cDNA fragment' + ('s' if fragments != 1 else ''))


def _gene_body_svg(body, title, color):
    """Show the spliced transcript when coding regions are unavailable."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np

    mean = np.asarray(body['mean_relative_coverage'], dtype=float)
    x = (np.arange(len(mean)) + .5) / len(mean)
    with plt.rc_context({'svg.fonttype': 'none', 'font.family': 'DejaVu Sans', 'font.size': 9,
                         'axes.spines.top': False, 'axes.spines.right': False}):
        fig, ax = plt.subplots(figsize=(7.6, 3.1))
        ax.plot(x, mean, color=color, linewidth=1.7)
        if body['contributing_genes'] >= 2:
            sem = np.asarray(body['sem_across_genes'], dtype=float)
            ax.fill_between(x, np.maximum(0, mean - sem), mean + sem, color=color, alpha=.18, linewidth=0)
        ax.axhline(1, color='#b8c0c8', linestyle=':', linewidth=.8)
        ax.grid(axis='y', alpha=.12)
        ax.set(xlim=(0, 1), ylim=(0, max(1.05, ax.get_ylim()[1])), ylabel='Relative coverage',
               xticks=[0, .25, .5, .75, 1], xticklabels=['5′', '25%', '50%', '75%', '3′'],
               xlabel='Position along spliced transcript')
        ax.set_title(title, loc='left', fontsize=10, pad=12)
        fig.text(.1, .025, 'Observed coverage · equal gene weighting · shading: SEM across genes',
                 fontsize=7, color='#526577')
        fig.subplots_adjust(left=.1, right=.98, bottom=.26, top=.85)
        buf = io.StringIO()
        fig.savefig(buf, format='svg', metadata={'Title': title})
        plt.close(fig)
    return buf.getvalue()


def figures(metrics):
    """Return downloadable, self-contained SVGs for scopes with exonic coverage."""
    scopes = metrics.get('scopes', {})
    images = {}
    for scope, profile in scopes.items():
        body = profile['gene_body']
        if not body['contributing_genes']:
            continue
        title = scope_label(scope) + ' · transcript coverage'
        color = {'all': POOLED, 'unassigned': '#738290'}.get(scope, scope_color(scope))
        if profile['CDS']['contributing_genes']:
            svg = profile_svg(profile, title, color)
            axis = '5′ UTR, CDS and 3′ UTR, each scaled to the same width.'
        else:
            svg = _gene_body_svg(body, title, color)
            axis = 'Relative position along spliced transcripts, from 5′ to 3′.'
        gene_count, fragment_count = support_counts(body)
        description = (f'{axis} {gene_count} and {fragment_count} '
                       'in the whole-transcript profile. The line shows mean relative coverage. '
                       'Shading shows one standard error across genes when available.')
        # Matplotlib puts its title in RDF metadata. Add the SVG accessible name
        # and description to the root for people opening the downloaded figure.
        start = svg.index('<svg ')
        end = svg.index('>', start) + 1
        svg = (svg[:start] + svg[start:end - 1] + ' role="img" aria-labelledby="plot-title plot-description">'
               + '<title id="plot-title">' + escape(title) + '</title><desc id="plot-description">'
               + escape(description) + '</desc>' + svg[end:])
        images[scope] = svg.encode('utf-8')
    return images
